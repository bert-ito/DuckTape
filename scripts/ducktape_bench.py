#!/usr/bin/env python3
"""
Resource monitoring for the DuckTape stage scripts (Linux only; reads /proc).

Subcommands
-----------
run     Run one command, wait for it and append its resource usage to a TSV
        file. CPU time, the largest process, page faults, context switches
        and I/O come from the kernel and include every child process the
        command started. The memory of the whole process tree is sampled
        while the command runs.
sample  Sample the memory of a process tree (a whole pipeline step) until
        stopped with SIGTERM, then write the peaks to a file.
step    Append the resource usage of one pipeline step to a TSV file.
total   Append the stage total to a TSV file and print a summary.

Memory columns
--------------
peak_mem_mb  Peak memory of the whole process tree: the sampled sum of PSS
             (proportional set size, which splits shared pages between the
             processes that share them), and never less than max_rss_mb.
peak_rss_mb  Sampled peak of the summed RSS. Pages shared between processes
             (for example forked R workers) are counted once per process, so
             this is an upper bound.
max_rss_mb   RSS of the largest single process, from the kernel. Exact, and
             not affected by the sampling interval.
"""

import argparse
import fcntl
import os
import shlex
import signal
import sys
import threading
import time

CLK_TCK = os.sysconf("SC_CLK_TCK")
MONITOR_MARKER = os.path.basename(__file__).encode()
HAVE_PSS = os.path.exists("/proc/self/smaps_rollup")
HAVE_CHILDREN = os.path.exists(f"/proc/self/task/{os.getpid()}/children")

COMMAND_COLUMNS = [
    "stage", "step", "label", "tool", "start", "end",
    "wall_s", "cpu_user_s", "cpu_sys_s", "cpu_s", "avg_cores",
    "peak_mem_mb", "peak_rss_mb", "max_rss_mb",
    "read_mb", "write_mb", "disk_read_mb", "disk_write_mb",
    "major_faults", "vol_ctx_switches", "invol_ctx_switches",
    "max_procs", "max_threads", "exit_status", "command",
]

STEP_COLUMNS = [
    "stage", "step", "start", "end",
    "wall_s", "cpu_user_s", "cpu_sys_s", "cpu_s", "avg_cores",
    "peak_mem_mb", "peak_rss_mb", "max_rss_mb",
    "read_mb", "write_mb", "disk_read_mb", "disk_write_mb",
    "max_procs", "max_threads", "n_commands",
]

IO_KEYS = ("rchar", "wchar", "read_bytes", "write_bytes")


# ---------------------------------------------------------------------------
# /proc readers
# ---------------------------------------------------------------------------

def _read(path):
    with open(path, "rb") as handle:
        return handle.read()


def _children(pid):
    """Direct children of pid, from /proc/<pid>/task/<tid>/children."""
    kids = []
    try:
        tids = os.listdir(f"/proc/{pid}/task")
    except OSError:
        return kids
    for tid in tids:
        try:
            kids.extend(int(k) for k in _read(f"/proc/{pid}/task/{tid}/children").split())
        except OSError:
            pass
    return kids


def _ppid_table():
    """Map of parent pid -> child pids, from a scan of /proc/*/stat."""
    table = {}
    for name in os.listdir("/proc"):
        if not name.isdigit():
            continue
        try:
            stat = _read(f"/proc/{name}/stat")
        except OSError:
            continue
        # Fields after "pid (comm) ": state, ppid, ...
        fields = stat[stat.rfind(b")") + 2:].split()
        table.setdefault(int(fields[1]), []).append(int(name))
    return table


def process_tree(root):
    """Return root and all of its descendants."""
    if HAVE_CHILDREN:
        children_of = _children
    else:
        table = _ppid_table()
        children_of = lambda pid: table.get(pid, ())  # noqa: E731
    tree, seen, todo = [], set(), [root]
    while todo:
        pid = todo.pop()
        if pid in seen:
            continue
        seen.add(pid)
        tree.append(pid)
        todo.extend(children_of(pid))
    return tree


def process_memory(pid):
    """Return (rss_kb, pss_kb, threads) for one process, or None if it is gone."""
    try:
        status = _read(f"/proc/{pid}/status")
    except OSError:
        return None
    rss = threads = 0
    for line in status.splitlines():
        if line.startswith(b"VmRSS:"):
            rss = int(line.split()[1])
        elif line.startswith(b"Threads:"):
            threads = int(line.split()[1])
    pss = rss
    if HAVE_PSS and rss:
        try:
            for line in _read(f"/proc/{pid}/smaps_rollup").splitlines():
                if line.startswith(b"Pss:"):
                    pss = int(line.split()[1])
                    break
        except (OSError, ValueError):
            pass
    return rss, pss, threads


def is_monitor(pid):
    """True for ducktape_bench.py processes, so steps do not count their memory."""
    try:
        return MONITOR_MARKER in _read(f"/proc/{pid}/cmdline")
    except OSError:
        return False


def read_io(pid):
    """Cumulative I/O counters of a process, including waited-for children."""
    try:
        data = _read(f"/proc/{pid}/io")
    except OSError:
        return None
    counters = {}
    for line in data.splitlines():
        key, _, value = line.partition(b":")
        counters[key.decode()] = int(value)
    return [counters.get(key) for key in IO_KEYS]


class TreeSampler:
    """Track the peak memory, process and thread counts of a process tree."""

    def __init__(self, root, skip_monitors=False):
        self.root = root
        self.skip_monitors = skip_monitors
        self.peak_rss = 0       # kB, summed over the tree
        self.peak_pss = 0       # kB, summed over the tree
        self.max_procs = 0
        self.max_threads = 0
        self.samples = 0

    def sample(self):
        rss = pss = threads = procs = 0
        for pid in process_tree(self.root):
            if self.skip_monitors and pid != self.root and is_monitor(pid):
                continue
            mem = process_memory(pid)
            if mem is None:
                continue
            rss += mem[0]
            pss += mem[1]
            threads += mem[2]
            procs += 1
        self.peak_rss = max(self.peak_rss, rss)
        self.peak_pss = max(self.peak_pss, pss)
        self.max_procs = max(self.max_procs, procs)
        self.max_threads = max(self.max_threads, threads)
        self.samples += 1


# ---------------------------------------------------------------------------
# Formatting and TSV output
# ---------------------------------------------------------------------------

def _num(value, digits=2):
    return "NA" if value is None else f"{value:.{digits}f}"


def _mb_from_kb(kb):
    return None if kb is None else kb / 1024


def _mb_from_bytes(nbytes):
    return None if nbytes is None else nbytes / 1048576


def _iso(epoch):
    return time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(epoch))


def _epoch(text):
    # bash's EPOCHREALTIME uses the locale's decimal separator
    return float(text.replace(",", "."))


def _clean(value):
    return str(value).replace("\t", " ").replace("\n", " ")


def append_row(path, columns, values):
    """Append one row; the header is written first if the file is empty."""
    line = "\t".join(_clean(v) for v in values) + "\n"
    try:
        with open(path, "a") as handle:
            fcntl.flock(handle, fcntl.LOCK_EX)
            if os.fstat(handle.fileno()).st_size == 0:
                handle.write("\t".join(columns) + "\n")
            handle.write(line)
    except OSError as exc:
        print(f"ducktape_bench.py: cannot write {path}: {exc}", file=sys.stderr)


def read_rows(path, stage=None, step=None):
    try:
        with open(path) as handle:
            lines = handle.read().splitlines()
    except OSError:
        return []
    if not lines:
        return []
    header = lines[0].split("\t")
    rows = [dict(zip(header, line.split("\t"))) for line in lines[1:] if line]
    return [
        r for r in rows
        if (stage is None or r.get("stage") == stage)
        and (step is None or r.get("step") == step)
    ]


def _max_field(rows, key):
    values = [float(r[key]) for r in rows if r.get(key, "NA") != "NA"]
    return max(values) if values else None


def _parse_counters(text, n):
    """Parse 'a,b,...' (trailing comma allowed) into n ints, or None."""
    parts = [p for p in (text or "").split(",") if p != ""]
    if len(parts) != n:
        return None
    return [int(p) for p in parts]


def _io_delta(before, after):
    if before is None or after is None:
        return [None] * len(IO_KEYS)
    return [
        None if a is None or b is None else _mb_from_bytes(a - b)
        for b, a in zip(before, after)
    ]


# ---------------------------------------------------------------------------
# Subcommands
# ---------------------------------------------------------------------------

def cmd_run(args):
    command = args.command
    if command and command[0] == "--":
        command = command[1:]
    if not command:
        sys.exit("ducktape_bench.py run: no command given")

    io_before = read_io("self")
    start = time.time()
    t0 = time.monotonic()
    try:
        # Python ignores SIGPIPE and SIGXFSZ; restore the defaults for the
        # child, as bash would, so pipelines behave normally.
        pid = os.posix_spawnp(
            command[0], command, os.environ,
            setsigdef=(signal.SIGPIPE, signal.SIGXFSZ),
        )
    except OSError as exc:
        print(f"{command[0]}: {exc.strerror}", file=sys.stderr)
        return 127

    def forward(signum, frame):
        try:
            os.kill(pid, signum)
        except ProcessLookupError:
            pass

    for sig in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
        if signal.getsignal(sig) != signal.SIG_IGN:
            signal.signal(sig, forward)

    sampler = TreeSampler(pid)
    stop = threading.Event()

    def sample_loop():
        while True:
            try:
                sampler.sample()
            except Exception:  # sampling must never break the pipeline
                pass
            if stop.wait(args.interval):
                break

    thread = threading.Thread(target=sample_loop, daemon=True)
    thread.start()

    # Wait for the child to exit without reaping it, so the sampler cannot
    # look at a recycled pid; then reap it to collect its resource usage.
    os.waitid(os.P_PID, pid, os.WEXITED | os.WNOWAIT)
    stop.set()
    thread.join()
    _, status, usage = os.wait4(pid, 0)
    wall = time.monotonic() - t0
    end = time.time()
    io_after = read_io("self")

    code = os.waitstatus_to_exitcode(status)
    exit_status = f"signal {-code}" if code < 0 else str(code)

    cpu = usage.ru_utime + usage.ru_stime
    max_rss = usage.ru_maxrss                     # kB on Linux
    peak_rss = max(sampler.peak_rss, max_rss)
    peak_mem = max(sampler.peak_pss, max_rss)
    tool = args.label.split()[0] if args.label.split() else command[0]

    append_row(args.log, COMMAND_COLUMNS, [
        args.stage, args.step, args.label, tool, _iso(start), _iso(end),
        _num(wall), _num(usage.ru_utime), _num(usage.ru_stime), _num(cpu),
        _num(cpu / wall if wall > 0 else None),
        _num(_mb_from_kb(peak_mem), 1), _num(_mb_from_kb(peak_rss), 1),
        _num(_mb_from_kb(max_rss), 1),
        *(_num(v, 1) for v in _io_delta(io_before, io_after)),
        usage.ru_majflt, usage.ru_nvcsw, usage.ru_nivcsw,
        sampler.max_procs, sampler.max_threads, exit_status,
        shlex.join(command),
    ])
    return 128 - code if code < 0 else code


def cmd_sample(args):
    sampler = TreeSampler(args.pid, skip_monitors=True)

    class Stop(Exception):
        pass

    def stop(signum, frame):
        raise Stop

    signal.signal(signal.SIGTERM, stop)
    try:
        # Stop if the pipeline shell has gone away (we would be reparented).
        while os.getppid() == args.pid:
            sampler.sample()
            time.sleep(args.interval)
    except Stop:
        pass
    signal.signal(signal.SIGTERM, signal.SIG_IGN)

    with open(args.out, "w") as handle:
        handle.write(
            f"{sampler.peak_rss}\t{sampler.peak_pss}\t"
            f"{sampler.max_procs}\t{sampler.max_threads}\t{sampler.samples}\n"
        )
    return 0


def _read_peak_file(path):
    try:
        fields = _read(path).split()
        return [int(f) for f in fields[:4]]
    except (OSError, ValueError):
        return [0, 0, 0, 0]


def _cpu_seconds(ticks):
    return None if ticks is None else [t / CLK_TCK for t in ticks]


def cmd_step(args):
    start, end = _epoch(args.start), _epoch(args.end)
    wall = end - start
    cpu0 = _parse_counters(args.cpu0, 2)
    cpu1 = _parse_counters(args.cpu1, 2)
    user = sys_ = cpu = None
    if cpu0 and cpu1:
        user, sys_ = _cpu_seconds([b - a for a, b in zip(cpu0, cpu1)])
        cpu = user + sys_
    io = _io_delta(_parse_counters(args.io0, 4), _parse_counters(args.io1, 4))

    peak_rss_kb, peak_pss_kb, max_procs, max_threads = _read_peak_file(args.peak_file)
    commands = read_rows(args.commands_tsv, args.stage, args.step)
    max_rss = _max_field(commands, "max_rss_mb")
    peak_mem = max(_mb_from_kb(peak_pss_kb), _max_field(commands, "peak_mem_mb") or 0)
    peak_rss = max(_mb_from_kb(peak_rss_kb), _max_field(commands, "peak_rss_mb") or 0)
    max_procs = max(max_procs, int(_max_field(commands, "max_procs") or 0))
    max_threads = max(max_threads, int(_max_field(commands, "max_threads") or 0))

    append_row(args.steps_tsv, STEP_COLUMNS, [
        args.stage, args.step, _iso(start), _iso(end),
        _num(wall), _num(user), _num(sys_), _num(cpu),
        _num(cpu / wall if cpu is not None and wall > 0 else None),
        _num(peak_mem, 1), _num(peak_rss, 1), _num(max_rss, 1),
        *(_num(v, 1) for v in io),
        max_procs, max_threads, len(commands),
    ])
    return 0


def cmd_total(args):
    start, end = _epoch(args.start), _epoch(args.end)
    wall = end - start
    user = sys_ = cpu = None
    ticks = _parse_counters(args.cpu, 2)
    if ticks:
        user, sys_ = _cpu_seconds(ticks)
        cpu = user + sys_
    io = _io_delta([0, 0, 0, 0], _parse_counters(args.io, 4))

    steps = [r for r in read_rows(args.steps_tsv, args.stage) if r["step"] != "TOTAL"]
    commands = read_rows(args.commands_tsv, args.stage)

    append_row(args.steps_tsv, STEP_COLUMNS, [
        args.stage, "TOTAL", _iso(start), _iso(end),
        _num(wall), _num(user), _num(sys_), _num(cpu),
        _num(cpu / wall if cpu is not None and wall > 0 else None),
        _num(_max_field(steps, "peak_mem_mb"), 1),
        _num(_max_field(steps, "peak_rss_mb"), 1),
        _num(_max_field(steps, "max_rss_mb"), 1),
        *(_num(v, 1) for v in io),
        int(_max_field(steps, "max_procs") or 0),
        int(_max_field(steps, "max_threads") or 0),
        len(commands),
    ])
    print_summary(args.stage, read_rows(args.steps_tsv, args.stage), commands,
                  args.steps_tsv, args.commands_tsv)
    return 0


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def _hms(seconds):
    if seconds in (None, "NA"):
        return "NA"
    seconds = int(round(float(seconds)))
    return f"{seconds // 3600}:{seconds % 3600 // 60:02d}:{seconds % 60:02d}"


def _gb(mb):
    if mb in (None, "NA"):
        return "NA"
    return f"{float(mb) / 1024:.2f} GB"


def print_summary(stage, steps, commands, steps_tsv, commands_tsv):
    print(f"Resource usage ({stage})")
    print(f"  {'step':<46} {'wall':>9} {'CPU':>9} {'cores':>6} "
          f"{'peak mem':>9} {'read':>9} {'written':>9}")
    for row in steps:
        name = row["step"]
        if name != "TOTAL":
            number, _, title = name.partition(" ")
            name = f"[{number}] {title}"
        if len(name) > 46:
            name = name[:45] + "…"
        print(f"  {name:<46} {_hms(row['wall_s']):>9} {_hms(row['cpu_s']):>9} "
              f"{row['avg_cores']:>6} {_gb(row['peak_mem_mb']):>9} "
              f"{_gb(row['read_mb']):>9} {_gb(row['write_mb']):>9}")

    tools = {}
    for row in commands:
        entry = tools.setdefault(row["tool"], {"calls": 0, "cpu": 0.0, "wall": 0.0, "mem": 0.0})
        entry["calls"] += 1
        entry["cpu"] += float(row["cpu_s"])
        entry["wall"] += float(row["wall_s"])
        entry["mem"] = max(entry["mem"], float(row["peak_mem_mb"]))
    if tools:
        print()
        print(f"  {'tool':<32} {'calls':>6} {'CPU':>9} {'run time':>9} {'peak mem':>9}")
        for name, entry in sorted(tools.items(), key=lambda kv: -kv[1]["cpu"]):
            print(f"  {name[:32]:<32} {entry['calls']:>6} {_hms(entry['cpu']):>9} "
                  f"{_hms(entry['wall']):>9} {_gb(entry['mem']):>9}")
        print("  (run time is summed over calls; parallel calls overlap in time)")
    print()
    print(f"  Per-step table   : {steps_tsv}")
    print(f"  Per-command table: {commands_tsv}")


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="mode", required=True)

    run = sub.add_parser("run", help="run and measure one command")
    run.add_argument("--log", required=True, help="TSV file to append to")
    run.add_argument("--stage", default="")
    run.add_argument("--step", default="")
    run.add_argument("--label", default="")
    run.add_argument("--interval", type=float, default=1.0,
                     help="memory sampling interval in seconds (default: 1)")
    run.add_argument("command", nargs=argparse.REMAINDER)

    smp = sub.add_parser("sample", help="sample a process tree until SIGTERM")
    smp.add_argument("--pid", type=int, required=True)
    smp.add_argument("--interval", type=float, default=1.0)
    smp.add_argument("--out", required=True)

    stp = sub.add_parser("step", help="append the row for one step")
    stp.add_argument("--steps-tsv", required=True)
    stp.add_argument("--commands-tsv", required=True)
    stp.add_argument("--stage", required=True)
    stp.add_argument("--step", required=True)
    stp.add_argument("--start", required=True)
    stp.add_argument("--end", required=True)
    stp.add_argument("--cpu0", default="")
    stp.add_argument("--cpu1", default="")
    stp.add_argument("--io0", default="")
    stp.add_argument("--io1", default="")
    stp.add_argument("--peak-file", required=True)

    tot = sub.add_parser("total", help="append the stage total and print a summary")
    tot.add_argument("--steps-tsv", required=True)
    tot.add_argument("--commands-tsv", required=True)
    tot.add_argument("--stage", required=True)
    tot.add_argument("--start", required=True)
    tot.add_argument("--end", required=True)
    tot.add_argument("--cpu", default="")
    tot.add_argument("--io", default="")

    args = parser.parse_args()
    handlers = {"run": cmd_run, "sample": cmd_sample, "step": cmd_step, "total": cmd_total}
    sys.exit(handlers[args.mode](args))


if __name__ == "__main__":
    main()
