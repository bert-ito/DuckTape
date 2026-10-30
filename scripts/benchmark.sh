# Resource monitoring helpers for the DuckTape stage scripts.
#
# Sourced by ducktape_stage1.sh and ducktape_stage2.sh once the conda
# environment is active. The measurements are made by ducktape_bench.py.
#
# Expects: BENCH_ON (1 or 0), BENCH_INTERVAL, BENCH_STAGE, BENCH_T_START,
#          OUTDIR, PIPELINE_DIR, and the warn() logging helper.
#
# Outputs (when BENCH_ON=1):
#   $OUTDIR/benchmark/commands.tsv  one row per measured command
#   $OUTDIR/benchmark/steps.tsv     one row per step, plus the stage TOTAL

BENCH_PY="$PIPELINE_DIR/ducktape_bench.py"
BENCH_DIR="$OUTDIR/benchmark"
BENCH_CMD_TSV="$BENCH_DIR/commands.tsv"
BENCH_STEP_TSV="$BENCH_DIR/steps.tsv"
BENCH_PEAK_FILE="$BENCH_DIR/.step_peak"
BENCH_STEP=""
BENCH_SAMPLER=""

export BENCH_ON BENCH_INTERVAL BENCH_STAGE BENCH_PY BENCH_CMD_TSV BENCH_STEP

# bench LABEL COMMAND [ARGS...]
#   Run COMMAND and append its resource usage, including every process it
#   starts, to commands.tsv. The first word of LABEL names the tool in the
#   summary. COMMAND runs unchanged when benchmarking is off.
bench() {
    local label=$1; shift
    if [[ "$BENCH_ON" != 1 ]]; then
        "$@"
        return
    fi
    python "$BENCH_PY" run --log "$BENCH_CMD_TSV" --stage "$BENCH_STAGE" \
        --step "$BENCH_STEP" --label "$label" --interval "$BENCH_INTERVAL" \
        -- "$@"
}
export -f bench

# Cumulative CPU ticks ("user,sys") and I/O bytes ("rchar,wchar,read,write")
# of this shell, including every child process it has waited for. Reads
# /proc directly so that no extra process is started.
bench_counters() {
    local stat key val
    local -a f
    read -r stat < /proc/$$/stat
    read -r -a f <<< "${stat##*) }"
    # Fields from state onwards: utime=11, stime=12, cutime=13, cstime=14
    BENCH_CPU="$(( f[11] + f[13] )),$(( f[12] + f[14] ))"
    BENCH_IO=""
    if [[ -r /proc/$$/io ]]; then
        while read -r key val; do
            case $key in
                rchar:|wchar:|read_bytes:|write_bytes:) BENCH_IO+="$val," ;;
            esac
        done < /proc/$$/io
    fi
}

bench_init() {
    [[ "$BENCH_ON" == 1 ]] || return 0
    rm -rf "$BENCH_DIR"
    mkdir -p "$BENCH_DIR"
    trap 'bench_cleanup; rm -f "$BENCH_PEAK_FILE"' EXIT
}

# Called by step_start: snapshot the counters and start sampling the memory
# of every process this shell runs during the step.
bench_step_begin() {
    printf -v BENCH_STEP '%02d %s' "$STEP_N" "$STEP_NAME"
    export BENCH_STEP
    [[ "$BENCH_ON" == 1 ]] || return 0
    BENCH_STEP_T0=${EPOCHREALTIME:-$(date +%s.%N)}
    bench_counters
    BENCH_CPU0=$BENCH_CPU
    BENCH_IO0=$BENCH_IO
    python "$BENCH_PY" sample --pid $$ --interval "$BENCH_INTERVAL" \
        --out "$BENCH_PEAK_FILE" &
    BENCH_SAMPLER=$!
}

# Called by step_end: take the end snapshot before the sampler is stopped,
# so the monitoring itself is not counted, then record the step.
bench_step_end() {
    [[ "$BENCH_ON" == 1 ]] || return 0
    local t1=${EPOCHREALTIME:-$(date +%s.%N)}
    bench_counters
    bench_cleanup
    python "$BENCH_PY" step --steps-tsv "$BENCH_STEP_TSV" \
        --commands-tsv "$BENCH_CMD_TSV" --stage "$BENCH_STAGE" \
        --step "$BENCH_STEP" --start "$BENCH_STEP_T0" --end "$t1" \
        --cpu0 "$BENCH_CPU0" --cpu1 "$BENCH_CPU" \
        --io0 "$BENCH_IO0" --io1 "$BENCH_IO" \
        --peak-file "$BENCH_PEAK_FILE" \
        || warn "could not record the resource usage of step $BENCH_STEP"
    rm -f "$BENCH_PEAK_FILE"
}

# Record the stage total (everything this shell and its children used since
# the script started) and print the summary tables.
bench_stage_end() {
    [[ "$BENCH_ON" == 1 ]] || return 0
    local t1=${EPOCHREALTIME:-$(date +%s.%N)}
    bench_counters
    python "$BENCH_PY" total --steps-tsv "$BENCH_STEP_TSV" \
        --commands-tsv "$BENCH_CMD_TSV" --stage "$BENCH_STAGE" \
        --start "$BENCH_T_START" --end "$t1" \
        --cpu "$BENCH_CPU" --io "$BENCH_IO" \
        || warn "could not record the stage resource usage"
}

# Stop the step sampler, if one is running.
bench_cleanup() {
    if [[ -n "$BENCH_SAMPLER" ]]; then
        kill -TERM "$BENCH_SAMPLER" 2>/dev/null || true
        wait "$BENCH_SAMPLER" 2>/dev/null || true
        BENCH_SAMPLER=""
    fi
}
