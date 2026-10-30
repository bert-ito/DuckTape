#!/usr/bin/env bash
############################################################
# STAGE 2: DuckTape NAT discovery
#
# Usage:
#   bash ducktape_stage2.sh -s1 stage1_results -b bam_dir -o stage2_results
#
# Resources: the script uses all CPUs available to the process (nproc) and
# the memory reported by /proc/meminfo. --threads and --mem-gb override both.
# Run with --help for the full list of options.
############################################################

############################################################
# Defaults
############################################################

# Start time, for the stage resource usage summary
BENCH_T_START=${EPOCHREALTIME:-$(date +%s.%N)}

ORPHAN_BUFFER=10000
ALPHA_FACTOR=4
INTERSECT_FRAC=0.2
TARGET_CLASS="all"
LIB_TYPE="RF"
READ_TYPE="PE"
MODE="analysis"
CONDA_ENV=""
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/scripts"

THREADS_OVERRIDE=""
MEM_GB_OVERRIDE=""
BENCH_ON=1
BENCH_INTERVAL=1

ANNOT_FILES=()

############################################################
# Help
############################################################

usage() {
cat << EOF

============================================================
DuckTape NAT Discovery - STAGE 2: NAT assembly
============================================================

Required arguments
------------------

  -s1, --stage1-dir DIR
        Output directory from the completed Stage 1 script. Must contain
        reference.gtf, annotation.gff3, orphans.bed, and featurecounts/*.sense.txt
        plus featurecounts/*.antisense.txt. Stage 1 / NASTIseq is NOT rerun.

  -b, --bam-dir DIR
        Directory containing *.bam files. Use the same BAM files as Stage 1.
        Bam files sorted by coordinate are preferred

  -o, --output DIR
        Output directory

Optional arguments
------------------

  --library-type RF|FR
        Strandedness of the RNA-seq library.

        RF  dUTP / fr-firststrand  (Illumina TruSeq stranded)
        FR  fr-secondstrand

        Default: RF

  --read-type PE|SE
        Paired-end or single-end reads

        Default: PE

  --target-class all|coding|noncoding
        Which orphan gene class should be used

        Default: all

  --orphan-buffer INT
        Extension around orphan loci before strand reversal
        Default: 10000

  --alpha-factor FLOAT
        Multiplier applied to estimated alpha
        Default: 4

  --intersect-frac FLOAT
        Minimum fractional overlap between reads and
        consensus NAT peaks
        Default: 0.2
        
  --mode analysis|test
        analysis: use annotation.gff3 and orphans.bed
        test:     use reference_test.gff3 and
                  TEST_KNOWNPAIRS_expression_filtered.bed

        Default: analysis

  --threads INT
        Number of CPUs to use
        Default: all CPUs available to the process (nproc)

  --mem-gb INT
        Memory available to the pipeline (GB)
        Default: MemAvailable from /proc/meminfo

  --env NAME|PATH
        Conda environment to activate before running.
        Default: use the currently active environment

  --pipeline-dir DIR
        Directory containing the helper scripts
        Default: scripts/ next to this script

  --no-benchmark
        Do not record resource usage (CPU, memory, I/O)

  --benchmark-interval SEC
        Memory sampling interval for the resource usage tables
        Default: 1

  -h, --help
        Show this help message and exit

Outputs
-------

  analysis mode:
    FINAL_NATS.gtf          Final NAT transcript models
    orphans.bed             NASTIseq orphan loci used as targets
    NAT_tracking_x.tsv      gffcompare antisense (class code x) assignments

  test mode:
    FINAL_QUERY.gtf         NAT models recovered for the removed query genes
    test_knownpairs.bed     Test genes used as targets
    NAT_tracking_x.tsv      gffcompare antisense (class code x) assignments
    final_queries.stats     gffcompare accuracy against the query genes

  both modes:
    benchmark/steps.tsv     Resource usage per step and for the stage
    benchmark/commands.tsv  Resource usage per command

Examples
--------

  # Standard analysis using completed Stage 1 results
  bash ducktape_stage2.sh \\
      --stage1-dir /path/to/stage1_results \\
      --bam-dir /path/to/bams \\
      --output /path/to/stage2_results \\
      --mode analysis

  # Test mode with query tracking and comparison metrics
  bash ducktape_stage2.sh \\
      --stage1-dir /path/to/stage1_results \\
      --bam-dir /path/to/bams \\
      --output /path/to/test_results \\
      --mode test \\
      --target-class all

  # Override library type, read type, and resource detection
  bash ducktape_stage2.sh \\
      --stage1-dir /path/to/stage1_results \\
      --bam-dir /path/to/bams \\
      --output /path/to/stage2_results \\
      --library-type FR \\
      --read-type SE \\
      --threads 16 \\
      --mem-gb 64 \\
      --env /path/to/conda/environment


============================================================

EOF
}

############################################################
# Parse args
############################################################

while [[ $# -gt 0 ]]
do
    case $1 in
        -h|--help)          usage; exit 0 ;;
        -s1|--stage1-dir)   STAGE1_DIR="$2"; shift 2 ;;
        -b|--bam-dir)       BAM_DIR="$2";        shift 2 ;;
        -o|--output)        OUTDIR="$2";         shift 2 ;;
        --library-type)     LIB_TYPE="${2^^}";   shift 2 ;;
        --read-type)        READ_TYPE="${2^^}";  shift 2 ;;
        --target-class)     TARGET_CLASS="$2";   shift 2 ;;
        --orphan-buffer)    ORPHAN_BUFFER="$2";  shift 2 ;;
        --alpha-factor)     ALPHA_FACTOR="$2";   shift 2 ;;
        --intersect-frac)   INTERSECT_FRAC="$2"; shift 2 ;;
        --mode)             MODE="${2,,}"; shift 2 ;;
        --threads)          THREADS_OVERRIDE="$2"; shift 2 ;;
        --mem-gb)           MEM_GB_OVERRIDE="$2";  shift 2 ;;
        --env)              CONDA_ENV="$2";      shift 2 ;;
        --pipeline-dir)     PIPELINE_DIR="$2";   shift 2 ;;
        --no-benchmark)     BENCH_ON=0;          shift ;;
        --benchmark-interval)
                            BENCH_INTERVAL="$2"; shift 2 ;;
        *)                  echo "Unknown argument: $1"; exit 1 ;;
    esac
done

############################################################
# Logging helpers
############################################################

STEP_N=0
STEP_NAME=""
STEP_T0=0
PIPELINE_T0=$(date +%s)
declare -a TIMINGS=()

fmt_time() {
    local s=$1
    if   (( s < 60 ));   then printf '%ds' "$s"
    elif (( s < 3600 )); then printf '%dm %02ds' $((s/60)) $((s%60))
    else                      printf '%dh %02dm %02ds' $((s/3600)) $(((s%3600)/60)) $((s%60))
    fi
}

log()  { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*"; }
warn() { printf '[%s] WARNING: %s\n' "$(date '+%H:%M:%S')" "$*" >&2; }
die()  { printf '\n[%s] ERROR: %s\n\n' "$(date '+%H:%M:%S')" "$*" >&2; exit 1; }

rule() { printf '%.0s─' $(seq 1 74); printf '\n'; }

step_start() {
    STEP_N=$((STEP_N+1))
    STEP_NAME="$1"
    STEP_T0=$(date +%s)
    printf '\n'
    printf '── [%02d] %s ' "$STEP_N" "$STEP_NAME"
    local pad=$(( 68 - ${#STEP_NAME} ))
    (( pad > 0 )) && printf '%.0s─' $(seq 1 $pad)
    printf '\n'
    bench_step_begin
}

step_end() {
    bench_step_end
    local dt=$(( $(date +%s) - STEP_T0 ))
    printf '     done in %s\n' "$(fmt_time "$dt")"
    TIMINGS+=("$(printf '%-46s %s' "[$(printf '%02d' $STEP_N)] $STEP_NAME" "$(fmt_time "$dt")")")
}

############################################################
# Validate required args
############################################################

if [[ -z "${STAGE1_DIR:-}" || -z "${BAM_DIR:-}" || -z "${OUTDIR:-}" ]]; then
    echo "ERROR: Missing required arguments"
    usage
    exit 1
fi

case "$LIB_TYPE" in
    RF|FR) ;;
    *) echo "ERROR: --library-type must be RF or FR"; exit 1 ;;
esac

case "$READ_TYPE" in
    PE|SE) ;;
    *) echo "ERROR: --read-type must be PE or SE"; exit 1 ;;
esac

case "$TARGET_CLASS" in
    all|coding|noncoding) ;;
    *) echo "ERROR: --target-class must be all, coding or noncoding"; exit 1 ;;
esac

case "$MODE" in
    analysis|test) ;;
    *) die "--mode must be analysis or test" ;;
esac

awk -v x="$BENCH_INTERVAL" \
    'BEGIN { exit !(x ~ /^[0-9]+([.][0-9]+)?$/ && x > 0) }' \
    || die "--benchmark-interval must be a positive number"

mkdir -p "$OUTDIR"

############################################################
# Environment
############################################################

set -euo pipefail

# Optionally activate a conda environment; otherwise use the current one.
# conda's activation scripts are not 'set -u' safe, so relax it briefly.
if [[ -n "$CONDA_ENV" ]]; then
    command -v conda >/dev/null 2>&1 \
        || die "--env was given but 'conda' is not on PATH"
    set +u
    eval "$(conda shell.bash hook)"
    conda activate "$CONDA_ENV" \
        || die "could not activate conda environment: $CONDA_ENV"
    set -u
fi

############################################################
# Resource detection
############################################################

detect_threads() {
    if [[ -n "$THREADS_OVERRIDE" ]]; then echo "$THREADS_OVERRIDE"; return; fi
    nproc 2>/dev/null || echo 1
}

detect_mem_mb() {
    if [[ -n "$MEM_GB_OVERRIDE" ]]; then echo $(( MEM_GB_OVERRIDE * 1024 )); return; fi
    awk '/MemAvailable/ {printf "%d", $2/1024}' /proc/meminfo 2>/dev/null || echo 8192
}

THREADS=$(detect_threads)
MEM_MB=$(detect_mem_mb)

[[ "$THREADS" =~ ^[0-9]+$ ]] || THREADS=1
(( THREADS < 1 )) && THREADS=1
[[ "$MEM_MB" =~ ^[0-9]+$ ]] || MEM_MB=8192

# Reserve 20% of memory as headroom
MEM_USABLE_MB=$(( MEM_MB * 80 / 100 ))

############################################################
# Adaptive scheduler
#
# schedule <n_jobs> <pref_threads> <cap_threads>
#
#   Maximises the number of concurrent jobs at the tool's preferred
#   thread width, then widens each job only if cores remain unused.
#   Sets SCHED_CONC (concurrent jobs) and SCHED_THREADS (per job).
############################################################

SCHED_CONC=1
SCHED_THREADS=1

schedule() {
    local njobs=$1 pref=$2 cap=$3
    local conc per

    (( njobs < 1 )) && njobs=1
    (( pref  < 1 )) && pref=1
    (( cap   < 1 )) && cap=1

    conc=$(( THREADS / pref ))
    (( conc < 1 ))      && conc=1
    (( conc > njobs ))  && conc=$njobs

    per=$(( THREADS / conc ))
    (( per > cap )) && per=$cap
    (( per < 1 ))   && per=1

    SCHED_CONC=$conc
    SCHED_THREADS=$per
}

# Memory budget for the current schedule, in MB per job
sched_mem_mb() { echo $(( MEM_USABLE_MB / SCHED_CONC )); }

# Run one exported bash function across a list of items.
#   run_parallel <func_name> <item> [item ...]
run_parallel() {
    local func=$1; shift
    (( $# == 0 )) && return 0
    printf '%s\n' "$@" \
        | xargs -P "$SCHED_CONC" -I '{}' \
              bash -euo pipefail -c "$func \"\$@\"" _ '{}'
}

############################################################
# Tool availability
############################################################

for tool in gffread gffcompare stringtie samtools bedtools \
            bamCoverage bigWigToBedGraph python xargs
do
    command -v "$tool" >/dev/null 2>&1 \
        || die "'$tool' not found on PATH; activate the DuckTape conda environment or use --env"
done

SAMTOOLS_VER=$(samtools --version | head -1 | awk '{print $2}')
SAMTOOLS_MAJ=${SAMTOOLS_VER%%.*}
SAMTOOLS_MIN=${SAMTOOLS_VER#*.}; SAMTOOLS_MIN=${SAMTOOLS_MIN%%.*}
if (( SAMTOOLS_MAJ < 1 )) || { (( SAMTOOLS_MAJ == 1 )) && (( SAMTOOLS_MIN < 12 )); }; then
    die "samtools >= 1.12 required (filter expressions); found $SAMTOOLS_VER"
fi

############################################################
# Helper scripts
############################################################

ISOLATE="$PIPELINE_DIR/isolate_gff3.py"
TOTAL_ALPHA="$PIPELINE_DIR/total_alpha.py"
PEARSON="$PIPELINE_DIR/pearson_alpha.py"
ALPHA="$PIPELINE_DIR/calculate_alpha.py"
REMOVE="$PIPELINE_DIR/remove_alpha.py"
BENCH_LIB="$PIPELINE_DIR/benchmark.sh"
for f in "$ISOLATE" "$TOTAL_ALPHA" "$PEARSON" \
         "$ALPHA" "$REMOVE" "$BENCH_LIB" "$PIPELINE_DIR/ducktape_bench.py"
do
    [[ -f "$f" ]] || die "helper script not found: $f
       Set --pipeline-dir to the directory containing them."
done

BENCH_STAGE="stage2"
source "$BENCH_LIB"


############################################################
# Validate inputs and select mode-specific files
############################################################

[[ -d "$STAGE1_DIR" ]] || die "Stage 1 directory not found: $STAGE1_DIR"
[[ -d "$BAM_DIR" ]] || die "BAM directory not found: $BAM_DIR"

COUNTS="$STAGE1_DIR/featurecounts"

case "$MODE" in
    analysis)
        ANNOT="$STAGE1_DIR/annotation.gff3"
        ORPHANS="$STAGE1_DIR/orphans.bed"
        ;;
    test)
        ANNOT="$STAGE1_DIR/reference_test.gff3"
        ORPHANS="$STAGE1_DIR/TEST_KNOWNPAIRS_expression_filtered.bed"
        ;;
esac

[[ -s "$ANNOT" ]] \
    || die "annotation not found or empty for mode '$MODE': $ANNOT"

[[ -s "$ORPHANS" ]] \
    || die "orphan/test BED not found or empty for mode '$MODE': $ORPHANS"

[[ -d "$COUNTS" ]] \
    || die "Stage 1 featurecounts directory not found: $COUNTS"

# The reference GTF is not used by the downstream analysis steps.
# Keep it as an optional informational path for analysis mode.
if [[ "$MODE" == "analysis" ]]; then
    GTF="$STAGE1_DIR/reference.gtf"
    [[ -s "$GTF" ]] || die "reference GTF not found or empty: $GTF"
else
    GTF=""
fi


############################################################
# Library-type dependent settings
############################################################

if [[ "$LIB_TYPE" == "RF" ]]; then
    FC_SENSE=2; FC_ANTI=1
    BW_PLUS="forward"; BW_MINUS="reverse"
    ST_STRAND="--rf"
    if [[ "$READ_TYPE" == "PE" ]]; then
        PLUS_EXPR='(flag.read1 && flag.reverse) || (flag.read2 && !flag.reverse)'
        MINUS_EXPR='(flag.read1 && !flag.reverse) || (flag.read2 && flag.reverse)'
    else
        PLUS_EXPR='flag.reverse'
        MINUS_EXPR='!flag.reverse'
    fi
else
    FC_SENSE=1; FC_ANTI=2
    BW_PLUS="reverse"; BW_MINUS="forward"
    ST_STRAND="--fr"
    if [[ "$READ_TYPE" == "PE" ]]; then
        PLUS_EXPR='(flag.read1 && !flag.reverse) || (flag.read2 && flag.reverse)'
        MINUS_EXPR='(flag.read1 && flag.reverse) || (flag.read2 && !flag.reverse)'
    else
        PLUS_EXPR='!flag.reverse'
        MINUS_EXPR='flag.reverse'
    fi
fi

if [[ "$READ_TYPE" == "PE" ]]; then
    FC_PAIRED=(-p --countReadPairs -B -C)
else
    FC_PAIRED=()
fi

# Discover BAM samples and initialize the sample count.
shopt -s nullglob
BAM_FILES=("$BAM_DIR"/*.bam)
shopt -u nullglob

((${#BAM_FILES[@]} > 0)) \
    || die "no .bam files found in BAM directory: $BAM_DIR"

SAMPLES=()
for bam in "${BAM_FILES[@]}"; do
    sample=$(basename "$bam" .bam)
    SAMPLES+=("$sample")
done

N_SAMPLES=${#SAMPLES[@]}

(( N_SAMPLES > 0 )) || die "no BAM samples detected"
log "Detected $N_SAMPLES BAM samples"

############################################################
# Working directories
############################################################

TMP="$OUTDIR/tmp_work"
rm -rf "$TMP"
mkdir -p "$TMP"

LOGS="$TMP/logs"
COUNTS="$STAGE1_DIR/featurecounts"
ALPHA_DIR="$TMP/alpha"
BW="$TMP/bigwigs"
BW_CORR="$TMP/corrected_bw"
PEAKDIR="$TMP/peaks"
ASSEMBLY="$TMP/assembly"
ISOLATED_COUNTS="$TMP/isolated_counts"
TOTAL_ISOLATED_COUNTS="$TMP/total_isolated_counts"
NAT_COUNTS="$TMP/nat_counts"
CORRELATION="$TMP/correlation"
REST_BAMS="$TMP/restricted_bams"
SAMPLE_DIR="$TMP/samples"

mkdir -p "$LOGS" "$ALPHA_DIR" "$BW" "$BW_CORR" \
         "$PEAKDIR" "$ASSEMBLY" "$ISOLATED_COUNTS" "$TOTAL_ISOLATED_COUNTS" \
         "$NAT_COUNTS" "$CORRELATION" "$REST_BAMS" "$SAMPLE_DIR"

# GNU sort tuning
SORT_MEM=$(( MEM_USABLE_MB / 4 ))
(( SORT_MEM > 8192 )) && SORT_MEM=8192
(( SORT_MEM < 256 ))  && SORT_MEM=256
SORT_OPTS=(--parallel="$THREADS" -S "${SORT_MEM}M" -T "$TMP")

############################################################
# Run header
############################################################



rule
printf 'DuckTape NAT Discovery - STAGE 2: NAT assembly\n'
rule
printf '  run started     : %s\n' "$(date '+%Y-%m-%d %H:%M:%S %Z')"
printf '  mode            : %s\n' "$MODE"
printf '  Stage 1 inputs  : %s\n' "$STAGE1_DIR"
printf '  selected GFF    : %s\n' "$(basename "$ANNOT")"
printf '  selected BED    : %s\n' "$(basename "$ORPHANS")"
printf '  reused counts   : %s\n' "$COUNTS"
printf '  BAM directory   : %s\n' "$BAM_DIR"
printf '  samples         : %d\n' "$N_SAMPLES"
printf '  output          : %s\n' "$OUTDIR"
printf '  library type    : %s (%s)\n' "$LIB_TYPE" "$READ_TYPE"
printf '  target class    : %s\n' "$TARGET_CLASS"
printf '  orphan buffer   : %s bp\n' "$ORPHAN_BUFFER"
printf '  alpha factor    : %s\n' "$ALPHA_FACTOR"
printf '  intersect frac  : %s\n' "$INTERSECT_FRAC"
printf '  CPUs            : %d\n' "$THREADS"
printf '  memory          : %d GB (%d GB usable)\n' \
       $(( MEM_MB / 1024 )) $(( MEM_USABLE_MB / 1024 ))
rule

# Show the computed execution plan
printf 'Execution plan (concurrent jobs x threads each)\n'
printf '  featureCounts       reused from Stage 1\n'
schedule "$N_SAMPLES"         4 8; printf '  BAM restriction     %2d x %-2d\n' "$SCHED_CONC" "$SCHED_THREADS"
schedule $(( N_SAMPLES * 2 )) 4 8; printf '  bamCoverage         %2d x %-2d\n' "$SCHED_CONC" "$SCHED_THREADS"
schedule "$N_SAMPLES"         1 1; printf '  alpha correction    %2d x %-2d\n' "$SCHED_CONC" "$SCHED_THREADS"
schedule "$N_SAMPLES"         1 2; printf '  peak calling        %2d x %-2d\n' "$SCHED_CONC" "$SCHED_THREADS"
schedule "$N_SAMPLES"         4 8; printf '  read extraction     %2d x %-2d\n' "$SCHED_CONC" "$SCHED_THREADS"
schedule 1                    8 8; printf '  merge / sort        %2d x %-2d\n' "$SCHED_CONC" "$SCHED_THREADS"
printf '  stringtie            1 x %-2d\n' "$THREADS"
rule
printf 'Tool versions\n'
printf '  samtools   %s\n' "$SAMTOOLS_VER"
printf '  stringtie  %s\n' "$(stringtie --version 2>&1 | head -1)"
rule
printf 'Working directory: %s\n' "$TMP"
printf 'Tool logs        : %s\n' "$LOGS"

############################################################
# Export everything the parallel workers need
############################################################

export TMP LOGS COUNTS ALPHA_DIR BW BW_CORR PEAKDIR ASSEMBLY
export ISOLATED_COUNTS TOTAL_ISOLATED_COUNTS NAT_COUNTS CORRELATION
export REST_BAMS SAMPLE_DIR BAM_DIR OUTDIR
export FC_SENSE FC_ANTI BW_PLUS BW_MINUS PLUS_EXPR MINUS_EXPR
export INTERSECT_FRAC ALPHA_FACTOR
export ISOLATE TOTAL_ALPHA PEARSON ALPHA REMOVE

bench_init

############################################################
# Reuse completed Stage 1 outputs (no annotation prep, featureCounts,
# or NASTIseq rerun)
############################################################

step_start "Loading completed Stage 1 results"

# ANNOT, GTF, ORPHANS and COUNTS were validated and set above.
# Verify the count tables required by this BAM set are present.
for sample in "${SAMPLES[@]}"; do
    [[ -s "$COUNTS/${sample}.sense.txt" ]] \
        || die "missing Stage 1 sense counts for sample $sample: $COUNTS/${sample}.sense.txt"
    [[ -s "$COUNTS/${sample}.antisense.txt" ]] \
        || die "missing Stage 1 antisense counts for sample $sample: $COUNTS/${sample}.antisense.txt"
done

if [[ "$MODE" == "analysis" ]]; then
    cp "$ORPHANS" "$OUTDIR/orphans.bed"
else
    cp "$ORPHANS" "$OUTDIR/test_knownpairs.bed"
fi

export ANNOT ORPHANS COUNTS

log "mode: $MODE"
log "reusing annotation: $ANNOT"
log "reusing featureCounts: $COUNTS"
log "reusing selected BED: $ORPHANS"
log "selected BED entries: $(wc -l < "$ORPHANS")"
step_end

############################################################
# STEP 2 - Target annotation and orphan sets
############################################################

step_start "Building target orphan set"

TARGET_GFF="$TMP/target.gff3"

case "$TARGET_CLASS" in
    all)       cp "$ANNOT" "$TARGET_GFF" ;;
    coding)    bench "gffread coding targets" \
                   gffread "$ANNOT" -C  --keep-genes -o "$TARGET_GFF" \
                   >"$LOGS/gffread_target.log" 2>&1 ;;
    noncoding) bench "gffread noncoding targets" \
                   gffread "$ANNOT" --nc --keep-genes -o "$TARGET_GFF" \
                   >"$LOGS/gffread_target.log" 2>&1 ;;
esac

awk -F'\t' '
$3=="gene"{
    id=""
    if(match($9,/gene_id=[^;]+/))       id=substr($9,RSTART+8,RLENGTH-8)
    else if(match($9,/geneID=[^;]+/))   id=substr($9,RSTART+7,RLENGTH-7)
    else if(match($9,/ID=gene:[^;]+/))  id=substr($9,RSTART+8,RLENGTH-8)
    else if(match($9,/ID=[^;]+/))       id=substr($9,RSTART+3,RLENGTH-3)
    if(id!="") print id
}' "$TARGET_GFF" > "$TMP/target_genes.txt"

awk 'FNR==NR{keep[$1]; next} $4 in keep' \
    "$TMP/target_genes.txt" "$ORPHANS" > "$TMP/target_orphans.bed"

cut -f4 "$TMP/target_orphans.bed" > "$TMP/orphan_ids.txt"

awk -v B="$ORPHAN_BUFFER" '
BEGIN{OFS="\t"}
{
    $2 -= B; if($2 < 0) $2 = 0
    $3 += B
    if($6=="+") $6="-"; else if($6=="-") $6="+"
    print
}' "$TMP/target_orphans.bed" > "$TMP/target_orphans_rev.bed"

awk '
BEGIN{FS=OFS="\t"}
function get_id(attr){
    if(match(attr,/gene_id=[^;]+/))  return substr(attr,RSTART+8,RLENGTH-8)
    if(match(attr,/geneID=[^;]+/))   return substr(attr,RSTART+7,RLENGTH-7)
    if(match(attr,/ID=gene:[^;]+/))  return substr(attr,RSTART+8,RLENGTH-8)
    if(match(attr,/ID=[^;]+/))       return substr(attr,RSTART+3,RLENGTH-3)
    return ""
}
FNR==NR{orphan[$1]; next}
/^#/{print; next}
$3=="gene"{ keep=0; if(get_id($9) in orphan) keep=1 }
keep{print}
' "$TMP/orphan_ids.txt" "$TARGET_GFF" > "$TMP/target_orphans.gff3"

log "target orphans: $(wc -l < "$TMP/target_orphans.bed")"
step_end

############################################################
# STEP 3 - Isolated genes and per-sample alpha
############################################################

step_start "Estimating antisense leakage (alpha)"

bench "isolate_gff3.py" python "$ISOLATE" -i "$ANNOT" -b 200 -o "$TMP/isolated.txt" \
    > "$LOGS/isolate.log" 2>&1 || die "isolate_gff3.py failed"

grep -Fvxf "$TMP/orphan_ids.txt" "$TMP/isolated.txt" \
    > "$TMP/isolated_no_orphans.txt" || true

log "isolated non-orphan genes: $(wc -l < "$TMP/isolated_no_orphans.txt")"

subset_counts() {
    local keep="$1" infile="$2" outfile="$3"
    awk '
    BEGIN{FS=OFS="\t"}
    NR==FNR {keep[$1]; next}
    FNR<=2  {print; next}
    { gene=$1; sub(/^gene:/,"",gene); if (gene in keep) print }
    ' "$keep" "$infile" > "$outfile"
}
export -f subset_counts

alpha_one() {
    local sample="$1"
    subset_counts "$TMP/isolated_no_orphans.txt" \
        "$COUNTS/${sample}.sense.txt" \
        "$ISOLATED_COUNTS/${sample}.sense.isolated.txt"
    subset_counts "$TMP/isolated_no_orphans.txt" \
        "$COUNTS/${sample}.antisense.txt" \
        "$ISOLATED_COUNTS/${sample}.antisense.isolated.txt"
    bench "calculate_alpha.py ${sample}" python "$ALPHA" \
        -s "$ISOLATED_COUNTS/${sample}.sense.isolated.txt" \
        -a "$ISOLATED_COUNTS/${sample}.antisense.isolated.txt" \
        -m 50 -r 0.5 \
        -o "$ALPHA_DIR/${sample}.alpha.txt" \
        > "$LOGS/alpha_${sample}.log" 2>&1
}
export -f alpha_one

schedule "$N_SAMPLES" 1 1
log "$SCHED_CONC concurrent jobs x $SCHED_THREADS thread"
run_parallel alpha_one "${SAMPLES[@]}" || die "alpha estimation failed"

step_end

############################################################
# STEP 4 - Pooled alpha vs expression
############################################################

step_start "Alpha-expression correlation"

merge_fc() {
    local outfile=$1; shift
    awk '
    BEGIN{OFS="\t"}
    FNR==1 && NR==1 {comment=$0; next}
    FNR==2 && NR==2 {header=$0;  next}
    FNR<=2 {next}
    {
        if (!($1 in meta))
            meta[$1]=$1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6
        sum[$1]+=$7
    }
    END{
        print comment
        print "Geneid","Chr","Start","End","Strand","Length","Counts"
        for (g in meta) print meta[g],sum[g]
    }' "$@" \
    | ( read -r h1; read -r h2; printf '%s\n%s\n' "$h1" "$h2"; sort -k1,1 ) \
    > "$outfile"
}

merge_fc "$TOTAL_ISOLATED_COUNTS/merged.sense.isolated.txt" \
         "$ISOLATED_COUNTS"/*.sense.isolated.txt
merge_fc "$TOTAL_ISOLATED_COUNTS/merged.antisense.isolated.txt" \
         "$ISOLATED_COUNTS"/*.antisense.isolated.txt

bench "total_alpha.py" python "$TOTAL_ALPHA" \
    --sense     "$TOTAL_ISOLATED_COUNTS/merged.sense.isolated.txt" \
    --antisense "$TOTAL_ISOLATED_COUNTS/merged.antisense.isolated.txt" \
    --output    "$CORRELATION/alpha_median_vs_expression.tsv" \
    > "$LOGS/total_alpha.log" 2>&1 || die "total_alpha.py failed"

bench "pearson_alpha.py" python "$PEARSON" \
    --input  "$CORRELATION/alpha_median_vs_expression.tsv" \
    --output "$CORRELATION/expression_alpha_correlation.txt" \
    > "$LOGS/pearson.log" 2>&1 || die "pearson_alpha.py failed"

sed 's/^/     /' "$CORRELATION/expression_alpha_correlation.txt"

step_end

############################################################
# STEP 5 - Restrict BAMs to reversed orphan regions
############################################################

step_start "Restricting BAMs to orphan regions"

schedule "$N_SAMPLES" 4 8
BAM_THREADS=$SCHED_THREADS
BAM_MEM=$(( $(sched_mem_mb) / (BAM_THREADS > 0 ? BAM_THREADS : 1) ))
(( BAM_MEM > 4096 )) && BAM_MEM=4096
(( BAM_MEM < 256 ))  && BAM_MEM=256
export BAM_THREADS BAM_MEM

log "$SCHED_CONC concurrent jobs x $BAM_THREADS threads (${BAM_MEM}M/thread sort)"

restrict_one() {
    local sample="$1"
    bench "samtools view restrict ${sample}" \
        samtools view -@ "$BAM_THREADS" -b \
        -L "$TMP/target_orphans_rev.bed" "$BAM_DIR/${sample}.bam" \
    | bench "samtools sort restrict ${sample}" \
        samtools sort -@ "$BAM_THREADS" -m "${BAM_MEM}M" \
        -T "$TMP/sort_${sample}" -o "$REST_BAMS/${sample}.bam" -
    bench "samtools index restrict ${sample}" \
        samtools index -@ "$BAM_THREADS" "$REST_BAMS/${sample}.bam"
}
export -f restrict_one

run_parallel restrict_one "${SAMPLES[@]}" || die "BAM restriction failed"
step_end

############################################################
# STEP 6 - Coverage tracks and alpha correction
############################################################

step_start "Building corrected coverage tracks"

schedule $(( N_SAMPLES * 2 )) 4 8
BW_THREADS=$SCHED_THREADS
export BW_THREADS
log "bamCoverage: $SCHED_CONC concurrent jobs x $BW_THREADS threads"

coverage_one() {
    local sample="$1"
    bench "bamCoverage forward ${sample}" bamCoverage -b "$REST_BAMS/${sample}.bam" \
        --filterRNAstrand "$BW_PLUS" --binSize 50 -p "$BW_THREADS" \
        -o "$BW/${sample}.forward.bw" \
        > "$LOGS/bamCoverage_${sample}_fwd.log" 2>&1
    bench "bamCoverage reverse ${sample}" bamCoverage -b "$REST_BAMS/${sample}.bam" \
        --filterRNAstrand "$BW_MINUS" --binSize 50 -p "$BW_THREADS" \
        -o "$BW/${sample}.reverse.bw" \
        > "$LOGS/bamCoverage_${sample}_rev.log" 2>&1
}
export -f coverage_one
run_parallel coverage_one "${SAMPLES[@]}" || die "bamCoverage failed"

# remove_alpha.py is single-threaded: run as many samples at once as we have cores
schedule "$N_SAMPLES" 1 1
log "alpha correction: $SCHED_CONC concurrent jobs x 1 thread"

correct_one() {
    local sample="$1"
    local a f
    a=$(cat "$ALPHA_DIR/${sample}.alpha.txt")
    f=$(awk -v a="$a" -v k="$ALPHA_FACTOR" 'BEGIN{print a*k}')
    bench "remove_alpha.py ${sample}" python "$REMOVE" \
        -p  "$BW/${sample}.forward.bw" \
        -m  "$BW/${sample}.reverse.bw" \
        -op "$BW_CORR/${sample}.forward.bw" \
        -om "$BW_CORR/${sample}.reverse.bw" \
        -a  "$f" \
        > "$LOGS/remove_alpha_${sample}.log" 2>&1
}
export -f correct_one
run_parallel correct_one "${SAMPLES[@]}" || die "alpha correction failed"
step_end

############################################################
# STEP 7 - Peak calling
############################################################

step_start "Calling consensus NAT peaks"

schedule "$N_SAMPLES" 1 2
PEAK_THREADS=$SCHED_THREADS
PEAK_SORT_MEM=$(( $(sched_mem_mb) ))
(( PEAK_SORT_MEM > 4096 )) && PEAK_SORT_MEM=4096
(( PEAK_SORT_MEM < 128 ))  && PEAK_SORT_MEM=128
export PEAK_THREADS PEAK_SORT_MEM
log "$SCHED_CONC concurrent jobs x $PEAK_THREADS threads"

peaks_one() {
    local sample="$1"
    local out="$PEAKDIR/${sample}.bed"
    : > "$out"

    bench "bigWigToBedGraph forward ${sample}" \
        bigWigToBedGraph "$BW_CORR/${sample}.forward.bw" stdout \
      | awk '$4>10' \
      | bench "sort peaks forward ${sample}" \
        sort -k1,1 -k2,2n --parallel="$PEAK_THREADS" -S "${PEAK_SORT_MEM}M" -T "$TMP" \
      | bench "bedtools merge peaks forward ${sample}" bedtools merge -i - \
      | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,".",".","+"}' >> "$out"

    bench "bigWigToBedGraph reverse ${sample}" \
        bigWigToBedGraph "$BW_CORR/${sample}.reverse.bw" stdout \
      | awk '$4>10' \
      | bench "sort peaks reverse ${sample}" \
        sort -k1,1 -k2,2n --parallel="$PEAK_THREADS" -S "${PEAK_SORT_MEM}M" -T "$TMP" \
      | bench "bedtools merge peaks reverse ${sample}" bedtools merge -i - \
      | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,".",".","-"}' >> "$out"
}
export -f peaks_one
run_parallel peaks_one "${SAMPLES[@]}" || die "peak calling failed"

cat "$PEAKDIR"/*.bed | bench "sort all peaks" sort -k1,1 -k2,2n "${SORT_OPTS[@]}" \
    > "$TMP/all_peaks.sorted.bed"

bench "bedtools merge all peaks" \
    bedtools merge -i "$TMP/all_peaks.sorted.bed" -s -c 6 -o distinct \
  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,".",".",$4}' \
  > "$TMP/merged_peaks.bed"

bench "bedtools intersect consensus peaks" bedtools intersect \
    -a "$TMP/merged_peaks.bed" \
    -b "$TMP/target_orphans_rev.bed" \
    -s -f 0.1 \
    > "$TMP/consensus_peaks.bed"

awk '$6=="+"' "$TMP/consensus_peaks.bed" > "$TMP/plus.bed"
awk '$6=="-"' "$TMP/consensus_peaks.bed" > "$TMP/minus.bed"

log "consensus peaks: $(wc -l < "$TMP/consensus_peaks.bed")"
step_end

############################################################
# STEP 8 - Strand-aware read extraction
############################################################

step_start "Extracting strand-specific reads"

schedule "$N_SAMPLES" 4 8
EXTRACT_THREADS=$SCHED_THREADS
export EXTRACT_THREADS
log "$SCHED_CONC concurrent jobs x $EXTRACT_THREADS threads"

extract_one() {
    local sample="$1"
    local sd="$SAMPLE_DIR/$sample"
    mkdir -p "$sd"

    # One pass per strand using a filter expression, replacing the
    # previous two view calls + samtools cat per strand.
    bench "samtools view plus ${sample}" \
        samtools view -@ "$EXTRACT_THREADS" -b -e "$PLUS_EXPR" \
        "$REST_BAMS/${sample}.bam" > "$sd/plus.bam"
    bench "samtools view minus ${sample}" \
        samtools view -@ "$EXTRACT_THREADS" -b -e "$MINUS_EXPR" \
        "$REST_BAMS/${sample}.bam" > "$sd/minus.bam"

    bench "bedtools intersect plus ${sample}" \
        bedtools intersect -abam "$sd/plus.bam"  -b "$TMP/plus.bed" \
        -f "$INTERSECT_FRAC" > "$sd/plus_hits.bam"
    bench "bedtools intersect minus ${sample}" \
        bedtools intersect -abam "$sd/minus.bam" -b "$TMP/minus.bed" \
        -f "$INTERSECT_FRAC" > "$sd/minus_hits.bam"

    bench "samtools cat ${sample}" \
        samtools cat -o "$sd/all.bam" "$sd/plus_hits.bam" "$sd/minus_hits.bam"

    # Free intermediates as soon as they are consumed
    rm -f "$sd/plus.bam" "$sd/minus.bam" "$sd/plus_hits.bam" "$sd/minus_hits.bam"
}
export -f extract_one
run_parallel extract_one "${SAMPLES[@]}" || die "read extraction failed"

MERGELIST="$TMP/merge_bams.txt"
: > "$MERGELIST"
for sample in "${SAMPLES[@]}"; do
    echo "$SAMPLE_DIR/$sample/all.bam" >> "$MERGELIST"
done
step_end

############################################################
# STEP 9 - Merge and assemble
############################################################

step_start "Merging BAMs and assembling transcripts"

schedule 1 8 8
MERGE_THREADS=$SCHED_THREADS
MERGE_MEM=$(( MEM_USABLE_MB / (MERGE_THREADS > 0 ? MERGE_THREADS : 1) ))
(( MERGE_MEM > 8192 )) && MERGE_MEM=8192
(( MERGE_MEM < 256 ))  && MERGE_MEM=256

log "samtools merge/sort: $MERGE_THREADS threads, ${MERGE_MEM}M per thread"

bench "samtools merge all samples" \
    samtools merge -f -@ "$MERGE_THREADS" -b "$MERGELIST" "$ASSEMBLY/all.bam"
bench "samtools sort merged" \
    samtools sort  -@ "$MERGE_THREADS" -m "${MERGE_MEM}M" -T "$TMP/sort_all" \
    -o "$ASSEMBLY/all.sorted.bam" "$ASSEMBLY/all.bam"
bench "samtools index merged" \
    samtools index -@ "$MERGE_THREADS" "$ASSEMBLY/all.sorted.bam"
rm -f "$ASSEMBLY/all.bam"

log "stringtie: $THREADS threads ($ST_STRAND)"
bench "stringtie" stringtie "$ASSEMBLY/all.sorted.bam" \
    "$ST_STRAND" --conservative -p "$THREADS" \
    -o "$ASSEMBLY/assembly.gtf" \
    > "$LOGS/stringtie.log" 2>&1 || die "stringtie failed"

log "assembled transcripts: $(awk -F'\t' '$3=="transcript"' "$ASSEMBLY/assembly.gtf" | wc -l)"
step_end

############################################################
# STEP 10 - gffcompare rounds
############################################################

step_start "Identifying antisense transcripts (gffcompare)"

bench "gffcompare against targets" \
    gffcompare -r "$TMP/target_orphans.gff3" "$ASSEMBLY/assembly.gtf" \
    -o "$ASSEMBLY/orphan" > "$LOGS/gffcompare_1.log" 2>&1

awk -F'\t' '$4=="x"{split($5,a,"|"); print a[2]}' \
    "$ASSEMBLY/orphan.tracking" | sort -u > "$TMP/x1.txt"

awk '
NR==FNR{keep[$1]; next}
match($0,/transcript_id "[^"]+"/){
    id=substr($0,RSTART+15,RLENGTH-16)
    if(id in keep) print
}' "$TMP/x1.txt" "$ASSEMBLY/assembly.gtf" > "$TMP/nats_round1.gtf"

log "round 1 antisense transcripts: $(wc -l < "$TMP/x1.txt")"

bench "gffcompare against annotation" \
    gffcompare -r "$ANNOT" "$TMP/nats_round1.gtf" \
    -o "$ASSEMBLY/final" > "$LOGS/gffcompare_2.log" 2>&1

awk -F'\t' '$4=="x"{split($5,a,"|"); print a[2]}' \
    "$ASSEMBLY/final.tracking" | sort -u > "$TMP/x2.txt"

awk -F'\t' '$4=="x"' "$ASSEMBLY/final.tracking" > "$OUTDIR/NAT_tracking_x.tsv"

awk -F'\t' '$4=="x"{split($3,a,"|"); sub(/^gene:/,"",a[1]); print a[1]}' \
    "$ASSEMBLY/final.tracking" | sort -u > "$TMP/nat_target_genes.txt"

log "round 2 antisense transcripts: $(wc -l < "$TMP/x2.txt")"
step_end


############################################################
# STEP 11 - Final GTF
############################################################

step_start "Writing final NAT annotation"

if [[ "$MODE" == "test" ]]; then
    FINAL_GTF="$OUTDIR/FINAL_QUERY.gtf"
else
    FINAL_GTF="$OUTDIR/FINAL_NATS.gtf"
fi

awk '
NR==FNR{keep[$1]; next}
match($0,/transcript_id "[^"]+"/){
    id=substr($0,RSTART+15,RLENGTH-16)
    if(id in keep) print
}' "$TMP/x2.txt" "$TMP/nats_round1.gtf" > "$FINAL_GTF"

N_FINAL=$(awk -F'\t' '$3=="transcript"' "$FINAL_GTF" | wc -l)

log "final NAT transcripts: $N_FINAL"
log "final annotation: $FINAL_GTF"

step_end


############################################################
# STEP 12 - Test-mode query tracking and metrics
############################################################

if [[ "$MODE" == "test" ]]; then

    step_start "Test-mode query tracking and gffcompare metrics"

    # Inputs from Stage 1
    QUERY_BED="$STAGE1_DIR/TEST_KNOWNPAIRS_QUERY_expression_filtered.bed"
    QUERY_SOURCE_GFF="$STAGE1_DIR/annotation.gff3"
    QUERY_ANNOT="$TMP/annotation_query.gff3"

    # Output directory for this gffcompare run
    FINAL_TRACKINGS="$TMP/final_trackings"
    mkdir -p "$FINAL_TRACKINGS"

    [[ -s "$QUERY_BED" ]] \
        || die "test query BED not found or empty: $QUERY_BED"

    [[ -s "$QUERY_SOURCE_GFF" ]] \
        || die "original Stage 1 annotation not found or empty: $QUERY_SOURCE_GFF"

    [[ -s "$OUTDIR/FINAL_QUERY.gtf" ]] \
        || die "final query annotation not found or empty: $OUTDIR/FINAL_QUERY.gtf"

    # Extract the query gene IDs from column 4 of the BED.
    cut -f4 "$QUERY_BED" | sort -u > "$TMP/query_gene_ids.txt"

    [[ -s "$TMP/query_gene_ids.txt" ]] \
        || die "no query gene IDs found in $QUERY_BED"

    # Build annotation_query.gff3 from the ORIGINAL Stage 1 annotation.
    # Retain each selected gene and its associated child features.
    awk '
    BEGIN {
        FS = OFS = "\t"
    }

    function get_id(attr) {
        if (match(attr, /gene_id=[^;]+/))
            return substr(attr, RSTART + 8, RLENGTH - 8)

        if (match(attr, /geneID=[^;]+/))
            return substr(attr, RSTART + 7, RLENGTH - 7)

        if (match(attr, /ID=gene:[^;]+/))
            return substr(attr, RSTART + 8, RLENGTH - 8)

        if (match(attr, /ID=[^;]+/))
            return substr(attr, RSTART + 3, RLENGTH - 3)

        return ""
    }

    FNR == NR {
        query[$1] = 1
        next
    }

    /^#/ {
        print
        next
    }

    $3 == "gene" {
        keep = (get_id($9) in query)
    }

    keep {
        print
    }
    ' "$TMP/query_gene_ids.txt" "$QUERY_SOURCE_GFF" > "$QUERY_ANNOT"

    [[ -s "$QUERY_ANNOT" ]] \
        || die "failed to create query annotation: $QUERY_ANNOT"

    # Keep a clearly named copy of the final NAT models for comparison.
    FINAL_QUERIES="$FINAL_TRACKINGS/FINAL_QUERY.gtf"
    cp "$OUTDIR/FINAL_QUERY.gtf" "$FINAL_QUERIES"

    # Compare the final NATs against the query reference annotation.
    bench "gffcompare against queries" gffcompare \
        -r "$QUERY_ANNOT" \
        "$FINAL_QUERIES" \
        -o "$FINAL_TRACKINGS/final_queries" \
        > "$LOGS/gffcompare_final_queries.log" 2>&1 \
        || die "gffcompare against query annotation failed"

    # Copy the summary metrics to the main output directory.
    QUERY_STATS="$FINAL_TRACKINGS/final_queries.stats"

    [[ -s "$QUERY_STATS" ]] \
        || die "gffcompare completed but stats file is missing: $QUERY_STATS"

    cp "$QUERY_STATS" "$OUTDIR/final_queries.stats"

    log "query BED: $QUERY_BED"
    log "query reference annotation: $QUERY_ANNOT"
    log "final NATs compared: $FINAL_QUERIES"
    log "gffcompare tracking directory: $FINAL_TRACKINGS"
    log "query metrics: $OUTDIR/final_queries.stats"

    step_end
fi

############################################################
# Summary
############################################################

TOTAL=$(( $(date +%s) - PIPELINE_T0 ))


printf 'Results\n'
rule

printf '  mode             : %s\n' "$MODE"
printf '  final annotation : %s (%s transcripts)\n' \
       "$(basename "$FINAL_GTF")" "$N_FINAL"
printf '  tracking file    : %s\n' "$OUTDIR/NAT_tracking_x.tsv"

if [[ "$MODE" == "analysis" ]]; then
    printf '  orphan BED       : %s (%s entries)\n' \
           "$OUTDIR/orphans.bed" "$(wc -l < "$OUTDIR/orphans.bed")"
else
    printf '  selected test BED: %s (%s entries)\n' \
           "$OUTDIR/test_knownpairs.bed" \
           "$(wc -l < "$OUTDIR/test_knownpairs.bed")"
    printf '  query metrics    : %s\n' "$OUTDIR/final_queries.stats"
fi

printf '  output directory : %s\n' "$OUTDIR"
printf '  tool logs        : %s\n' "$LOGS"
printf '  elapsed time     : %s\n' "$(fmt_time "$TOTAL")"
rule
bench_stage_end
rule
printf 'Finished.\n'

