#!/usr/bin/env bash
############################################################
# STAGE 1: NASTIseq - Orphan Detection & FeatureCounts
#
# Usage:
#   bash ducktape_stage1.sh -g annotation.gff3 -b bam_dir -o stage1_results
#
# Run with --help for the full list of options.
############################################################

############################################################
# Defaults
############################################################

# Start time, for the stage resource usage summary
BENCH_T_START=${EPOCHREALTIME:-$(date +%s.%N)}

CONDA_ENV=""
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/scripts"

THREADS_OVERRIDE=""
MEM_GB_OVERRIDE=""
TEST_EXPRESSION="20"
QUERY_FRACTION="0.05"
BENCH_ON=1
BENCH_INTERVAL=1
ANNOT_FILES=()

############################################################
# Help
############################################################

usage() {
cat << EOF

============================================================
DuckTape NAT Discovery - STAGE 1: NASTIseq
============================================================

Required arguments
------------------

  -g, --annotation FILE
        Genome annotation. GTF or GFF3, optionally gzipped.
        May be given MULTIPLE times to merge annotations.

  -b, --bam-dir DIR
        Directory containing *.bam files.

  -o, --output DIR
        Output directory for Stage 1 results

Optional arguments
------------------

  --read-type PE|SE
        Paired-end or single-end reads
        Default: PE

  --library-type RF|FR
        Strandedness of the RNA-seq library.
        RF  dUTP / fr-firststrand  (Illumina TruSeq stranded)
        FR  fr-secondstrand
        Default: RF

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

  --test-expression FLOAT
        Minimum expression required for BOTH genes in a pair.
        Default: 20

  --query-fraction FLOAT
        Maximum query/test expression ratio for retaining a pair.
        If the test gene is the lower-expressed gene, the roles
        are swapped. Default: 0.05

  --no-benchmark
        Do not record resource usage (CPU, memory, I/O)

  --benchmark-interval SEC
        Memory sampling interval for the resource usage tables
        Default: 1

  -h, --help
        Show this help message and exit

Outputs
-------

  annotation.gff3                                 Merged/normalized annotation (GFF3)
  reference.gtf                                   Merged/normalized annotation (GTF)
  featurecounts/<sample>.sense.txt                Sense strand counts per sample
  featurecounts/<sample>.antisense.txt            Antisense strand counts per sample
  orphans.bed                                     NASTIseq orphan loci
  TEST_KNOWNPAIRS_expression_filtered.bed         Test genes of the known pairs
  TEST_KNOWNPAIRS_QUERY_expression_filtered.bed   Query genes of the known pairs
  reference_test.gff3                             annotation.gff3 without the query genes
  benchmark/steps.tsv                             Resource usage per step and for the stage
  benchmark/commands.tsv                          Resource usage per command

Examples
--------

  bash ducktape_stage1.sh \\
      -g Arabidopsis.gff3 \\
      -b hisat_bams \\
      -o stage1_results \\
      --threads 32 --mem-gb 200

============================================================

EOF
}

############################################################
# Parse args
############################################################

while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)          usage; exit 0 ;;
        -g|--annotation|--gff3)
                            ANNOT_FILES+=("$2"); shift 2 ;;
        -b|--bam-dir)       BAM_DIR="$2";        shift 2 ;;
        -o|--output)        OUTDIR="$2";         shift 2 ;;
        --library-type)     LIB_TYPE="${2^^}";   shift 2 ;;
        --read-type)        READ_TYPE="${2^^}";  shift 2 ;;
        --threads)          THREADS_OVERRIDE="$2"; shift 2 ;;
        --mem-gb)           MEM_GB_OVERRIDE="$2";  shift 2 ;;
        --env)              CONDA_ENV="$2";      shift 2 ;;
        --pipeline-dir)     PIPELINE_DIR="$2";   shift 2 ;;
        --test-expression)  TEST_EXPRESSION="$2"; shift 2 ;;
        --query-fraction)   QUERY_FRACTION="$2";  shift 2 ;;
        --no-benchmark)     BENCH_ON=0;           shift ;;
        --benchmark-interval)
                            BENCH_INTERVAL="$2";  shift 2 ;;
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

if [[ ${#ANNOT_FILES[@]} -eq 0 || -z "${BAM_DIR:-}" || -z "${OUTDIR:-}" ]]; then
    echo "ERROR: Missing required arguments"
    usage
    exit 1
fi

# Set defaults
LIB_TYPE="${LIB_TYPE:-RF}"
READ_TYPE="${READ_TYPE:-PE}"

case "$LIB_TYPE" in
    RF|FR) ;;
    *) echo "ERROR: --library-type must be RF or FR"; exit 1 ;;
esac

case "$READ_TYPE" in
    PE|SE) ;;
    *) echo "ERROR: --read-type must be PE or SE"; exit 1 ;;
esac

[[ "$TEST_EXPRESSION" =~ ^[0-9]+([.][0-9]+)?$ ]] \
    || die "--test-expression must be a non-negative number"

awk -v x="$QUERY_FRACTION" \
    'BEGIN { exit !(x ~ /^[0-9]+([.][0-9]+)?$/ && x >= 0 && x <= 1) }' \
    || die "--query-fraction must be a number between 0 and 1"

awk -v x="$BENCH_INTERVAL" \
    'BEGIN { exit !(x ~ /^[0-9]+([.][0-9]+)?$/ && x > 0) }' \
    || die "--benchmark-interval must be a positive number"

mkdir -p "$OUTDIR"

############################################################
# Environment setup
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

MEM_USABLE_MB=$(( MEM_MB * 80 / 100 ))

############################################################
# Adaptive scheduler
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

sched_mem_mb() { echo $(( MEM_USABLE_MB / SCHED_CONC )); }

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

for tool in gffread featureCounts samtools python Rscript xargs \
            agat agat_convert_sp_gxf2gxf.pl \
            agat_convert_sp_gff2gtf.pl \
            agat_sp_merge_annotations.pl; do
    command -v "$tool" >/dev/null 2>&1 \
        || die "'$tool' not found on PATH; activate the DuckTape conda environment or use --env"
done

SAMTOOLS_VER=$(samtools --version | head -1 | awk '{print $2}')
SAMTOOLS_MAJ=${SAMTOOLS_VER%%.*}
SAMTOOLS_MIN=${SAMTOOLS_VER#*.}; SAMTOOLS_MIN=${SAMTOOLS_MIN%%.*}
if (( SAMTOOLS_MAJ < 1 )) || { (( SAMTOOLS_MAJ == 1 )) && (( SAMTOOLS_MIN < 12 )); }; then
    die "samtools >= 1.12 required; found $SAMTOOLS_VER"
fi

############################################################
# Helper scripts
############################################################

NASTI_R="$PIPELINE_DIR/nastiseq_pipeline.R"
BENCH_LIB="$PIPELINE_DIR/benchmark.sh"

for f in "$NASTI_R" "$BENCH_LIB" "$PIPELINE_DIR/ducktape_bench.py"; do
    [[ -f "$f" ]] || die "helper script not found: $f
       Set --pipeline-dir to the directory containing it."
done

BENCH_STAGE="stage1"
source "$BENCH_LIB"


############################################################
# Validate inputs
############################################################

for a in "${ANNOT_FILES[@]}"; do
    [[ -f "$a" ]] || die "annotation file not found: $a"
done

[[ -d "$BAM_DIR" ]] || die "BAM directory not found: $BAM_DIR"

mapfile -t BAM_FILES < <(find "$BAM_DIR" -maxdepth 1 -name '*.bam' | sort)
N_SAMPLES=${#BAM_FILES[@]}
(( N_SAMPLES > 0 )) || die "no *.bam files found in $BAM_DIR"

declare -a SAMPLES=()
for bam in "${BAM_FILES[@]}"; do
    SAMPLES+=("$(basename "$bam" .bam)")
    [[ -f "${bam}.bai" || -f "${bam%.bam}.bai" ]] \
        || warn "no index found for $(basename "$bam")"
done

############################################################
# Library-type dependent settings
############################################################

if [[ "$LIB_TYPE" == "RF" ]]; then
    FC_SENSE=2; FC_ANTI=1
else
    FC_SENSE=1; FC_ANTI=2
fi

if [[ "$READ_TYPE" == "PE" ]]; then
    FC_PAIRED=(-p --countReadPairs -B -C)
else
    FC_PAIRED=()
fi

############################################################
# Working directories
############################################################

TMP="$OUTDIR/tmp_work"
rm -rf "$TMP"
mkdir -p "$TMP"

LOGS="$TMP/logs"
COUNTS="$OUTDIR/featurecounts"
NASTI="$TMP/nastiseq"

mkdir -p "$LOGS" "$COUNTS" "$NASTI/results"

# GNU sort tuning
SORT_MEM=$(( MEM_USABLE_MB / 4 ))
(( SORT_MEM > 8192 )) && SORT_MEM=8192
(( SORT_MEM < 256 ))  && SORT_MEM=256
SORT_OPTS=(--parallel="$THREADS" -S "${SORT_MEM}M" -T "$TMP")

############################################################
# Run header
############################################################

rule
printf 'DuckTape NAT Discovery - STAGE 1: NASTIseq\n'
rule
printf '  annotation      : %s\n' "${ANNOT_FILES[0]}"
for (( i=1; i<${#ANNOT_FILES[@]}; i++ )); do
    printf '                    %s\n' "${ANNOT_FILES[$i]}"
done
printf '  BAM directory   : %s\n' "$BAM_DIR"
printf '  samples         : %d\n' "$N_SAMPLES"
printf '  output          : %s\n' "$OUTDIR"
printf '  library type    : %s (%s)\n' "$LIB_TYPE" "$READ_TYPE"
printf '  CPUs            : %d\n' "$THREADS"
printf '  memory          : %d GB (%d GB usable)\n' \
       $(( MEM_MB / 1024 )) $(( MEM_USABLE_MB / 1024 ))
rule

############################################################
# Export for parallel workers
############################################################

export TMP LOGS COUNTS
export FC_SENSE FC_ANTI
export NASTI_R

bench_init

############################################################
# STEP 1 - Prepare annotation with AGAT + gffread
############################################################

step_start "Preparing annotation with AGAT + gffread"

# ----------------------------------------------------------
# 1. Normalize EVERY input annotation to GFF3 using AGAT
# ----------------------------------------------------------

AGAT_DIR="$TMP/AGAT"

mkdir -p "$AGAT_DIR"

AGAT_FILES=()
idx=0

for a in "${ANNOT_FILES[@]}"; do
    idx=$((idx + 1))

    # AGAT can read compressed annotation files, but decompress
    # explicitly to a temporary file for predictable handling.
    src="$a"

    if [[ "$a" == *.gz ]]; then
        src="$TMP/agat_input_${idx}.gff"
        log "Decompressing annotation: $(basename "$a")"
        bench "gzip ${a##*/}" gzip -cd "$a" > "$src" \
            || die "could not decompress annotation: $a"
    fi

    out="$AGAT_DIR/annotation_${idx}.gff3"

    log "AGAT preprocessing [$idx/${#ANNOT_FILES[@]}]: $(basename "$a")"

    bench "agat_convert_sp_gxf2gxf.pl ${a##*/}" agat_convert_sp_gxf2gxf.pl \
        --g "$src" \
        --output "$out" \
        --cpu "${THREADS}" \
        >"$LOGS/agat_${idx}.log" 2>&1 \
        || die "AGAT failed on $a; see $LOGS/agat_${idx}.log"

    [[ -s "$out" ]] \
        || die "AGAT produced an empty output for $a"

    AGAT_FILES+=("$out")
done


ANNOT="$TMP/annotation_merged.gff3"

if (( ${#AGAT_FILES[@]} == 1 )); then
    log "Only one annotation supplied; using normalized GFF3"
    cp "${AGAT_FILES[0]}" "$ANNOT"
else
    log "Merging ${#AGAT_FILES[@]} normalized annotations with AGAT"

    MERGE_ARGS=()
    for f in "${AGAT_FILES[@]}"; do
        MERGE_ARGS+=(--gff "$f")
    done

    bench "agat_sp_merge_annotations.pl" agat_sp_merge_annotations.pl \
        "${MERGE_ARGS[@]}" \
        --out "$ANNOT" \
        --cpu "${THREADS}" \
        >"$LOGS/agat_merge.log" 2>&1 \
        || die "AGAT annotation merge failed; see $LOGS/agat_merge.log"
fi

[[ -s "$ANNOT" ]] || die "merged annotation is empty"

# ----------------------------------------------------------
# 2. Create the final reference GFF3
# ----------------------------------------------------------


GTF="$OUTDIR/reference.gtf"

log "Converting normalized GFF3 to GTF"

bench "gffread GFF3 to GTF" gffread "$ANNOT" -T --keep-genes -o "$GTF" \
    >"$LOGS/gffread_gtf.log" 2>&1 \
    || die "gffread failed building reference GTF"

[[ -s "$GTF" ]] || die "gffread produced an empty GTF"

log "Converting GTF back to GFF3"

# Without -T, gffread writes GFF3 output.
bench "gffread GTF to GFF3" gffread "$GTF" --keep-genes -o "$ANNOT" \
    >"$LOGS/gffread_gff.log" 2>&1 \
    || die "gffread failed building reference GFF3"

[[ -s "$ANNOT" ]] || die "gffread produced an empty GFF3"

# Validate the GTF using GTF-style attributes.
awk -F '\t' '
    $0 !~ /^#/ && $3 == "exon" {
        total++
        if ($9 !~ /gene_id "/) missing_gene++
        if ($9 !~ /transcript_id "/) missing_tx++
    }
    END {
        printf "GTF exon records: %d\n", total+0
        printf "Exons missing gene_id: %d\n", missing_gene+0
        printf "Exons missing transcript_id: %d\n", missing_tx+0
        if (total == 0 || missing_gene > 0 || missing_tx > 0)
            exit 1
    }
' "$GTF" || die "GTF validation failed"

# Validate the GFF3 using GFF3-style attributes.
awk -F '\t' '
    $0 !~ /^#/ && NF != 9 { bad_columns++; next }
    $0 !~ /^#/ && NF == 9 {
        total++
        if ($3 == "gene" && $9 !~ /(^|;)ID=/) missing_id++
    }
    END {
        printf "GFF3 feature records: %d\n", total+0
        printf "Gene records missing ID: %d\n", missing_id+0
        if (total == 0 || bad_columns > 0 || missing_id > 0)
            exit 1
    }
' "$ANNOT" || die "GFF3 validation failed"

# Export the final annotation that downstream steps will actually use.
cp "$ANNOT" "$OUTDIR/annotation.gff3"

export ANNOT GTF


log "annotation ready: $(grep -vc '^#' "$ANNOT") GFF3 feature lines"
log "GFF3: $OUTDIR/annotation.gff3"
log "GTF:  $GTF"

step_end




############################################################
# STEP 2 - featureCounts
############################################################

step_start "Counting reads (featureCounts)"

schedule $(( N_SAMPLES * 2 )) 4 8
FC_THREADS=$SCHED_THREADS
export FC_THREADS FC_PAIRED_STR="${FC_PAIRED[*]:-}" BAM_DIR

log "$SCHED_CONC concurrent jobs x $FC_THREADS threads"

fc_one() {
    local sample="$1"
    local bam="$BAM_DIR/${sample}.bam"
    local -a paired=()

    [[ -n "${FC_PAIRED_STR:-}" ]] &&
        read -r -a paired <<< "$FC_PAIRED_STR"

    bench "featureCounts sense ${sample}" featureCounts -T "$FC_THREADS" -a "$GTF" \
        -o "$COUNTS/${sample}.sense.txt" \
        -g gene_id -t exon -s "$FC_SENSE" \
        "${paired[@]}" -O "$bam" \
        > "$LOGS/featureCounts_${sample}_sense.log" 2>&1

    bench "featureCounts antisense ${sample}" featureCounts -T "$FC_THREADS" -a "$GTF" \
        -o "$COUNTS/${sample}.antisense.txt" \
        -g gene_id -t exon -s "$FC_ANTI" \
        "${paired[@]}" -O "$bam" \
        > "$LOGS/featureCounts_${sample}_antisense.log" 2>&1
}
export -f fc_one

run_parallel fc_one "${SAMPLES[@]}" \
    || die "featureCounts failed; inspect $LOGS"

# Informational strandedness sanity check.
first="${SAMPLES[0]}"

pct() {
    grep -oE 'Successfully assigned alignments : [0-9]+ \([0-9.]+%\)' "$1" 2>/dev/null \
        | grep -oE '\([0-9.]+%' \
        | tr -d '(%' \
        | head -1 || true
}

S_PCT=$(pct "$LOGS/featureCounts_${first}_sense.log")
A_PCT=$(pct "$LOGS/featureCounts_${first}_antisense.log")
S_PCT=${S_PCT:-0}
A_PCT=${A_PCT:-0}

log "Assignment rate on $first: sense ${S_PCT}% / antisense ${A_PCT}%"

if awk -v s="$S_PCT" -v a="$A_PCT" 'BEGIN { exit !(a > s) }'; then
    warn "Antisense assignment exceeds sense assignment; verify --library-type $LIB_TYPE"
fi

step_end


############################################################
# STEP 3 - NASTIseq and paired BED creation
############################################################

step_start "NASTIseq known-pairs creation"

bench "NASTIseq" Rscript "$NASTI_R" \
    -i "$COUNTS" \
    -g "$GTF" \
    -c "$THREADS" \
    -v 0.2 \
    -o "$NASTI/results" \
    || die "NASTIseq failed"

KNOWNPAIRS=$(find "$NASTI/results" -name '*_knownpairs.bed' -print -quit)
[[ -n "$KNOWNPAIRS" && -s "$KNOWNPAIRS" ]] \
    || die "NASTIseq did not produce a non-empty knownpairs BED"

# This split assumes NASTIseq writes each pair as two consecutive rows:
# odd-numbered rows are query genes; even-numbered rows are test genes.
# Confirm this ordering matches the output format of your NASTIseq version.
PAIR_TEST_ORIGINAL="$TMP/TEST_KNOWNPAIRS.bed"
PAIR_QUERY_ORIGINAL="$TMP/TEST_KNOWNPAIRS_QUERY.bed"

awk 'NR % 2 == 0' "$KNOWNPAIRS" > "$PAIR_TEST_ORIGINAL"
awk 'NR % 2 == 1' "$KNOWNPAIRS" > "$PAIR_QUERY_ORIGINAL"

[[ -s "$PAIR_TEST_ORIGINAL" ]] \
    || die "test known-pairs BED is empty"
[[ -s "$PAIR_QUERY_ORIGINAL" ]] \
    || die "query known-pairs BED is empty"

N_TEST=$(wc -l < "$PAIR_TEST_ORIGINAL")
N_QUERY=$(wc -l < "$PAIR_QUERY_ORIGINAL")


ORPHANS_SOURCE="$NASTI/results/reference_orphans.bed"
ORPHANS_OUT="$OUTDIR/orphans.bed"

[[ -s "$ORPHANS_SOURCE" ]] \
    || die "NASTIseq reference orphan BED not found or empty: $ORPHANS_SOURCE"

cp "$ORPHANS_SOURCE" "$ORPHANS_OUT"
[[ -s "$ORPHANS_OUT" ]] \
    || die "failed to copy reference orphan BED to: $ORPHANS_OUT"

log "Copied NASTIseq reference orphans: $ORPHANS_SOURCE -> $ORPHANS_OUT"
log "Original test rows:  $N_TEST"
log "Original query rows: $N_QUERY"

step_end


############################################################
# STEP 4 - Expression-filtered known-pairs BED files
############################################################

step_start "Filtering known pairs by expression"

EXPR_TABLE="$TMP/knownpairs_gene_expression.tsv"

# Sum sense counts across all samples.
# This makes TEST_EXPRESSION a cohort-wide total, not a per-sample threshold.
awk '
BEGIN { FS=OFS="\t" }
FNR > 2 {
    id=$1
    sub(/^gene:/, "", id)
    count[id] += $7
}
END {
    for (id in count) print id, count[id]
}' "$COUNTS"/*.sense.txt > "$EXPR_TABLE"

[[ -s "$EXPR_TABLE" ]] || die "expression table is empty"

TEST_KNOWNPAIRS="$OUTDIR/TEST_KNOWNPAIRS_expression_filtered.bed"
TEST_KNOWNPAIRS_QUERY="$OUTDIR/TEST_KNOWNPAIRS_QUERY_expression_filtered.bed"

bench "known_pairs_filter (inline python)" python - \
    "$EXPR_TABLE" \
    "$PAIR_TEST_ORIGINAL" \
    "$PAIR_QUERY_ORIGINAL" \
    "$TEST_KNOWNPAIRS" \
    "$TEST_KNOWNPAIRS_QUERY" \
    "$TEST_EXPRESSION" \
    "$QUERY_FRACTION" <<'PYCODE'
import sys

(
    expr_file,
    test_file,
    query_file,
    out_test,
    out_query,
    min_expression_arg,
    query_fraction_arg,
) = sys.argv[1:]

min_expression = float(min_expression_arg)
query_fraction = float(query_fraction_arg)


def norm_id(value):
    value = value.strip()
    return value[5:] if value.startswith("gene:") else value


expression = {}
with open(expr_file) as handle:
    for line in handle:
        fields = line.rstrip("\n").split("\t")
        if len(fields) >= 2:
            expression[norm_id(fields[0])] = float(fields[1])


def read_bed(path):
    with open(path) as handle:
        return [
            line for line in handle
            if line.strip() and not line.startswith("#")
        ]


test_rows = read_bed(test_file)
query_rows = read_bed(query_file)


kept = 0
swapped = 0
excluded = 0

with open(out_test, "w") as test_out, open(out_query, "w") as query_out:
    for line_test, line_query in zip(test_rows, query_rows):
        fields_test = line_test.rstrip("\n").split("\t")
        fields_query = line_query.rstrip("\n").split("\t")

        if len(fields_test) < 4 or len(fields_query) < 4:
            raise SystemExit("ERROR: known-pairs BED row has fewer than 4 columns")

        id_test = norm_id(fields_test[3])
        id_query = norm_id(fields_query[3])

        missing = [
            gene_id for gene_id in (id_test, id_query)
            if gene_id not in expression
        ]
        if missing:
            raise SystemExit(
                "ERROR: no featureCounts sense-count row for gene(s): "
                + ", ".join(missing)
                + ". Check the BED IDs against featureCounts gene IDs."
            )

        expr_test = expression[id_test]
        expr_query = expression[id_query]

        # Both genes must meet the minimum total expression threshold.
        if expr_test < min_expression or expr_query < min_expression:
            excluded += 1
            continue

        # Keep original orientation if query is sufficiently lower expressed.
        if expr_query <= query_fraction * expr_test:
            test_out.write(line_test)
            query_out.write(line_query)
            kept += 1

        # Swap roles if the original test gene is sufficiently lower expressed.
        elif expr_test <= query_fraction * expr_query:
            test_out.write(line_query)
            query_out.write(line_test)
            swapped += 1

        else:
            excluded += 1

print(f"Minimum expression threshold: {min_expression:g}")
print(f"Maximum query/test expression fraction: {query_fraction:g}")
print(f"Pairs retained (query lower): {kept}")
print(f"Pairs retained (roles swapped): {swapped}")
print(f"Pairs excluded: {excluded}")
PYCODE

[[ -s "$TEST_KNOWNPAIRS" && -s "$TEST_KNOWNPAIRS_QUERY" ]] \
    || die "no known pairs passed expression filtering"

N_FILTERED_TEST=$(wc -l < "$TEST_KNOWNPAIRS")
N_FILTERED_QUERY=$(wc -l < "$TEST_KNOWNPAIRS_QUERY")


log "Expression-filtered test pairs:  $N_FILTERED_TEST"
log "Expression-filtered query pairs: $N_FILTERED_QUERY"

step_end


############################################################
# STEP 5 - Create reference_test.gff3
############################################################

step_start "Creating reference_test.gff3"

REFERENCE_TEST="$OUTDIR/reference_test.gff3"

bench "reference_test_builder (inline python)" python - \
    "$ANNOT" \
    "$TEST_KNOWNPAIRS_QUERY" \
    "$REFERENCE_TEST" <<'PYCODE'
import sys
from collections import defaultdict

annotation_path, query_bed_path, output_path = sys.argv[1:]


def norm_id(value):
    value = value.strip()
    return value[5:] if value.startswith("gene:") else value


def parse_attributes(text):
    attrs = {}
    for item in text.strip().split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value
    return attrs


# Read query IDs from column 4 of the expression-filtered query BED.
query_ids = set()
with open(query_bed_path) as handle:
    for line in handle:
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 4:
            raise SystemExit("ERROR: query BED contains a row with fewer than 4 columns")
        query_ids.add(norm_id(fields[3]))

if not query_ids:
    raise SystemExit("ERROR: no query gene IDs found in filtered query BED")

headers = []
features = []
parent_map = defaultdict(list)
gene_ids = set()

with open(annotation_path) as handle:
    for line in handle:
        if line.startswith("##FASTA"):
            # Do not copy embedded FASTA into a feature-only subset GFF3.
            break

        if line.startswith("#") or not line.strip():
            headers.append(line)
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9:
            raise SystemExit(
                f"ERROR: malformed GFF3 record with {len(fields)} columns"
            )

        attrs = parse_attributes(fields[8])
        feature_id = norm_id(attrs.get("ID", ""))
        parents = [
            norm_id(parent)
            for parent in attrs.get("Parent", "").split(",")
            if parent
        ]

        record = {
            "line": line,
            "type": fields[2],
            "id": feature_id,
            "parents": parents,
        }
        features.append(record)

        if fields[2] == "gene" and feature_id:
            gene_ids.add(feature_id)

        if feature_id:
            parent_map[feature_id].extend(parents)

# Validate exact ID matches before removing anything.
matched = query_ids.intersection(gene_ids)
unmatched = query_ids - gene_ids

if unmatched:
    preview = ", ".join(sorted(unmatched)[:10])
    raise SystemExit(
        f"ERROR: {len(unmatched)} query ID(s) did not match gene IDs "
        f"in the annotation. Examples: {preview}. "
        "Check ID prefixes and the BED/GFF3 ID conventions."
    )

# Find each feature's own ID and all its ancestors.
def ancestry(record):
    seen = set()
    stack = []

    if record["id"]:
        stack.append(record["id"])
    stack.extend(record["parents"])

    while stack:
        current = stack.pop()
        if not current or current in seen:
            continue
        seen.add(current)
        stack.extend(parent_map.get(current, []))

    return seen


removed_records = 0

with open(output_path, "w") as out:
    out.writelines(headers)

    for record in features:
        if ancestry(record).intersection(query_ids):
            removed_records += 1
            continue
        out.write(record["line"])

print(f"Query gene IDs supplied: {len(query_ids)}")
print(f"Query gene IDs matched:  {len(matched)}")
print(f"Annotation records removed, including descendants: {removed_records}")
print(f"Created: {output_path}")
PYCODE

[[ -s "$REFERENCE_TEST" ]] \
    || die "reference_test.gff3 was not created or is empty"

log "Created reference: $REFERENCE_TEST"
log "Remaining GFF3 feature lines: $(grep -vc '^#' "$REFERENCE_TEST")"

step_end

rule
log "Stage 1 processing complete"
log "Output directory: $OUTDIR"
log "Main outputs:"
log "  reference.gtf"
log "  annotation.gff3"
log "  TEST_KNOWNPAIRS_expression_filtered.bed"
log "  TEST_KNOWNPAIRS_QUERY_expression_filtered.bed"
log "  reference_test.gff3"
rule

bench_stage_end
