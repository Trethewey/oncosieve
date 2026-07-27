#!/usr/bin/env bash
# =============================================================================
# run_oncosieve.sh
#
# Orchestrates the full pan-cancer whitelist pipeline.
#
# Steps:
#   1. Check dependencies
#   2. Install Python requirements
#   3. Build the whitelist (all parsers -> merge -> filter -> tier -> output)
#   4. Index output VCF with tabix
#   5. Post-process whitelist (add genome_version, transcript_id, protein_change)
#   6. Print summary statistics
#   7. Post-pipeline (REVEL annotation, high-confidence filter, Excel, HTML report)
#   8. (Optional) Run Mutect2 post-filter rescue
#
# Usage:
#   bash run_oncosieve.sh [data_dir] [options]
#
# Arguments:
#   data_dir               Optional path to reference data directory.
#                          Overrides relative data paths in config.yaml.
#                          e.g. bash run_oncosieve.sh /path/to/reference/
#
# Options:
#   --skip-sources STR     Comma-separated sources to skip, e.g. "genie,cbioportal"
#   --from-intermediates   Skip parsers; load saved intermediates and re-run merge/filter/output
#   --intermediate-only    Run parsers and save intermediates; skip merge
#   --rescue               Run Mutect2 post-filter rescue after building whitelist
#   --mutect2-vcf PATH     Mutect2 FilterMutectCalls VCF (required if --rescue)
#   --tumor-sample NAME    Tumour sample name in the Mutect2 VCF
#   --rescue-output PATH   Output path for rescued VCF
#   --help                 Show this message and exit
# =============================================================================

set -euo pipefail

# Resolve DATA_DIR to absolute path BEFORE cd changes working directory
DATA_DIR=""
if [[ $# -gt 0 && ! "$1" =~ ^-- ]]; then
    DATA_DIR="$(realpath "$1")"
    shift
fi

# Always run relative to the script's own directory
cd "$(dirname "$(realpath "$0")")"

# =============================================================================
# DATA DIRECTORY — optional first positional argument
# =============================================================================

if [[ -n "$DATA_DIR" ]]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Data directory: $DATA_DIR"
fi

# =============================================================================
# USER CONFIGURATION — set paths here before running
# =============================================================================

PYTHON="python3"
CONFIG="config.yaml"

OUTPUT_DIR="output"
INTERMEDIATE_DIR="intermediate"
PREFIX="pan_cancer_whitelist_GRCh38"

REFERENCE_FASTA="data/reference/GRCh38.fa"

# =============================================================================
# ARGUMENT PARSING
# =============================================================================

SKIP_SOURCES=""
INTERMEDIATE_ONLY=false
FROM_INTERMEDIATES=false
DO_RESCUE=false
MUTECT2_VCF=""
TUMOR_SAMPLE=""
RESCUE_OUTPUT=""

usage() {
    grep '^#' "$0" | grep -v '#!/' | sed 's/^# \{0,3\}//'
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --skip-sources)      SKIP_SOURCES="$2";       shift 2 ;;
        --intermediate-only) INTERMEDIATE_ONLY=true;  shift   ;;
        --from-intermediates) FROM_INTERMEDIATES=true; shift   ;;
        --rescue)            DO_RESCUE=true;           shift   ;;
        --mutect2-vcf)       MUTECT2_VCF="$2";        shift 2 ;;
        --tumor-sample)      TUMOR_SAMPLE="$2";        shift 2 ;;
        --rescue-output)     RESCUE_OUTPUT="$2";       shift 2 ;;
        --help|-h)           usage ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

# =============================================================================
# HELPERS
# =============================================================================

log()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
die()  { echo "ERROR: $*" >&2; exit 1; }
warn() { echo "WARNING: $*" >&2; }

# =============================================================================
# LOGGING
# =============================================================================

LOG_DIR="logs"
mkdir -p "$LOG_DIR"
LOG_FILE="${LOG_DIR}/run_$(date '+%Y%m%d_%H%M%S').log"
exec > >(tee -a "$LOG_FILE") 2>&1
log "Log file: $LOG_FILE"

# =============================================================================
# BANNER
# =============================================================================

echo ""
echo "  ╔══════════════════════════════════════════════════════════════╗"
echo "  ║                                                              ║"
echo "  ║                                                              ║"
echo "  ║                      O N C O S I E V E                       ║"
echo "  ║                                                              ║"
echo "  ╚══════════════════════════════════════════════════════════════╝"
echo "         Pan-Cancer Somatic Variant Whitelist Curation Tool"
echo "                  github.com/Trethewey/oncosieve "
echo ""
log "Running pre-run audit..."
echo ""
"$PYTHON" pre_check.py --config "$CONFIG" ${DATA_DIR:+--data-dir "$DATA_DIR"} ${SKIP_SOURCES:+--skip-sources "$SKIP_SOURCES"} || die "Pre-run audit failed. Fix errors above before running pipeline."

# =============================================================================
# STEP 1b: CHECK DEPENDENCIES
# =============================================================================

log "Checking dependencies..."

"$PYTHON" --version  || die "python3 not found. Install Python >= 3.10."
tabix --version 2>&1 | head -1 || die "tabix not found. Install htslib."
bcftools --version | head -1   || die "bcftools not found."

BCFTOOLS_VERSION=$(bcftools --version | head -1 | grep -oP '\d+\.\d+' | head -1)

log "  Python   : $("$PYTHON" --version)"
log "  bcftools : $BCFTOOLS_VERSION"
log "  tabix    : $(tabix --version 2>&1 | head -1)"

# =============================================================================
# STEP 2: INSTALL PYTHON REQUIREMENTS
# =============================================================================

log "Installing Python requirements..."
"$PYTHON" -m pip install --quiet -r requirements.txt \
    || die "pip install failed. Check requirements.txt."
log "  Done."

# =============================================================================
# STEP 3: BUILD WHITELIST
# =============================================================================

WL_RAW="${OUTPUT_DIR}/${PREFIX}.vcf.gz"

if [[ -f "$WL_RAW" ]] && ! $FROM_INTERMEDIATES; then
    log "Whitelist VCF already exists at $WL_RAW — skipping build step."
    log "  Delete $WL_RAW to force a full rebuild."
else
    log "Building pan-cancer whitelist..."

    BUILD_ARGS="--config config.yaml${DATA_DIR:+ --data-dir $DATA_DIR}"
    [[ -n "$SKIP_SOURCES" ]] && BUILD_ARGS="$BUILD_ARGS --skip-sources $SKIP_SOURCES"
    $INTERMEDIATE_ONLY       && BUILD_ARGS="$BUILD_ARGS --intermediate-only"
    $FROM_INTERMEDIATES      && BUILD_ARGS="$BUILD_ARGS --from-intermediates"

    "$PYTHON" build_whitelist.py $BUILD_ARGS \
        || die "build_whitelist.py failed."
fi

if $INTERMEDIATE_ONLY; then
    log "Intermediate-only mode: stopping before merge."
    log "Intermediate files written to: $INTERMEDIATE_DIR/"
    exit 0
fi

WL_RAW="${OUTPUT_DIR}/${PREFIX}.vcf.gz"
[[ -f "$WL_RAW" ]] || die "Expected output VCF not found: $WL_RAW"

# =============================================================================
# STEP 4: INDEX VCF
# =============================================================================

WL_FINAL="$WL_RAW"
log "Indexing whitelist VCF..."
tabix -f -p vcf "$WL_FINAL" || die "tabix indexing failed."
log "  Indexed: ${WL_FINAL}.tbi"

# =============================================================================
# STEP 5: POST-PROCESS WHITELIST (add genome_version, transcript_id, protein_change)
# =============================================================================

WL_TSV="${OUTPUT_DIR}/${PREFIX}.tsv.gz"
WL_TSV_ANNOTATED="${OUTPUT_DIR}/${PREFIX}_annotated.tsv.gz"

log "Running post-processing..."
"$PYTHON" post_process_whitelist.py \
    --whitelist "$WL_TSV" \
    --out       "$WL_TSV_ANNOTATED" \
    || die "post_process_whitelist.py failed."
log "  Annotated TSV: $WL_TSV_ANNOTATED"

# =============================================================================
# STEP 6: SUMMARY STATISTICS
# =============================================================================

log "Summary statistics:"
N_TOTAL=$(bcftools view -H "$WL_FINAL" | wc -l | tr -d ' ')
N_TIER1=$(bcftools view -H -i 'INFO/WL_TIER=1' "$WL_FINAL" | wc -l | tr -d ' ')
N_TIER2=$(bcftools view -H -i 'INFO/WL_TIER=2' "$WL_FINAL" | wc -l | tr -d ' ')
N_TIER3=$(bcftools view -H -i 'INFO/WL_TIER=3' "$WL_FINAL" | wc -l | tr -d ' ')

log "  Total whitelist variants : $N_TOTAL"
log "  Tier 1 (highest conf)    : $N_TIER1"
log "  Tier 2                   : $N_TIER2"
log "  Tier 3 (min threshold)   : $N_TIER3"
log ""
log "  TSV (annotated) : $WL_TSV_ANNOTATED"
log "  TSV (base)      : $WL_TSV"
log "  VCF             : ${WL_FINAL}"

# =============================================================================
# STEP 7: POST-PIPELINE (REVEL, high-confidence, Excel, HTML report)
# =============================================================================

REVEL_FILE="${DATA_DIR:-data}/REVEL/revel_with_transcript_ids"
PRIMATEAI_FILE="${DATA_DIR:-data}/PRIMATE_AI/PrimateAI-3D.hg38.txt.gz"

if [[ -f "$REVEL_FILE" ]]; then
    log "Running post-pipeline processing..."

    PAI_ARG=""
    if [[ -f "$PRIMATEAI_FILE" ]]; then
        PAI_ARG="--primateai $PRIMATEAI_FILE"
        log "  PrimateAI-3D: $PRIMATEAI_FILE"
    else
        log "  PrimateAI-3D file not found — skipping PrimateAI-3D annotation."
    fi

    "$PYTHON" tools/post_pipeline.py \
        --whitelist "$WL_TSV_ANNOTATED" \
        --vcf       "$WL_FINAL" \
        --revel     "$REVEL_FILE" \
        $PAI_ARG \
        --out-dir   "$OUTPUT_DIR" \
        || warn "post_pipeline.py failed (non-fatal)."
else
    warn "REVEL file not found at $REVEL_FILE — skipping post-pipeline."
    warn "Download REVEL v1.3 to enable REVEL annotation, Excel export, and HTML report."
fi

# =============================================================================
# STEP 8: MUTECT2 POST-FILTER RESCUE (optional)
# =============================================================================

if $DO_RESCUE; then
    [[ -n "$MUTECT2_VCF" ]] || die "--rescue requires --mutect2-vcf PATH"
    [[ -f "$MUTECT2_VCF"  ]] || die "Mutect2 VCF not found: $MUTECT2_VCF"

    if [[ -z "$RESCUE_OUTPUT" ]]; then
        RESCUE_OUTPUT="${MUTECT2_VCF%.vcf*}.rescued.vcf.gz"
    fi

    TUMOR_ARG=""
    [[ -n "$TUMOR_SAMPLE" ]] && TUMOR_ARG="--tumor-sample $TUMOR_SAMPLE"

    log "Running Mutect2 post-filter rescue..."
    log "  Input  : $MUTECT2_VCF"
    log "  Output : $RESCUE_OUTPUT"

    "$PYTHON" mutect2_rescue.py \
        --input     "$MUTECT2_VCF" \
        --whitelist "$WL_FINAL" \
        --output    "$RESCUE_OUTPUT" \
        --config    config.yaml \
        $TUMOR_ARG \
        || die "mutect2_rescue.py failed."

    tabix -f -p vcf "$RESCUE_OUTPUT" \
        || warn "tabix indexing of rescued VCF failed."

    log "Rescue complete: $RESCUE_OUTPUT"
fi

log "Pipeline complete."
