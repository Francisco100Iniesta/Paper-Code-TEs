#!/usr/bin/env bash
#
###############################################################################
# Regional DNA methylation profiling around TE insertion sites
#
# PURPOSE
# -------
# This script calculates DNA methylation across predefined genomic regions
# surrounding non-reference transposable element (TE) insertion sites.
#
# In this study, the analyzed regions corresponded to 500-bp windows located
# immediately upstream (left flank) and downstream (right flank) of TE insertion
# breakpoints identified by RetroInspector.
#
# The script uses a haplotagged nanopore BAM and runs Modkit with `--phased`.
# This allows methylation to be evaluated both at the regional level and with
# haplotype resolution (HP1 and HP2), enabling downstream comparison of the
# TE-containing and insertion-free haplotypes.
#
# IMPORTANT
# ---------
# This script does NOT generate the flanking regions. The left and right BED
# files must already contain the genomic windows to be analyzed.
#
# For the analyses reported in the paper:
#
#   LEFT flank  = 500 bp upstream of the RetroInspector insertion breakpoint
#   RIGHT flank = 500 bp downstream of the RetroInspector insertion breakpoint
#
# INPUTS
# ------
# 1. Haplotagged nanopore BAM file
#      - Reads should contain haplotype assignments (HP tags).
#      - The BAM should also retain nanopore base-modification information.
#
# 2. GRCh38 reference FASTA
#
# 3. Left-flank BED file
#      chromosome    start    end
#
# 4. Right-flank BED file
#      chromosome    start    end
#
# OUTPUT
# ------
# Modkit pileup output for the left and right flanking regions.
#
# Because `--phased` is enabled, the resulting methylation output can be used
# for both:
#
#   - combined regional methylation
#   - haplotype-resolved methylation (HP1 / HP2)
#
# These outputs are subsequently used to compare methylation between samples
# and, for heterozygous TE insertions, between the TE-containing haplotype and
# the insertion-free haplotype within the same sample.
#
# SOFTWARE REQUIREMENTS
# ---------------------
# bash
# modkit
#
# The exact software environment used for the paper should be provided in the
# repository as `environment.yml` (and optionally an explicit lock file).
#
###############################################################################

###############################################################################
# SLURM SETTINGS
#
# Adapt these directives to your computing environment or remove them when
# running outside a SLURM cluster.
###############################################################################

#SBATCH -p nadal-q
#SBATCH --cpus-per-task=36
#SBATCH --mem=100G
#SBATCH --mail-type=END
#SBATCH -J regional_methylation
#SBATCH -o regional_methylation_%j.out


set -Eeuo pipefail

trap 'echo "[ERROR] Failure at line $LINENO: $BASH_COMMAND" >&2' ERR


###############################################################################
# ENVIRONMENT
###############################################################################

source ~/.bashrc
mamba activate nanopore

CORES="${SLURM_CPUS_PER_TASK:-1}"


###############################################################################
# CONFIGURATION
#
# Replace the paths below with those corresponding to your analysis.
###############################################################################

BAM="/path/to/haplotagged_sample.bam"

REFERENCE="/path/to/GRCh38.fa"

LEFT_BED="/path/to/left_flank.bed"
RIGHT_BED="/path/to/right_flank.bed"

SAMPLE_NAME="sample_name"

OUTDIR="./regional_methylation_results"


###############################################################################
# INITIAL CHECKS
###############################################################################

mkdir -p "$OUTDIR"

test -s "$BAM"
test -s "$REFERENCE"
test -s "$LEFT_BED"
test -s "$RIGHT_BED"

echo "[INFO] Starting regional methylation analysis"
echo "[INFO] BAM: $BAM"
echo "[INFO] Reference: $REFERENCE"
echo "[INFO] Left-flank BED: $LEFT_BED"
echo "[INFO] Right-flank BED: $RIGHT_BED"
echo "[INFO] Output directory: $OUTDIR"
echo "[INFO] Threads: $CORES"


###############################################################################
# 1) LEFT-FLANK METHYLATION
###############################################################################

echo "[INFO] Running Modkit pileup for left flanking regions"

modkit pileup \
    -t "$CORES" \
    --modified-bases 5mC 5hmC \
    --reference "$REFERENCE" \
    --include-bed "$LEFT_BED" \
    --phased \
    --prefix "${OUTDIR}/${SAMPLE_NAME}_left" \
    "$BAM" \
    "${OUTDIR}/${SAMPLE_NAME}_left"


###############################################################################
# 2) RIGHT-FLANK METHYLATION
###############################################################################

echo "[INFO] Running Modkit pileup for right flanking regions"

modkit pileup \
    -t "$CORES" \
    --modified-bases 5mC 5hmC \
    --reference "$REFERENCE" \
    --include-bed "$RIGHT_BED" \
    --phased \
    --prefix "${OUTDIR}/${SAMPLE_NAME}_right" \
    "$BAM" \
    "${OUTDIR}/${SAMPLE_NAME}_right"


###############################################################################
# FINAL
###############################################################################

echo "[INFO] Regional methylation analysis completed successfully"
echo "[INFO] Left-flank output prefix: ${OUTDIR}/${SAMPLE_NAME}_left"
echo "[INFO] Right-flank output prefix: ${OUTDIR}/${SAMPLE_NAME}_right"
