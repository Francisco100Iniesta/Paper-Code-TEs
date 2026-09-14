#!/bin/bash
#
###############################################################################
# Internal methylation profiling of non-reference transposable element insertions
#
# PURPOSE
# -------
# This script quantifies DNA methylation within non-reference transposable
# element (TE) insertions detected with RetroInspector.
#
# Because these insertions are absent from the standard GRCh38 reference,
# methylation cannot be profiled directly across the inserted sequence using
# reference-based coordinates. For each TE insertion, this workflow therefore:
#
#   1. Reads the insertion breakpoint and inserted sequence reported by
#      RetroInspector.
#   2. Extracts genomic sequence flanking the insertion site from GRCh38.
#   3. Builds an insertion-specific custom contig:
#
#          [left reference flank] + [TE insertion] + [right reference flank]
#
#   4. Combines all insertion-specific contigs into a custom reference FASTA.
#   5. Extracts nanopore reads covering the corresponding regions from the
#      original BAM.
#   6. Converts the selected reads to FASTQ while retaining auxiliary tags.
#   7. Realigns the reads against the custom reference with minimap2.
#   8. Runs Modkit on the coordinates corresponding to the inserted TE sequence.
#   9. Produces BedMethyl output and regional methylation statistics.
#
# INPUTS
# ------
# 1. TE insertion table (`bedsecuencias`), tab-separated, with at least:
#
#      column 1: insertion ID
#      column 2: chromosome
#      column 3: insertion breakpoint
#      column 4: inserted TE sequence
#
#    The insertion breakpoint convention is controlled by INPUT_START_IS_BED0:
#      0 = column 3 is already 1-based (faidx/IGV-like coordinate)
#      1 = column 3 is BED-style 0-based and is converted to 1-based
#
# 2. GRCh38 reference FASTA (`reference`), with or without an existing .fai.
#
# 3. Nanopore BAM (`sample`) containing base-modification information.
#
# OUTPUTS
# -------
# - A custom FASTA containing one contig per TE insertion.
# - A BAM containing reads realigned to the custom reference.
# - BED files defining the original and custom-reference regions.
# - A table mapping original GRCh38 coordinates to custom-reference coordinates.
# - Modkit BedMethyl output.
# - Modkit regional methylation statistics for each inserted TE.
#
# SOFTWARE REQUIREMENTS
# ---------------------
# bash
# samtools
# minimap2
# seqkit
# modkit
# bgzip
# tabix
#
# IMPORTANT IMPLEMENTATION NOTES
# ------------------------------
# - `window=100000` defines the amount of reference sequence extracted on each
#   side of the insertion to provide alignment context. It is NOT the 500-bp
#   flanking window used elsewhere in the manuscript for flanking methylation
#   analyses.
#
# - `samtools fastq -T '*'` and `minimap2 -y` are used so auxiliary read tags
#   can be carried through the FASTQ realignment step. This is important when
#   working with nanopore base-modification tags.
#
# - Internal TE methylation is evaluated with a minimum coverage of 10 valid
#   modification calls per cytosine in the Modkit statistics step.
#
# - The SLURM directives and paths below reflect the original computing
#   environment used for the analysis and should be adapted when running on
#   another system.
#
# REPRODUCIBILITY
# ---------------
# For publication, store the software environment used for this analysis in the
# repository (for example as `environment.yml`) together with this script.
#
###############################################################################
#SBATCH -p nadal-q
#SBATCH --chdir=/home/finiesta/resultados_retroinspector/Modkit_Internal
#SBATCH -J internal_TE_methylation
#SBATCH --cpus-per-task=32
#SBATCH --mem=125G
#SBATCH -o /home/finiesta/resultados_retroinspector/Modkit_Internal/salida_ejecucion/%j.out

set -Ee -o pipefail

trap 'echo "[ERROR] Failure at line $LINENO: $BASH_COMMAND" >&2' ERR

source ~/.bashrc
mamba activate samtools_env

set -u

###############################################################################
# CONFIGURATION
###############################################################################

base_dir="/home/finiesta/resultados_retroinspector/Modkit_Internal"

bedsecuencias="${base_dir}/300_custom.bed"
reference="/pool_sas/crh/ref/hg38.fa"
sample="/pool_sas/crh/analysis/LSA-2026-0003/wgs/LSA-2026-0003_OM300-buffy/minimap2_clair3/alignments_LSA-2026-0003_OM300-buffy_sorted.bam"

sample_name="$(basename "$sample" .bam)"
outdir="${base_dir}/${sample_name}"

window=100000

# IMPORTANT:
# If column 3 of 300_custom.bed uses true BED coordinates (0-based), set this to 1.
# If column 3 is already a faidx/IGV-style coordinate (1-based), keep this at 0.
INPUT_START_IS_BED0=0

CORES=${SLURM_CPUS_PER_TASK:-1}

# Threads assigned to each stage of the workflow.
VIEW_THREADS=4
FASTQ_THREADS=4
SORT_THREADS=8
INDEX_THREADS=8

MAP_THREADS=$(( CORES - VIEW_THREADS - FASTQ_THREADS - SORT_THREADS ))

if (( MAP_THREADS < 1 )); then
    MAP_THREADS=1
fi

SORT_MEM="6G"

###############################################################################
# DIRECTORIES
###############################################################################

contigs_dir="${outdir}/contigs"
allcontigs_dir="${contigs_dir}/allcontigs"
customref_dir="${contigs_dir}/custom_reference"
bamcustom_dir="${outdir}/bamcustom"
modkit="${outdir}/Modkit"

mkdir -p \
    "$contigs_dir" \
    "$allcontigs_dir" \
    "$customref_dir" \
    "$bamcustom_dir" \
    "$modkit"

###############################################################################
# OUTPUT FILES
###############################################################################

custom_fa="${customref_dir}/custom.fa"

samtools_bed="${bamcustom_dir}/samtools.bed"
modkit_bed="${modkit}/modkit.bed"

regions_map="${modkit}/regions_map_${sample_name}.tsv"
original_regions_bed="${modkit}/original_regions_${sample_name}.bed"

sorted_bam="${bamcustom_dir}/sorted_${sample_name}_custom.bam"
sorted_tmp="${bamcustom_dir}/.tmp_sorted_${sample_name}_custom.bam"

raw_modkit_bedgz="${modkit}/custom_${sample_name}.raw.bed.gz"
final_modkit_bedgz="${modkit}/custom_${sample_name}.bed.gz"
stats_out="${modkit}/results_${sample_name}.tsv"

###############################################################################
# INITIAL CHECKS
###############################################################################

echo "[INFO] Job started"
echo "[INFO] Date: $(date)"
echo "[INFO] Host: $(hostname)"
echo "[INFO] PWD: $(pwd)"
echo "[INFO] CORES: $CORES"
echo "[INFO] VIEW_THREADS: $VIEW_THREADS"
echo "[INFO] FASTQ_THREADS: $FASTQ_THREADS"
echo "[INFO] MAP_THREADS: $MAP_THREADS"
echo "[INFO] SORT_THREADS: $SORT_THREADS"
echo "[INFO] INDEX_THREADS: $INDEX_THREADS"
echo "[INFO] SORT_MEM: $SORT_MEM"
echo "[INFO] sample: $sample"
echo "[INFO] reference: $reference"
echo "[INFO] bedsecuencias: $bedsecuencias"
echo "[INFO] outdir: $outdir"

test -s "$bedsecuencias"
test -s "$reference"
test -s "$sample"

samtools quickcheck -v "$sample"

if [[ ! -s "${reference}.fai" ]]; then
    echo "[INFO] Creating FASTA index: ${reference}.fai"
    samtools faidx "$reference"
fi

###############################################################################
# REMOVE PREVIOUS OUTPUTS
###############################################################################

echo "[INFO] Removing previous outputs"

rm -f "$samtools_bed"
rm -f "$modkit_bed"
rm -f "$regions_map"
rm -f "$original_regions_bed"

rm -f "$custom_fa" "${custom_fa}.fai"

rm -f "$sorted_bam" "${sorted_bam}.bai"
rm -f "$sorted_tmp" "${sorted_tmp}.bai"

rm -f "$raw_modkit_bedgz" "${raw_modkit_bedgz}.tbi"
rm -f "$final_modkit_bedgz" "${final_modkit_bedgz}.tbi"
rm -f "$stats_out"

rm -f "$allcontigs_dir"/*.fa
rm -f "$bamcustom_dir"/sort_tmp_*.bam
rm -f "$bamcustom_dir"/*.tmp.*.bam

###############################################################################
# ID / COORDINATE MAPPING TABLE
###############################################################################

printf "id\tchr\tinput_start\tfaidx_pos_1based\tcustom_chr\tcustom_insert_start_0based\tcustom_insert_end_0based\tleft_start_1based\tleft_end_1based\tright_start_1based\tright_end_1based\tsamtools_bed_start_0based\tsamtools_bed_end_0based\tinsert_len\n" \
    > "$regions_map"

###############################################################################
# 1) BUILD CUSTOM CONTIGS AND CUSTOM REFERENCE
###############################################################################

echo "[INFO] Building custom contigs"

while IFS=$'\t' read -r -a columnas; do

    [[ ${#columnas[@]} -lt 4 ]] && continue
    [[ -z "${columnas[0]}" ]] && continue
    [[ "${columnas[0]:0:1}" == "#" ]] && continue

    id="${columnas[0]}"
    chr="${columnas[1]}"
    input_start="${columnas[2]}"
    secuencia="${columnas[3]}"

    # Remove possible Windows carriage-return characters.
    id="${id//$'\r'/}"
    chr="${chr//$'\r'/}"
    input_start="${input_start//$'\r'/}"
    secuencia="${secuencia//$'\r'/}"

    if [[ ! "$input_start" =~ ^[0-9]+$ ]]; then
        echo "[ERROR] Non-numeric coordinate for id=${id}: ${input_start}" >&2
        exit 1
    fi

    Svlen=${#secuencia}

    if (( Svlen == 0 )); then
        echo "[WARNING] Empty insertion sequence for id=${id}; skipping." >&2
        continue
    fi

    # Convert the input coordinate to a 1-based faidx coordinate.
    if (( INPUT_START_IS_BED0 == 1 )); then
        pos_faidx=$(( input_start + 1 ))
    else
        pos_faidx=$input_start
    fi

    if (( pos_faidx < 1 )); then
        echo "[ERROR] Invalid faidx coordinate for id=${id}: ${pos_faidx}" >&2
        exit 1
    fi

    chrom_len=$(awk -v c="$chr" '$1 == c {print $2}' "${reference}.fai")

    if [[ -z "${chrom_len:-}" ]]; then
        echo "[ERROR] Chromosome ${chr} not found in ${reference}.fai for id=${id}" >&2
        exit 1
    fi

    if (( pos_faidx > chrom_len )); then
        echo "[ERROR] Coordinate ${pos_faidx} exceeds chromosome length ${chr}=${chrom_len} for id=${id}" >&2
        exit 1
    fi

    ###########################################################################
    # Define flanking regions
    #
    # Left flank: exactly `window` bases unless truncated by the chromosome start.
    # Right flank: exactly `window` bases unless truncated by the chromosome end.
    ###########################################################################

    Leftstart=$(( pos_faidx - window + 1 ))
    Leftend=$pos_faidx

    Rightstart=$(( pos_faidx + 1 ))
    Rightend=$(( pos_faidx + window ))

    if (( Leftstart < 1 )); then
        Leftstart=1
    fi

    if (( Rightend > chrom_len )); then
        Rightend=$chrom_len
    fi

    actual_left_len=$(( Leftend - Leftstart + 1 ))

    ###########################################################################
    # 0-based BED coordinates within the custom contig
    #
    # modkit_bed must contain EXACTLY 3 columns:
    # chr_custom    start    end
    ###########################################################################

    modkit_start=$actual_left_len
    modkit_end=$(( actual_left_len + Svlen ))

    custom_chr="chr_${id}"

    left_fa="${contigs_dir}/${id}_left.fa"
    right_fa="${contigs_dir}/${id}_right.fa"
    center_fa="${contigs_dir}/${id}_center.fa"

    left2_fa="${contigs_dir}/${id}_left2.fa"
    right2_fa="${contigs_dir}/${id}_right2.fa"
    merge_fa="${contigs_dir}/${id}_merge.fa"
    final_contig_fa="${allcontigs_dir}/contig_${id}.fa"

    ###########################################################################
    # Extract the left flank
    ###########################################################################

    samtools faidx "$reference" "${chr}:${Leftstart}-${Leftend}" > "$left_fa"

    ###########################################################################
    # Extract the right flank, if present
    ###########################################################################

    if (( Rightstart <= chrom_len )); then
        samtools faidx "$reference" "${chr}:${Rightstart}-${Rightend}" > "$right_fa"
    else
        : > "$right_fa"
        Rightstart=$(( chrom_len + 1 ))
        Rightend=$chrom_len
    fi

    ###########################################################################
    # Inserted/custom central sequence
    ###########################################################################

    printf "%s\n" "$secuencia" > "$center_fa"

    ###########################################################################
    # Prepare the custom-contig FASTA
    ###########################################################################

    seqkit seq -w 0 "$left_fa" \
        | sed "s/^>.*/>${custom_chr}/" \
        > "$left2_fa"

    if [[ -s "$right_fa" ]]; then
        seqkit seq -w 0 "$right_fa" \
            | awk '!/^>/' \
            > "$right2_fa"
    else
        : > "$right2_fa"
    fi

    cat "$left2_fa" "$center_fa" "$right2_fa" > "$merge_fa"

    seqkit seq -w 60 "$merge_fa" > "$final_contig_fa"

    ###########################################################################
    # BED for Modkit: EXACTLY 3 COLUMNS
    ###########################################################################

    printf "%s\t%d\t%d\n" \
        "$custom_chr" "$modkit_start" "$modkit_end" \
        >> "$modkit_bed"

    ###########################################################################
    # BED for `samtools view -L`: also kept as 3 columns
    ###########################################################################

    samtools_bed_start=$(( Leftstart - 1 ))
    samtools_bed_end=$Rightend

    printf "%s\t%d\t%d\n" \
        "$chr" "$samtools_bed_start" "$samtools_bed_end" \
        >> "$samtools_bed"

    ###########################################################################
    # Original-reference BED with 3 columns
    ###########################################################################

    printf "%s\t%d\t%d\n" \
        "$chr" "$samtools_bed_start" "$samtools_bed_end" \
        >> "$original_regions_bed"

    ###########################################################################
    # Auxiliary table containing insertion IDs and coordinate mappings.
    # This file is NOT used as the Modkit BED.
    ###########################################################################

    printf "%s\t%s\t%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" \
        "$id" \
        "$chr" \
        "$input_start" \
        "$pos_faidx" \
        "$custom_chr" \
        "$modkit_start" \
        "$modkit_end" \
        "$Leftstart" \
        "$Leftend" \
        "$Rightstart" \
        "$Rightend" \
        "$samtools_bed_start" \
        "$samtools_bed_end" \
        "$Svlen" \
        >> "$regions_map"

    rm -f \
        "$merge_fa" \
        "$left2_fa" \
        "$center_fa" \
        "$right2_fa" \
        "$left_fa" \
        "$right_fa"

done < "$bedsecuencias"

###############################################################################
# SORT BED FILES AND BUILD THE CUSTOM REFERENCE
###############################################################################

echo "[INFO] Sorting BED files"

test -s "$modkit_bed"
test -s "$samtools_bed"

LC_ALL=C sort -k1,1 -k2,2n "$modkit_bed" > "${modkit_bed}.tmp"
mv -f "${modkit_bed}.tmp" "$modkit_bed"

LC_ALL=C sort -k1,1 -k2,2n "$samtools_bed" > "${samtools_bed}.tmp"
mv -f "${samtools_bed}.tmp" "$samtools_bed"

LC_ALL=C sort -k1,1 -k2,2n "$original_regions_bed" > "${original_regions_bed}.tmp"
mv -f "${original_regions_bed}.tmp" "$original_regions_bed"

echo "[INFO] Building custom reference: $custom_fa"

cat "$allcontigs_dir"/*.fa > "$custom_fa"

test -s "$custom_fa"

samtools faidx "$custom_fa"

test -s "${custom_fa}.fai"

echo "[INFO] First lines of modkit.bed:"
head "$modkit_bed"

echo "[INFO] First lines of samtools.bed:"
head "$samtools_bed"

###############################################################################
# 2) EXTRACT READS, CONVERT TO FASTQ, REALIGN, AND SORT
###############################################################################

echo "[INFO] Checking disk space"

df -h "$bamcustom_dir" || true
df -ih "$bamcustom_dir" || true

if [[ -n "${SLURM_TMPDIR:-}" && -d "$SLURM_TMPDIR" ]]; then
    tmp_sort_parent="$SLURM_TMPDIR"
else
    tmp_sort_parent="$bamcustom_dir"
fi

tmp_sort_dir="${tmp_sort_parent}/samtools_sort_${SLURM_JOB_ID:-$$}"
mkdir -p "$tmp_sort_dir"

echo "[INFO] tmp_sort_dir: $tmp_sort_dir"
echo "[INFO] Temporary sorted BAM: $sorted_tmp"
echo "[INFO] Final sorted BAM: $sorted_bam"

rm -f "$sorted_tmp" "${sorted_tmp}.bai" "$sorted_bam" "${sorted_bam}.bai"

samtools view -@ "$VIEW_THREADS" -h -L "$samtools_bed" "$sample" \
    | samtools fastq -@ "$FASTQ_THREADS" -T '*' - \
    | minimap2 -k17 -t "$MAP_THREADS" -ax map-ont -y "$custom_fa" - \
    | samtools sort \
        -@ "$SORT_THREADS" \
        -m "$SORT_MEM" \
        -T "${tmp_sort_dir}/sort_tmp" \
        -o "$sorted_tmp" \
        -

samtools quickcheck -v "$sorted_tmp"

mv -f "$sorted_tmp" "$sorted_bam"

samtools quickcheck -v "$sorted_bam"

echo "[INFO] Indexing custom BAM"

samtools index -@ "$INDEX_THREADS" "$sorted_bam"

test -s "$sorted_bam"
test -s "${sorted_bam}.bai"

###############################################################################
# 3) MODKIT PILEUP
###############################################################################

echo "[INFO] Running modkit pileup"

modkit pileup \
    -t "$CORES" \
    --modified-bases 5mC 5hmC \
    --include-bed "$modkit_bed" \
    --reference "$custom_fa" \
    --bgzf \
    "$sorted_bam" \
    "$raw_modkit_bedgz"

test -s "$raw_modkit_bedgz"

###############################################################################
# 4) SORT AND INDEX BEDMETHYL OUTPUT
###############################################################################

echo "[INFO] Sorting Modkit output"

zcat "$raw_modkit_bedgz" \
    | LC_ALL=C sort -k1,1 -k2,2n \
    | bgzip -c > "$final_modkit_bedgz"

test -s "$final_modkit_bedgz"

tabix -f -p bed "$final_modkit_bedgz"

test -s "${final_modkit_bedgz}.tbi"

###############################################################################
# 5) MODKIT STATS
###############################################################################

echo "[INFO] Running modkit stats"

modkit stats \
    "$final_modkit_bedgz" \
    -t "$CORES" \
    --regions "$modkit_bed" \
    --mod-codes m,h \
    --min-coverage 10 \
    --out-table "$stats_out"

test -s "$stats_out"

###############################################################################
# FINAL SUMMARY
###############################################################################

echo "[INFO] Completed successfully"
echo "[INFO] Sorted custom BAM: $sorted_bam"
echo "[INFO] BAM index: ${sorted_bam}.bai"
echo "[INFO] Custom reference: $custom_fa"
echo "[INFO] 3-column Modkit BED: $modkit_bed"
echo "[INFO] 3-column samtools BED: $samtools_bed"
echo "[INFO] 3-column original-region BED: $original_regions_bed"
echo "[INFO] ID/coordinate mapping table: $regions_map"
echo "[INFO] Final BedMethyl: $final_modkit_bedgz"
echo "[INFO] Stats: $stats_out"
echo "[INFO] End date: $(date)"
