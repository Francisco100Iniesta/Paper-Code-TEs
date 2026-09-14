#!/usr/bin/env Rscript

###############################################################################
# Summarize sequencing quality-control metrics
#
# PURPOSE
# -------
# This script combines selected sequencing/alignment quality-control metrics
# into a compact summary table for downstream reporting.
#
# It was used as part of the sequencing quality-control workflow associated
# with the manuscript section:
#
#   Methods -> "Quality analysis of the sequencing samples"
#
# The manuscript describes quality assessment using NanoPlot, SAMtools stats,
# and Mosdepth. This script does not calculate those metrics itself; instead,
# it extracts selected values from precomputed QC output files and combines
# them into a single summary table.
#
# INPUTS
# ------
# The script expects three command-line arguments:
#
#   1. PATH1
#      A whitespace-delimited QC summary file containing:
#        - rows whose first column (V1) may end in "_region"
#        - a row containing "total" in V1
#        - the average coverage value in column V4 of that "total" row
#
#      In the original analysis this input was used to obtain the global
#      average coverage metric.
#
#   2. PATH2
#      A tab-delimited QC statistics file, without a header, containing one or
#      more rows with "average" in column V2.
#
#      These rows are retained and combined with the average-coverage metric.
#
#   3. OUT
#      Output directory.
#
# OUTPUT
# ------
#   <OUT>/stats_summary.tsv
#
# The output table contains:
#   - rows selected from PATH2 where V2 contains "average"
#   - one additional row named "Average_cover", obtained from PATH1
#
# USAGE
# -----
#   Rscript summarize_qc_metrics.R PATH1 PATH2 OUTPUT_DIR
#
# Example:
#
#   Rscript summarize_qc_metrics.R \
#       sample.mosdepth.summary.txt \
#       stats.txt \
#       qc_summary/
#
# REQUIREMENTS
# ------------
# R
# dplyr
#
# NOTES
# -----
# - This script summarizes already-computed QC metrics; it does not run
#   NanoPlot, SAMtools, or Mosdepth.
# - Column positions are inherited from the original analysis files.
# - If the format of the upstream QC files changes, the column selections
#   below may need to be adapted.
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
})

###############################################################################
# COMMAND-LINE ARGUMENTS
###############################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    paste0(
      "Usage: Rscript summarize_qc_metrics.R ",
      "<PATH1> <PATH2> <OUTPUT_DIR>"
    )
  )
}

PATH1 <- args[1]
PATH2 <- args[2]
OUT   <- args[3]

###############################################################################
# INPUT CHECKS
###############################################################################

if (!file.exists(PATH1)) {
  stop("PATH1 does not exist: ", PATH1)
}

if (!file.exists(PATH2)) {
  stop("PATH2 does not exist: ", PATH2)
}

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# 1) EXTRACT GLOBAL AVERAGE COVERAGE FROM PATH1
###############################################################################

coverage_table <- read.table(
  PATH1,
  header = FALSE,
  stringsAsFactors = FALSE
)

# Remove region-specific rows, retaining the global summary entries.
coverage_clean <- coverage_table %>%
  filter(!grepl("_region$", V1))

# Select the global "total" row used in the original analysis.
coverage_total <- coverage_clean %>%
  filter(grepl("total", V1))

if (nrow(coverage_total) == 0) {
  stop("No row containing 'total' was found in column V1 of PATH1.")
}

average_coverage <- coverage_total$V4

coverage_row <- data.frame(
  V2 = "Average_cover",
  V3 = average_coverage,
  stringsAsFactors = FALSE
)

###############################################################################
# 2) EXTRACT AVERAGE QC METRICS FROM PATH2
###############################################################################

stats <- read.delim(
  PATH2,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  comment.char = "#"
)

stats_average <- stats %>%
  filter(grepl("average", V2))

if (nrow(stats_average) == 0) {
  warning("No rows containing 'average' were found in column V2 of PATH2.")
}

# Retain the same columns used in the original analysis.
stats_average$V1 <- NULL
stats_average$V4 <- NULL

###############################################################################
# 3) COMBINE METRICS
###############################################################################

summary_table <- bind_rows(
  stats_average,
  coverage_row
)

###############################################################################
# 4) WRITE OUTPUT
###############################################################################

output_file <- file.path(
  OUT,
  "stats_summary.tsv"
)

write.table(
  summary_table,
  file = output_file,
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)

message("[INFO] QC summary written to: ", output_file)
