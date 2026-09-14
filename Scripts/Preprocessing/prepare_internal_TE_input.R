#!/usr/bin/env Rscript

###############################################################################
# Prepare input for internal TE methylation profiling
#
# PURPOSE
# -------
# Convert RetroInspector-derived insertion information into the four-column
# input required by `01_internal_TE_methylation.sh`.
#
# The original exploratory scripts used:
#   - a candidate table containing an insertion locus ("chr:position") and a
#     RetroInspector sequence identifier; and
#   - a RetroInspector insertion table containing `seqId` and `vcf_alt`.
#
# This cleaned version joins both sources and writes:
#
#   TE_ID    chromosome    insertion_position    insertion_sequence
#
# The output intentionally has NO header because `01_internal_TE_methylation.sh`
# reads it as a raw four-column tab-separated file.
#
# USAGE
# -----
#   Rscript prepare_internal_TE_input.R \
#       candidate_insertions.csv \
#       retroinspector_insertions.rds \
#       internal_TE_input.tsv
#
# INPUT 1: candidate table
# ------------------------
# Expected to contain at least two columns:
#   column 1: locus in the form chr:position
#   column 2: RetroInspector sequence ID
#
# A header is not required. Duplicate loci are reduced to the first occurrence,
# matching the behavior of the original analysis script.
#
# INPUT 2: RetroInspector insertion table
# ---------------------------------------
# Supported formats:
#   .rds              R object/data frame
#   .tsv / .txt       tab-delimited table with header
#   .csv              comma-delimited table with header
#
# Required columns:
#   seqId
#   vcf_alt
#
# OUTPUT
# ------
# Four-column, tab-separated file without a header:
#
#   TE_ID    chromosome    insertion_position    insertion_sequence
#
# REQUIREMENTS
# ------------
# R
# dplyr
# readr
# stringr
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    paste0(
      "Usage: Rscript prepare_internal_TE_input.R ",
      "<candidate_table> <retroinspector_table> <output.tsv>"
    )
  )
}

candidate_path <- args[1]
retroinspector_path <- args[2]
output_path <- args[3]

if (!file.exists(candidate_path)) {
  stop("Candidate table does not exist: ", candidate_path)
}

if (!file.exists(retroinspector_path)) {
  stop("RetroInspector table does not exist: ", retroinspector_path)
}

read_candidate_table <- function(path) {
  ext <- tolower(tools::file_ext(path))

  if (ext == "csv") {
    x <- read.csv(path, header = FALSE, stringsAsFactors = FALSE)
  } else {
    x <- read.delim(
      path,
      header = FALSE,
      sep = "\t",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  if (ncol(x) < 2) {
    stop("Candidate table must contain at least two columns.")
  }

  x <- x[, 1:2, drop = FALSE]
  colnames(x) <- c("locus", "seq_id")
  x
}

read_retroinspector_table <- function(path) {
  ext <- tolower(tools::file_ext(path))

  if (ext == "rds") {
    x <- readRDS(path)
  } else if (ext == "csv") {
    x <- read.csv(path, header = TRUE, stringsAsFactors = FALSE)
  } else {
    x <- read.delim(
      path,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  required <- c("seqId", "vcf_alt")
  missing <- setdiff(required, colnames(x))

  if (length(missing) > 0) {
    stop(
      "RetroInspector table is missing required column(s): ",
      paste(missing, collapse = ", ")
    )
  }

  x %>%
    select(seqId, vcf_alt) %>%
    distinct(seqId, .keep_all = TRUE)
}

candidate <- read_candidate_table(candidate_path)

# Keep one entry per locus, reproducing the de-duplication used in the
# exploratory analysis.
candidate <- candidate %>%
  group_by(locus) %>%
  slice(1) %>%
  ungroup()

locus_parts <- str_split_fixed(candidate$locus, ":", 2)

candidate <- candidate %>%
  mutate(
    chromosome = locus_parts[, 1],
    insertion_position = suppressWarnings(
      as.numeric(str_replace_all(locus_parts[, 2], "\\s+", ""))
    )
  )

if (anyNA(candidate$insertion_position)) {
  stop("At least one locus could not be parsed as 'chromosome:position'.")
}

retroinspector <- read_retroinspector_table(retroinspector_path)

internal_input <- candidate %>%
  left_join(retroinspector, by = c("seq_id" = "seqId")) %>%
  transmute(
    TE_ID = seq_id,
    chromosome = chromosome,
    insertion_position = insertion_position,
    insertion_sequence = vcf_alt
  )

missing_sequence <- is.na(internal_input$insertion_sequence) |
  internal_input$insertion_sequence == ""

if (any(missing_sequence)) {
  missing_ids <- internal_input$TE_ID[missing_sequence]
  stop(
    "No insertion sequence was found for TE ID(s): ",
    paste(missing_ids, collapse = ", ")
  )
}

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

write.table(
  internal_input,
  file = output_path,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

message("[INFO] Internal TE input written to: ", output_path)
message("[INFO] Insertions retained: ", nrow(internal_input))
