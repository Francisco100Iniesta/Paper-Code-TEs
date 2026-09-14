#!/usr/bin/env Rscript

###############################################################################
# Pairwise methylation processing and descriptive statistics
#
# PURPOSE
# -------
# Process two methylation TSV files (sample A and sample B), retain loci that
# pass the selected quality threshold in BOTH samples, match the same loci
# across samples, calculate methylation-rate differences, and generate the
# descriptive statistics and threshold tables used for downstream reporting.
#
# This script replaces exploratory code that mixed filtering, statistics,
# annotation, and plotting in the same R session. Plotting is intentionally not
# included here; the output tables are designed to be reusable for figures.
#
# EXPECTED TSV COLUMNS
# --------------------
# Both sample files must contain:
#
#   #chrom
#   percent_m
#   count_valid_m
#
# For mode = "regional", they must also contain:
#
#   start
#   end
#
# These are the column names used by the TSV files processed in the original
# analysis.
#
# USAGE
# -----
# Regional methylation (run once for LEFT and once for RIGHT):
#
#   Rscript pairwise_methylation_statistics.R \
#       regional \
#       sample_A_left.tsv \
#       sample_B_left.tsv \
#       results/DT1_left \
#       A \
#       B \
#       100
#
# Internal TE methylation:
#
#   Rscript pairwise_methylation_statistics.R \
#       internal \
#       sample_A_internal.tsv \
#       sample_B_internal.tsv \
#       results/DT1_internal \
#       A \
#       B \
#       100
#
# ARGUMENTS
# ---------
# 1. mode
#      "regional" or "internal"
#
# 2. sample_A.tsv
# 3. sample_B.tsv
# 4. output_prefix
# 5. sample_A_label
# 6. sample_B_label
# 7. min_valid_calls
#
# QUALITY FILTER
# --------------
# A locus is retained only when `count_valid_m >= min_valid_calls` in BOTH
# samples.
#
# For regional 500-bp windows, the paper analysis used a minimum of 100 pooled
# valid modification calls per region.
#
# For internal TE methylation, this cleaned script leaves the threshold as an
# explicit argument so the exact final-analysis criterion is visible and
# reproducible.
#
# FINAL PAPER CUTOFF BASED ON QUALITY THRESHOLDS:
# If the final internal-TE analysis uses 100 valid calls, run this script with:
#
#     min_valid_calls = 100
#
# OUTPUTS
# -------
# <prefix>_paired_loci.tsv
#     One row per matched locus with sample A/B methylation, signed difference,
#     and absolute difference.
#
# <prefix>_summary.tsv
#     Number of compared loci, sample medians/IQRs, paired-difference summaries,
#     median absolute difference, and paired Hodges-Lehmann location shift with
#     95% confidence interval.
#
# <prefix>_thresholds.tsv
#     Counts and percentages of loci with |Delta MR| >= 5, 10, 15, and 20
#     percentage points, including direction of change.
#
# <prefix>_top50.tsv
#     The 50 loci with the largest absolute methylation-rate differences.
#
# STATISTICAL INTERPRETATION
# --------------------------
# Statistics are calculated at the matched-locus level. The paired
# Hodges-Lehmann estimator is obtained from a paired Wilcoxon signed-rank
# calculation with confidence intervals.
#
# As in the manuscript, these summaries should be interpreted descriptively:
# loci within the same individual are not independent biological replicates.
#
# REQUIREMENTS
# ------------
# R
# dplyr
# readr
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 7) {
  stop(
    paste0(
      "Usage: Rscript pairwise_methylation_statistics.R ",
      "<regional|internal> <sample_A.tsv> <sample_B.tsv> ",
      "<output_prefix> <sample_A_label> <sample_B_label> ",
      "<min_valid_calls>"
    )
  )
}

mode <- args[1]
sample_a_path <- args[2]
sample_b_path <- args[3]
output_prefix <- args[4]
sample_a_label <- args[5]
sample_b_label <- args[6]
min_valid_calls <- suppressWarnings(as.numeric(args[7]))

if (!mode %in% c("regional", "internal")) {
  stop("mode must be either 'regional' or 'internal'.")
}

if (!file.exists(sample_a_path)) {
  stop("Sample A TSV does not exist: ", sample_a_path)
}

if (!file.exists(sample_b_path)) {
  stop("Sample B TSV does not exist: ", sample_b_path)
}

if (is.na(min_valid_calls) || min_valid_calls < 0) {
  stop("min_valid_calls must be a non-negative number.")
}

read_methylation_tsv <- function(path, mode) {
  x <- read_tsv(
    path,
    show_col_types = FALSE,
    name_repair = "minimal"
  )

  required <- c("#chrom", "percent_m", "count_valid_m")

  if (mode == "regional") {
    required <- c(required, "start", "end")
  }

  missing <- setdiff(required, colnames(x))

  if (length(missing) > 0) {
    stop(
      "File ", path, " is missing required column(s): ",
      paste(missing, collapse = ", ")
    )
  }

  x
}

sample_a_raw <- read_methylation_tsv(sample_a_path, mode)
sample_b_raw <- read_methylation_tsv(sample_b_path, mode)

# Apply the same quality criterion independently to both samples before
# matching loci.
sample_a_filtered <- sample_a_raw %>%
  filter(count_valid_m >= min_valid_calls)

sample_b_filtered <- sample_b_raw %>%
  filter(count_valid_m >= min_valid_calls)

if (mode == "regional") {

  sample_a <- sample_a_filtered %>%
    transmute(
      chrom = `#chrom`,
      start = start,
      end = end,
      methylation_A = percent_m,
      valid_calls_A = count_valid_m
    )

  sample_b <- sample_b_filtered %>%
    transmute(
      chrom = `#chrom`,
      start = start,
      end = end,
      methylation_B = percent_m,
      valid_calls_B = count_valid_m
    )

  paired <- inner_join(
    sample_a,
    sample_b,
    by = c("chrom", "start", "end")
  )

} else {

  sample_a <- sample_a_filtered %>%
    transmute(
      TE_ID = `#chrom`,
      methylation_A = percent_m,
      valid_calls_A = count_valid_m
    )

  sample_b <- sample_b_filtered %>%
    transmute(
      TE_ID = `#chrom`,
      methylation_B = percent_m,
      valid_calls_B = count_valid_m
    )

  paired <- inner_join(
    sample_a,
    sample_b,
    by = "TE_ID"
  )
}

if (nrow(paired) == 0) {
  stop("No loci remained after filtering and pairwise matching.")
}

paired <- paired %>%
  mutate(
    delta_methylation = methylation_A - methylation_B,
    abs_delta_methylation = abs(delta_methylation),
    higher_methylation = case_when(
      delta_methylation > 0 ~ sample_a_label,
      delta_methylation < 0 ~ sample_b_label,
      TRUE ~ "equal"
    )
  ) %>%
  arrange(desc(abs_delta_methylation))

###############################################################################
# Paired Hodges-Lehmann location-shift estimator
###############################################################################

hl_test <- suppressWarnings(
  wilcox.test(
    paired$methylation_A,
    paired$methylation_B,
    paired = TRUE,
    conf.int = TRUE,
    conf.level = 0.95,
    exact = FALSE
  )
)

hl_estimate <- if (!is.null(hl_test$estimate)) {
  unname(hl_test$estimate)
} else {
  NA_real_
}

hl_conf_low <- if (!is.null(hl_test$conf.int)) {
  unname(hl_test$conf.int[1])
} else {
  NA_real_
}

hl_conf_high <- if (!is.null(hl_test$conf.int)) {
  unname(hl_test$conf.int[2])
} else {
  NA_real_
}

###############################################################################
# Summary table
###############################################################################

iqr_string <- function(x) {
  q <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE, names = FALSE)
  paste0(q[1], " - ", q[2])
}

summary_table <- tibble(
  mode = mode,
  sample_A = sample_a_label,
  sample_B = sample_b_label,
  min_valid_calls = min_valid_calls,
  loci_compared = nrow(paired),

  median_methylation_A = median(paired$methylation_A, na.rm = TRUE),
  IQR_methylation_A = iqr_string(paired$methylation_A),

  median_methylation_B = median(paired$methylation_B, na.rm = TRUE),
  IQR_methylation_B = iqr_string(paired$methylation_B),

  median_signed_difference = median(
    paired$delta_methylation,
    na.rm = TRUE
  ),
  IQR_signed_difference = iqr_string(
    paired$delta_methylation
  ),

  median_absolute_difference = median(
    paired$abs_delta_methylation,
    na.rm = TRUE
  ),

  hodges_lehmann_shift_A_minus_B = hl_estimate,
  hodges_lehmann_CI95_low = hl_conf_low,
  hodges_lehmann_CI95_high = hl_conf_high
)

###############################################################################
# Threshold table
###############################################################################

thresholds <- c(5, 10, 15, 20)

threshold_table <- bind_rows(
  lapply(thresholds, function(threshold) {

    selected <- paired %>%
      filter(abs_delta_methylation >= threshold)

    n_selected <- nrow(selected)

    n_a_higher <- sum(
      selected$delta_methylation > 0,
      na.rm = TRUE
    )

    n_b_higher <- sum(
      selected$delta_methylation < 0,
      na.rm = TRUE
    )

    n_equal <- sum(
      selected$delta_methylation == 0,
      na.rm = TRUE
    )

    tibble(
      threshold_abs_delta = threshold,
      loci = n_selected,
      percent_of_compared = 100 * n_selected / nrow(paired),
      sample_A_higher = n_a_higher,
      sample_A_higher_percent = if (n_selected > 0) {
        100 * n_a_higher / n_selected
      } else {
        NA_real_
      },
      sample_B_higher = n_b_higher,
      sample_B_higher_percent = if (n_selected > 0) {
        100 * n_b_higher / n_selected
      } else {
        NA_real_
      },
      equal = n_equal
    )
  })
)

###############################################################################
# Top-50 table
###############################################################################

top50 <- paired %>%
  slice_head(n = min(50, nrow(paired)))

###############################################################################
# Write outputs
###############################################################################

output_directory <- dirname(output_prefix)

if (!identical(output_directory, ".")) {
  dir.create(
    output_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

paired_path <- paste0(output_prefix, "_paired_loci.tsv")
summary_path <- paste0(output_prefix, "_summary.tsv")
threshold_path <- paste0(output_prefix, "_thresholds.tsv")
top50_path <- paste0(output_prefix, "_top50.tsv")

write_tsv(paired, paired_path)
write_tsv(summary_table, summary_path)
write_tsv(threshold_table, threshold_path)
write_tsv(top50, top50_path)

message("[INFO] Mode: ", mode)
message("[INFO] Quality cutoff: count_valid_m >= ", min_valid_calls)
message("[INFO] Loci compared: ", nrow(paired))
message("[INFO] Paired loci: ", paired_path)
message("[INFO] Summary: ", summary_path)
message("[INFO] Threshold table: ", threshold_path)
message("[INFO] Top 50: ", top50_path)
