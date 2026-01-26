#!/usr/bin/env Rscript


library(cowplot)
library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(viridis)
setDTthreads(1)


# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)


if (length(args) < 2) {
  cat("Usage: Rscript panel_c.r <output_prefix> <merged_table1:condition1> [merged_table2:condition2] ...\n")
  cat("  <output_prefix>: Prefix for output files\n")
  cat("  <merged_tableN:conditionN>: Merged table path and condition name separated by colon\n")
  cat("                              Example: merged1.tsv:1X merged2.tsv:2X\n")
  quit(status = 1)
}


output_prefix <- args[1]
condition_args <- args[2:length(args)]


cat("=== QV Fraction Multi-Condition Comparison (STACKED) ===\n")
cat(sprintf("Output prefix: %s\n", output_prefix))
cat(sprintf("Number of conditions: %d\n", length(condition_args)))


# Parse condition information
condition_info <- lapply(condition_args, function(x) {
  parts <- strsplit(x, ":", fixed = TRUE)[[1]]
  if (length(parts) != 2) {
    stop(sprintf("Invalid format: %s. Expected format: merged_table:condition_name", x))
  }
  list(merged_file = parts[1], condition = parts[2])
})


cat("\nConditions:\n")
for (i in seq_along(condition_info)) {
  cat(sprintf("  %s: %s\n", 
              condition_info[[i]]$condition,
              condition_info[[i]]$merged_file))
}


# Function to calculate error rate from QV (for backward compatibility if needed)
qv_to_error_rate <- function(qv) {
  error_rate <- ifelse(is.infinite(qv) & qv < 0, 1, 10^(-qv/10))
  return(error_rate)
}


# Function to process merged QV data
process_qv_data <- function(merged_path, condition_name) {
  cat(sprintf("\nProcessing condition: %s\n", condition_name))

  # Read the merged table
  merged <- fread(merged_path)

  # Select and rename relevant columns
  # Map: QV_1_pred -> qv1, QV_2_pred -> qv2, QV_1_best -> best_qv1, QV_2_best -> best_qv2
  result <- merged[, .(
    sample = sample,
    region = region,
    gene = gene_name,
    qv1 = QV_1_pred,
    qv2 = QV_2_pred,
    best_qv1 = QV_1_best,
    best_qv2 = QV_2_best
  )]

  # Calculate QV fraction
  result[, qv_fraction := round((qv1 + qv2) / (best_qv1 + best_qv2), 2)]

  # Calculate error rates
  result[, error_rate1 := qv_to_error_rate(qv1)]
  result[, error_rate2 := qv_to_error_rate(qv2)]
  result[, avg_error_rate := (error_rate1 + error_rate2) / 2]
  result[, qv_sum := qv1 + qv2]
  result[, condition := condition_name]

  cat(sprintf("  Processed %d calls, %d unique samples, %d unique genes\n", 
              nrow(result),
              length(unique(result$sample)), 
              length(unique(result$gene))))

  return(result[, .(sample, region, gene, qv1, qv2, qv_sum, qv_fraction, 
                     error_rate1, error_rate2, avg_error_rate, condition)])
}


# Process all conditions
cat("\n=== Processing Data ===\n")
all_data <- list()


for (i in seq_along(condition_info)) {
  info <- condition_info[[i]]

  tryCatch({
    df <- process_qv_data(info$merged_file, info$condition)
    all_data[[i]] <- df
  }, error = function(e) {
    cat(sprintf("Error processing condition %s: %s\n", info$condition, e$message))
  })
}


if (length(all_data) == 0) {
  cat("\nError: No valid data was processed. Exiting.\n")
  quit(status = 1)
}


# Combine all data
qv_data <- rbindlist(all_data)
condition_names <- sapply(condition_info, function(x) x$condition)


cat(sprintf("\nTotal rows in merged data: %d\n", nrow(qv_data)))
cat(sprintf("Conditions: %s\n", paste(condition_names, collapse = ", ")))


# Save the full data table
error_rate_file <- paste0(output_prefix, ".error_rates.tsv")
cat(sprintf("\nSaving error rate table to: %s\n", error_rate_file))
fwrite(qv_data, file = error_rate_file, sep = "\t", quote = FALSE, 
       col.names = TRUE, row.names = FALSE)


# Calculate summary statistics per condition - using n_calls instead of n_samples
error_rate_summary <- qv_data %>%
  group_by(condition) %>%
  summarise(
    n_calls = n(),
    n_unique_samples = n_distinct(sample),
    n_genes = n_distinct(gene),
    median_error_rate = median(avg_error_rate, na.rm = TRUE),
    mean_error_rate = mean(avg_error_rate, na.rm = TRUE),
    min_error_rate = min(avg_error_rate, na.rm = TRUE),
    max_error_rate = max(avg_error_rate, na.rm = TRUE),
    median_qv_fraction = median(qv_fraction, na.rm = TRUE),
    mean_qv_fraction = mean(qv_fraction, na.rm = TRUE),
    .groups = "drop"
  )


error_rate_summary_file <- paste0(output_prefix, ".error_rate_summary.tsv")
cat(sprintf("Saving error rate summary to: %s\n", error_rate_summary_file))
fwrite(error_rate_summary, file = error_rate_summary_file, sep = "\t", 
       quote = FALSE, col.names = TRUE, row.names = FALSE)


# Create QV fraction plot - STACKED BARS like Panel B
cat("\n=== Creating Stacked QV Fraction Plot ===\n")


# Create quintile bins
qv_data <- qv_data %>%
  mutate(
    quintile = cut(
      qv_fraction,
      breaks = seq(0, 1, 0.2),
      include.lowest = TRUE,
      labels = c("Q1: 0.0-0.2", "Q2: 0.2-0.4", "Q3: 0.4-0.6", "Q4: 0.6-0.8", "Q5: 0.8-1.0")
    ),
    condition = factor(condition, levels = condition_names)
  )


# Use same colorblind-friendly palette (viridis mako) - matching the quintile order
quintile_palette <- rev(viridis(n = 5, option = "mako", direction = -1, begin = 0.2, end = 0.9))
names(quintile_palette) <- c("Q5: 0.8-1.0", "Q4: 0.6-0.8", "Q3: 0.4-0.6", "Q2: 0.2-0.4", "Q1: 0.0-0.2")


# Calculate percentages by condition and quintile - count = number of calls
qv_summary <- qv_data %>%
  count(condition, quintile) %>%
  group_by(condition) %>%
  mutate(
    total_calls = sum(n),
    percentage = n / total_calls * 100
  ) %>%
  ungroup() %>%
  filter(!is.na(quintile))  # Remove any NA quintiles


qv_summary_file <- paste0(output_prefix, ".qv_fraction_summary.tsv")
cat(sprintf("Saving QV fraction summary to: %s\n", qv_summary_file))
fwrite(qv_summary, file = qv_summary_file, sep = "\t", quote = FALSE, 
       col.names = TRUE, row.names = FALSE)


# Create STACKED HORIZONTAL BAR plot like Panel B
# One row per condition, stacked bars showing quintile distribution
cat("\nCreating stacked barplot...\n")


qv_fraction_plot <- ggplot(qv_summary, aes(y = condition, x = percentage, fill = quintile)) +
  geom_col(stat = "identity", position = "stack") + 
  scale_fill_manual(
    values = quintile_palette,
    name = expression(QV[frac])
  ) +
  labs(
    x = expression("% of"~QV[frac]),
    y = ""
  ) +
  theme_bw(base_family = "Helvetica") +
  theme(
    axis.title = element_text(size = 18, family = "Helvetica"),
    axis.text.x = element_text(size = 18, family = "Helvetica"),
    axis.text.y = element_text(size = 18, family = "Helvetica"),
    strip.text.y = element_text(size = 18, family = "Helvetica", angle = 0,
                                margin = margin(t = 5, b = 5, l = 5, r = 5, unit = "pt")),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    panel.spacing.y = unit(0.3, "lines"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 18, family = "Helvetica"),
    legend.title = element_text(size = 18, family = "Helvetica")
  ) +
  guides(fill = guide_legend(nrow = 1)) +
  scale_x_continuous(
    breaks = seq(0, 100, by = 25),
    expand = expansion(mult = c(0, 0.05))
  )


# Calculate dimensions
num_conditions <- length(condition_names)
plot_width <- 10
plot_height <- 4


cat(sprintf("Plot dimensions: %.1f x %.1f inches\n", plot_width, plot_height))


# Save in multiple formats
output_png <- paste0(output_prefix, ".qv_fraction_stacked.png")
cat(sprintf("\nSaving PNG: %s\n", output_png))
ggsave(output_png, plot = qv_fraction_plot, 
       width = plot_width, height = plot_height, 
       dpi = 600, limitsize = FALSE)


output_pdf <- paste0(output_prefix, ".qv_fraction_stacked.pdf")
cat(sprintf("Saving PDF: %s\n", output_pdf))
ggsave(output_pdf, plot = qv_fraction_plot, 
       width = plot_width, height = plot_height, 
       limitsize = FALSE)


output_svg <- paste0(output_prefix, ".qv_fraction_stacked.svg")
cat(sprintf("Saving SVG: %s\n", output_svg))
ggsave(output_svg, plot = qv_fraction_plot, 
       width = plot_width, height = plot_height, 
       device = "svg", limitsize = FALSE)


output_eps <- paste0(output_prefix, ".qv_fraction_stacked.eps")
cat(sprintf("Saving EPS: %s\n", output_eps))
ggsave(output_eps, plot = qv_fraction_plot, 
       width = plot_width, height = plot_height, 
       device = "eps", limitsize = FALSE)


cat("\n=== Done! ===\n")

