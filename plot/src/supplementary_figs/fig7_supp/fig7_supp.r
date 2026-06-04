#!/usr/bin/env Rscript


library(tidyverse)
library(cowplot)
library(data.table)


# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)


# Check if correct number of arguments provided
if (length(args) < 2) {
  cat("Usage: Rscript fig4_supp.r <output_prefix> <table_file> [max_bars_per_row]\n")
  cat("  <output_prefix>: Prefix for output files (without extension)\n")
  cat("  <table_file>: Path to merged table (TSV)\n")
  cat("  [max_bars_per_row]: Optional, maximum bars per row (default: 30)\n")
  quit(status = 1)
}


output_prefix <- args[1]
table_file <- args[2]
max_bars_per_row <- ifelse(length(args) >= 3, as.numeric(args[3]), 30)


cat("=== Single Dataset QV Comparison (Multi-row) ===\n")
cat(sprintf("Table file: %s\n", table_file))
cat(sprintf("Output prefix: %s\n", output_prefix))
cat(sprintf("Max bars per row: %d\n", max_bars_per_row))


# Color-blind friendly palette (Tol palette)
# For QV quality categories - earth tone gradient
qv_colors <- c(
  "very low: <= 17" = "#DDCC77",  # Sand/beige
  "low: >17, <= 23" = "#999933",  # Olive
  "mid: >23, <=33" = "#117733",   # Forest green
  "high: >33" = "#44AA99"         # Teal
)


# For error rate difference bars
# Orange (locityper) and Blue (cosigt) from Wong 2011 palette
performance_colors <- c(
  "locityper_better" = "#E69F00",  # Orange (locityper better)
  "cosigt_better" = "#0072B2"      # Blue (cosigt better)
)


# Read data
cat("\nReading data...\n")
if (!file.exists(table_file)) {
  stop(sprintf("Table file not found: %s", table_file))
}
merged_data <- fread(table_file)


cat(sprintf("Total rows: %d, Unique genes: %d\n", 
            nrow(merged_data), length(unique(merged_data$gene_name))))


# Function to categorize QV values (valid categories only)
categorize_qv <- function(qv_values) {
  factor(
    case_when(
      qv_values > 33         ~ "high: >33",
      qv_values > 23         ~ "mid: >23, <=33",
      qv_values > 17         ~ "low: >17, <= 23",
      qv_values <= 17        ~ "very low: <= 17",
      TRUE                   ~ NA_character_
    ),
    levels = c("very low: <= 17", "low: >17, <= 23", "mid: >23, <=33", "high: >33")
  )
}


# Process COSIGT QV data
cat("\n=== Processing COSIGT QV data ===\n")
qv_cols_cosigt <- c("QV_1_cosigt", "QV_2_cosigt")


qv_data_cosigt <- merged_data %>%
  select(sample, region, gene_name, all_of(qv_cols_cosigt)) %>%
  pivot_longer(
    cols = all_of(qv_cols_cosigt),
    names_to = "metric",
    values_to = "metric.values"
  ) %>%
  mutate(
    quality = categorize_qv(metric.values)
  ) %>%
  filter(!is.na(quality))  # Remove any NA quality values


# Summarize QV by gene
qv_summary <- qv_data_cosigt %>%
  count(gene_name, quality) %>%
  group_by(gene_name) %>%
  mutate(
    total = sum(n),
    percent = n / total * 100
  ) %>%
  ungroup()


# Order genes by percentage of high quality QVs
gene_order_df <- qv_summary %>%
  filter(quality == "high: >33") %>%
  arrange(desc(percent))


all_genes <- unique(qv_summary$gene_name)
missing_genes <- setdiff(all_genes, gene_order_df$gene_name)
if (length(missing_genes) > 0) {
  missing_df <- data.frame(gene_name = missing_genes, percent = 0)
  gene_order_df <- bind_rows(gene_order_df, missing_df)
}


gene_order <- gene_order_df$gene_name


qv_summary <- qv_summary %>%
  mutate(gene_name = factor(gene_name, levels = gene_order))


cat(sprintf("Processed %d genes\n", length(gene_order)))


# Calculate sample counts per gene
cat("\n=== Calculating sample counts ===\n")
gene_sample_counts <- merged_data %>%
  group_by(gene_name) %>%
  summarise(n_samples = n_distinct(sample), .groups = "drop") %>%
  mutate(gene_name = factor(gene_name, levels = gene_order))


cat(sprintf("Sample counts calculated for %d genes\n", nrow(gene_sample_counts)))


# Process error rate data
cat("\n=== Processing error rate data ===\n")


error_data <- merged_data %>%
  group_by(gene_name) %>%
  summarise(
    avg_error_rate_diff = mean(avg_error_rate_locityper - avg_error_rate_cosigt, na.rm = TRUE),
    n_calls = n(),
    .groups = "drop"
  ) %>%
  mutate(
    gene_name = factor(gene_name, levels = gene_order),
    abs_error_rate_diff = abs(avg_error_rate_diff),
    performance = ifelse(avg_error_rate_diff > 0, "cosigt_better", "locityper_better")
  )


# Calculate global max for background bars
global_max_diff <- max(error_data$abs_error_rate_diff, na.rm = TRUE)


error_data <- error_data %>%
  mutate(max_diff = global_max_diff)


# Calculate layout
num_genes <- length(gene_order)
num_rows <- ceiling(num_genes / max_bars_per_row)
bars_per_row <- ceiling(num_genes / num_rows)


cat(sprintf("\nPlotting %d genes in %d rows (%d bars per row)\n", 
            num_genes, num_rows, bars_per_row))


# Determine middle row for y-axis labels
middle_row_qv <- ceiling(num_rows / 2)
middle_row_error <- ceiling(num_rows / 2)


# Create plots row by row
qv_plots <- list()
error_plots <- list()


for (i in 1:num_rows) {
  start_idx <- (i-1) * bars_per_row + 1
  end_idx <- min(i * bars_per_row, num_genes)
  if (start_idx > num_genes) break
  
  genes_in_row <- gene_order[start_idx:end_idx]
  
  # Filter data for this row
  row_qv_data <- qv_summary %>%
    filter(gene_name %in% genes_in_row) %>%
    mutate(gene_name = factor(gene_name, levels = genes_in_row))
  
  row_sample_counts <- gene_sample_counts %>%
    filter(gene_name %in% genes_in_row) %>%
    mutate(gene_name = factor(gene_name, levels = genes_in_row))
  
  # QV plot for this row - only middle row gets y-axis label
  qv_plot <- ggplot(row_qv_data, aes(x = gene_name, y = percent, fill = quality)) +
    geom_bar(stat = "identity", position = "stack", width = 0.8) +
    geom_text(
      data = row_sample_counts,
      aes(x = gene_name, y = 101, label = n_samples),
      inherit.aes = FALSE,
      angle = 90,
      hjust = -0.1,
      vjust = 0.5,
      size = 6,
      family = "Helvetica"
    ) +
    scale_fill_manual(
      values = qv_colors,
      name = expression(QV[pred]),
      drop = FALSE
    ) +
    labs(
      x = if(i == num_rows) "Gene" else "",
      y = if(i == middle_row_qv) expression("% of"~QV[pred]) else ""
    ) +
    theme_bw(base_family = "Helvetica") +
    theme(
      axis.title = element_text(size = 24, family = "Helvetica"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 24, family = "Helvetica"),
      axis.text.y = element_text(size = 24, family = "Helvetica"),
      legend.position = "none",
      plot.margin = margin(t = 5, r = 20, b = 5, l = 5, unit = "pt")
    ) +
    guides(fill = guide_legend(nrow = 1)) +
    scale_y_continuous(
      limits = c(0, 115),
      breaks = seq(0, 100, by = 25),
      expand = expansion(mult = c(0, 0))
    )
  
  qv_plots[[i]] <- qv_plot
  
  # Error plot for this row - only middle row gets y-axis label, no sample counts
  row_error_data <- error_data %>%
    filter(gene_name %in% genes_in_row) %>%
    mutate(gene_name = factor(gene_name, levels = genes_in_row))
  
  error_plot <- ggplot(row_error_data, aes(x = gene_name)) +
    # Background bar at max height
    geom_bar(aes(y = max_diff), stat = "identity", 
             fill = "gray90", color = "gray60", linetype = "dashed", 
             linewidth = 0.1, width = 0.8) +
    # Colored bar at actual value
    geom_bar(aes(y = abs_error_rate_diff, fill = performance), 
             stat = "identity", linewidth = 0.1, width = 0.8) +
    scale_fill_manual(
      values = performance_colors,
      labels = c(
        "locityper_better" = "locityper better",
        "cosigt_better" = "cosigt better"
      ),
      name = "Performance",
      drop = FALSE
    ) +
    labs(
      x = if(i == num_rows) "Gene" else "",
      y = if(i == middle_row_error) "Average error rate difference\n|locityper - cosigt|" else ""
    ) +
    theme_bw(base_family = "Helvetica") +
    theme(
      axis.title = element_text(size = 24, family = "Helvetica"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 24, family = "Helvetica"),
      axis.text.y = element_text(size = 24, family = "Helvetica"),
      legend.position = "none",
      plot.margin = margin(t = 5, r = 5, b = 5, l = 20, unit = "pt")
    ) +
    guides(fill = guide_legend(nrow = 1)) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.05))
    )
  
  error_plots[[i]] <- error_plot
}


# Extract legends - only valid QV categories
cat("\n=== Extracting legends ===\n")

# Create dummy QV data with all valid categories (no failed/unknown)
dummy_qv_data <- data.frame(
  gene_name = rep("dummy", 4),
  quality = factor(
    c("very low: <= 17", "low: >17, <= 23", "mid: >23, <=33", "high: >33"),
    levels = c("very low: <= 17", "low: >17, <= 23", "mid: >23, <=33", "high: >33")
  ),
  percent = rep(1, 4)
)

qv_plot_for_legend <- ggplot(dummy_qv_data, aes(x = gene_name, y = percent, fill = quality)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    values = qv_colors,
    name = expression(QV[pred]),
    drop = FALSE
  ) +
  theme_bw(base_family = "Helvetica") +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 24, family = "Helvetica"),
    legend.title = element_text(size = 24, family = "Helvetica")
  ) +
  guides(fill = guide_legend(nrow = 1))

qv_legend <- get_legend(qv_plot_for_legend)


# Create dummy error data with both performance categories
dummy_error_data <- data.frame(
  gene_name = rep("dummy", 2),
  performance = factor(
    c("locityper_better", "cosigt_better"),
    levels = c("locityper_better", "cosigt_better")
  ),
  abs_error_rate_diff = rep(1, 2)
)

error_plot_for_legend <- ggplot(dummy_error_data, aes(x = gene_name, y = abs_error_rate_diff, fill = performance)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = performance_colors,
    labels = c(
      "locityper_better" = "locityper better",
      "cosigt_better" = "cosigt better"
    ),
    name = "Performance",
    drop = FALSE
  ) +
  theme_bw(base_family = "Helvetica") +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 24, family = "Helvetica"),
    legend.title = element_text(size = 24, family = "Helvetica")
  ) +
  guides(fill = guide_legend(nrow = 1))

error_legend <- get_legend(error_plot_for_legend)


# Combine legends horizontally
combined_legend <- plot_grid(
  qv_legend, 
  error_legend, 
  ncol = 2, 
  rel_widths = c(1.25, 1)
)


# Combine plots for each row
cat("\n=== Creating combined plots ===\n")
row_combined_plots <- list()


for (i in 1:num_rows) {
  row_combined <- plot_grid(
    qv_plots[[i]],
    error_plots[[i]],
    ncol = 2,
    align = 'h',
    axis = 'tb',
    rel_widths = c(1.25, 1)
  )
  row_combined_plots[[i]] <- row_combined
}


# Stack all rows vertically
all_plots_stacked <- plot_grid(
  plotlist = row_combined_plots,
  ncol = 1,
  align = 'v',
  axis = 'lr'
)


# Add legend at the very bottom
final_plot <- plot_grid(
  all_plots_stacked,
  combined_legend,
  ncol = 1,
  rel_heights = c(1, 0.08)
)


# Calculate plot dimensions
plot_width <- max(22, bars_per_row * 0.55 + 7)
plot_height <- max(8, num_rows * 5.5)


cat(sprintf("Plot dimensions: %.1f x %.1f inches\n", plot_width, plot_height))


# Save in multiple formats
output_png <- paste0(output_prefix, ".multi_row_comparison.png")
cat(sprintf("Saving PNG: %s\n", output_png))
ggsave(output_png, plot = final_plot, 
       width = plot_width, height = plot_height, 
       dpi = 600, limitsize = FALSE)


output_pdf <- paste0(output_prefix, ".multi_row_comparison.pdf")
cat(sprintf("Saving PDF: %s\n", output_pdf))
ggsave(output_pdf, plot = final_plot, 
       width = plot_width, height = plot_height, 
       limitsize = FALSE)


output_svg <- paste0(output_prefix, ".multi_row_comparison.svg")
cat(sprintf("Saving SVG: %s\n", output_svg))
ggsave(output_svg, plot = final_plot, 
       width = plot_width, height = plot_height, 
       device = "svg", limitsize = FALSE)


output_eps <- paste0(output_prefix, ".multi_row_comparison.eps")
cat(sprintf("Saving EPS: %s\n", output_eps))
ggsave(output_eps, plot = final_plot, 
       width = plot_width, height = plot_height, 
       device = "eps", limitsize = FALSE)


# === Calculate summary statistics ===
cat("\n=== Summary Statistics ===\n")


# QV statistics for BOTH cosigt and locityper
qv_cols_locityper <- c("QV_1_locityper", "QV_2_locityper")


qv_data_locityper <- merged_data %>%
  select(sample, region, gene_name, all_of(qv_cols_locityper)) %>%
  pivot_longer(
    cols = all_of(qv_cols_locityper),
    names_to = "metric",
    values_to = "metric.values"
  ) %>%
  mutate(
    quality = categorize_qv(metric.values),
    tool = "locityper"
  ) %>%
  filter(!is.na(quality))


qv_data_cosigt <- qv_data_cosigt %>%
  mutate(tool = "cosigt")


combined_qv <- bind_rows(qv_data_cosigt, qv_data_locityper)


qv_stats <- combined_qv %>%
  group_by(tool) %>%
  summarise(
    total_calls = n(),
    pct_high = 100 * sum(quality == "high: >33") / n(),
    pct_mid_or_higher = 100 * sum(quality %in% c("high: >33", "mid: >23, <=33")) / n(),
    pct_low_or_worse = 100 * sum(quality %in% c("low: >17, <= 23", "very low: <= 17")) / n(),
    .groups = "drop"
  )


# Save QV statistics
qv_stats_file <- paste0(output_prefix, ".qv_statistics.tsv")
write_tsv(qv_stats, qv_stats_file)
cat(sprintf("Saved: %s\n", qv_stats_file))


# Error rate statistics
error_stats <- error_data %>%
  summarise(
    n_genes = n(),
    cosigt_better = sum(performance == "cosigt_better"),
    locityper_better = sum(performance == "locityper_better"),
    pct_cosigt_better = 100 * sum(performance == "cosigt_better") / n(),
    pct_locityper_better = 100 * sum(performance == "locityper_better") / n(),
    median_abs_diff = median(abs_error_rate_diff, na.rm = TRUE),
    mean_abs_diff = mean(abs_error_rate_diff, na.rm = TRUE)
  )


# Save error statistics
error_stats_file <- paste0(output_prefix, ".error_statistics.tsv")
write_tsv(error_stats, error_stats_file)
cat(sprintf("Saved: %s\n", error_stats_file))


# Save QV summary
qv_summary_file <- paste0(output_prefix, ".qv_summary.tsv")
write_tsv(qv_summary, qv_summary_file)
cat(sprintf("Saved: %s\n", qv_summary_file))


# Save error data
error_data_file <- paste0(output_prefix, ".error_data.tsv")
write_tsv(error_data, error_data_file)
cat(sprintf("Saved: %s\n", error_data_file))


# Save sample counts
sample_counts_file <- paste0(output_prefix, ".sample_counts.tsv")
write_tsv(gene_sample_counts, sample_counts_file)
cat(sprintf("Saved: %s\n", sample_counts_file))


cat("\n=== Done! ===\n")

