#!/usr/bin/env Rscript


library(tidyverse)
library(cowplot)
library(data.table)
library(viridis)


# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)


# Check if correct number of arguments provided
if (length(args) < 2) {
  cat("Usage: Rscript fig7_supp.r <output_prefix> <table_file> [max_bars_per_row]\n")
  cat("  <output_prefix>: Prefix for output files (without extension)\n")
  cat("  <table_file>: Path to merged table (TSV)\n")
  cat("  [max_bars_per_row]: Optional, maximum bars per row (default: 30)\n")
  quit(status = 1)
}


output_prefix <- args[1]
table_file <- args[2]
max_bars_per_row <- ifelse(length(args) >= 3, as.numeric(args[3]), 30)


cat("=== Single Dataset QV Fraction Comparison (Multi-row) ===\n")
cat(sprintf("Table file: %s\n", table_file))
cat(sprintf("Output prefix: %s\n", output_prefix))
cat(sprintf("Max bars per row: %d\n", max_bars_per_row))


# Quintile color palette (viridis mako)
quintile_palette <- rev(viridis(n = 5, option = "mako", direction = -1, begin = 0.2, end = 0.9))
names(quintile_palette) <- c("Q5: 0.8-1.0", "Q4: 0.6-0.8", "Q3: 0.4-0.6", "Q2: 0.2-0.4", "Q1: 0.0-0.2")


# Read data
cat("\nReading data...\n")
if (!file.exists(table_file)) {
  stop(sprintf("Table file not found: %s", table_file))
}
merged_data <- fread(table_file)


cat(sprintf("Total rows: %d, Unique genes: %d\n", 
            nrow(merged_data), length(unique(merged_data$gene_name))))


# Process QV fraction data
cat("\n=== Processing QV fraction data ===\n")


qv_data <- merged_data %>%
  mutate(
    qv_sum_pred = QV_1_pred + QV_2_pred,
    qv_sum_best = QV_1_best + QV_2_best,
    qv_fraction = round(qv_sum_pred / qv_sum_best, 2),
    quintile = cut(
      qv_fraction,
      breaks = seq(0, 1, 0.2),
      include.lowest = TRUE,
      labels = c("Q1: 0.0-0.2", "Q2: 0.2-0.4", "Q3: 0.4-0.6", "Q4: 0.6-0.8", "Q5: 0.8-1.0")
    )
  ) %>%
  filter(!is.na(quintile))  # Remove any NA quintiles


# Summarize by gene and quintile
qv_summary <- qv_data %>%
  count(gene_name, quintile) %>%
  group_by(gene_name) %>%
  mutate(
    total = sum(n),
    percent = n / total * 100
  ) %>%
  ungroup()


# Order genes by percentage of Q5 (highest quintile)
gene_order_df <- qv_summary %>%
  filter(quintile == "Q5: 0.8-1.0") %>%
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
gene_sample_counts <- qv_data %>%
  group_by(gene_name) %>%
  summarise(n_samples = n_distinct(sample), .groups = "drop") %>%
  mutate(gene_name = factor(gene_name, levels = gene_order))


cat(sprintf("Sample counts calculated for %d genes\n", nrow(gene_sample_counts)))


# Calculate layout
num_genes <- length(gene_order)
num_rows <- ceiling(num_genes / max_bars_per_row)
bars_per_row <- ceiling(num_genes / num_rows)


cat(sprintf("\nPlotting %d genes in %d rows (%d bars per row)\n", 
            num_genes, num_rows, bars_per_row))


# Determine middle row for y-axis label
middle_row <- ceiling(num_rows / 2)


# Create plots row by row
plots <- list()


for (i in 1:num_rows) {
  start_idx <- (i-1) * bars_per_row + 1
  end_idx <- min(i * bars_per_row, num_genes)
  if (start_idx > num_genes) break
  
  genes_in_row <- gene_order[start_idx:end_idx]
  
  # Filter data for this row
  row_data <- qv_summary %>%
    filter(gene_name %in% genes_in_row) %>%
    mutate(gene_name = factor(gene_name, levels = genes_in_row))
  
  row_sample_counts <- gene_sample_counts %>%
    filter(gene_name %in% genes_in_row) %>%
    mutate(gene_name = factor(gene_name, levels = genes_in_row))
  
  # Create plot for this row
  p <- ggplot(row_data, aes(x = gene_name, y = percent, fill = quintile)) +
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
      values = quintile_palette,
      name = expression(QV[frac]),
      drop = FALSE
    ) +
    labs(
      x = if(i == num_rows) "Gene" else "",
      y = if(i == middle_row) expression("% of"~QV[frac]) else ""
    ) +
    theme_bw(base_family = "Helvetica") +
    theme(
      axis.title = element_text(size = 24, family = "Helvetica"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 24, family = "Helvetica"),
      axis.text.y = element_text(size = 24, family = "Helvetica"),
      legend.position = "none",
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")
    ) +
    guides(fill = guide_legend(nrow = 1)) +
    scale_y_continuous(
      limits = c(0, 115),
      breaks = seq(0, 100, by = 25),
      expand = expansion(mult = c(0, 0))
    )
  
  plots[[i]] <- p
}


# Extract legend - create dummy plot with all quintiles
cat("\n=== Extracting legend ===\n")

dummy_data <- data.frame(
  gene_name = rep("dummy", 5),
  quintile = factor(
    c("Q1: 0.0-0.2", "Q2: 0.2-0.4", "Q3: 0.4-0.6", "Q4: 0.6-0.8", "Q5: 0.8-1.0"),
    levels = c("Q1: 0.0-0.2", "Q2: 0.2-0.4", "Q3: 0.4-0.6", "Q4: 0.6-0.8", "Q5: 0.8-1.0")
  ),
  percent = rep(1, 5)
)

plot_for_legend <- ggplot(dummy_data, aes(x = gene_name, y = percent, fill = quintile)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    values = quintile_palette,
    name = expression(QV[frac]),
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

legend <- get_legend(plot_for_legend)


# Stack all rows vertically
cat("\n=== Creating combined plot ===\n")
all_plots_stacked <- plot_grid(
  plotlist = plots,
  ncol = 1,
  align = 'v',
  axis = 'lr'
)


# Add legend at the very bottom
final_plot <- plot_grid(
  all_plots_stacked,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.08)
)


# Calculate plot dimensions
plot_width <- max(22, bars_per_row * 0.55 + 7)
plot_height <- max(8, num_rows * 5.5)


cat(sprintf("Plot dimensions: %.1f x %.1f inches\n", plot_width, plot_height))


# Save in multiple formats
output_png <- paste0(output_prefix, ".qv_fraction_genes.png")
cat(sprintf("Saving PNG: %s\n", output_png))
ggsave(output_png, plot = final_plot, 
       width = plot_width, height = plot_height, 
       dpi = 600, limitsize = FALSE)


output_pdf <- paste0(output_prefix, ".qv_fraction_genes.pdf")
cat(sprintf("Saving PDF: %s\n", output_pdf))
ggsave(output_pdf, plot = final_plot, 
       width = plot_width, height = plot_height, 
       limitsize = FALSE)


output_svg <- paste0(output_prefix, ".qv_fraction_genes.svg")
cat(sprintf("Saving SVG: %s\n", output_svg))
ggsave(output_svg, plot = final_plot, 
       width = plot_width, height = plot_height, 
       device = "svg", limitsize = FALSE)


output_eps <- paste0(output_prefix, ".qv_fraction_genes.eps")
cat(sprintf("Saving EPS: %s\n", output_eps))
ggsave(output_eps, plot = final_plot, 
       width = plot_width, height = plot_height, 
       device = "eps", limitsize = FALSE)


# === Calculate summary statistics ===
cat("\n=== Summary Statistics ===\n")


# Quintile statistics
quintile_stats <- qv_summary %>%
  group_by(quintile) %>%
  summarise(
    n_calls = sum(n),
    n_genes = n_distinct(gene_name),
    .groups = "drop"
  ) %>%
  mutate(
    total_calls = sum(n_calls),
    pct_of_total = 100 * n_calls / total_calls
  )


# Save quintile statistics
quintile_stats_file <- paste0(output_prefix, ".quintile_statistics.tsv")
write_tsv(quintile_stats, quintile_stats_file)
cat(sprintf("Saved: %s\n", quintile_stats_file))


# Per-gene statistics
gene_stats <- qv_summary %>%
  group_by(gene_name) %>%
  summarise(
    total_calls = sum(n),
    pct_Q5 = 100 * sum(n[quintile == "Q5: 0.8-1.0"]) / total_calls,
    pct_Q4_or_higher = 100 * sum(n[quintile %in% c("Q5: 0.8-1.0", "Q4: 0.6-0.8")]) / total_calls,
    pct_Q3_or_higher = 100 * sum(n[quintile %in% c("Q5: 0.8-1.0", "Q4: 0.6-0.8", "Q3: 0.4-0.6")]) / total_calls,
    .groups = "drop"
  ) %>%
  left_join(gene_sample_counts, by = "gene_name") %>%
  arrange(desc(pct_Q5))


# Save per-gene statistics
gene_stats_file <- paste0(output_prefix, ".gene_statistics.tsv")
write_tsv(gene_stats, gene_stats_file)
cat(sprintf("Saved: %s\n", gene_stats_file))


# Save QV summary
qv_summary_file <- paste0(output_prefix, ".qv_fraction_summary.tsv")
write_tsv(qv_summary, qv_summary_file)
cat(sprintf("Saved: %s\n", qv_summary_file))


# Save sample counts
sample_counts_file <- paste0(output_prefix, ".sample_counts.tsv")
write_tsv(gene_sample_counts, sample_counts_file)
cat(sprintf("Saved: %s\n", sample_counts_file))


cat("\n=== Done! ===\n")

