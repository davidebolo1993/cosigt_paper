#!/usr/bin/env Rscript


library(tidyverse)
library(cowplot)
library(data.table)


# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)


# Check if correct number of arguments provided
if (length(args) < 3) {
  cat("Usage: Rscript panel_d.r <gene_list_file> <output_prefix> <table1:exp1_name> [table2:exp2_name] ...\n")
  cat("  <gene_list_file>: Text file with one gene name per line\n")
  cat("  <output_prefix>: Prefix for output files (without extension)\n")
  cat("  <tableN:expN_name>: Table path and experiment name separated by colon\n")
  quit(status = 1)
}


gene_list_file <- args[1]
output_prefix <- args[2]
exp_args <- args[3:length(args)]


# Parse experiment names and file paths
experiment_info <- lapply(exp_args, function(x) {
  parts <- strsplit(x, ":", fixed = TRUE)[[1]]
  if (length(parts) != 2) {
    stop(sprintf("Invalid format: %s. Expected format: path:name", x))
  }
  list(file = parts[1], name = parts[2])
})


merged_files <- sapply(experiment_info, function(x) x$file)
experiment_names <- sapply(experiment_info, function(x) x$name)


cat("=== Multi-Experiment QV Comparison ===\n")
cat(sprintf("Gene list file: %s\n", gene_list_file))
cat(sprintf("Output prefix: %s\n", output_prefix))
cat(sprintf("Number of experiments: %d\n", length(merged_files)))
cat("\nExperiments:\n")
for (i in seq_along(merged_files)) {
  cat(sprintf("  %s: %s\n", experiment_names[i], merged_files[i]))
}


# Color-blind friendly palette (Tol palette)
# For QV quality categories - earth tone gradient
qv_colors <- c(
  "failed" = "#000000",           # Black
  "unknown" = "#999999",          # Gray
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


# Read gene list
cat("\nReading gene list...\n")
if (!file.exists(gene_list_file)) {
  stop(sprintf("Gene list file not found: %s", gene_list_file))
}
gene_list <- read_lines(gene_list_file) %>%
  str_trim() %>%
  .[. != ""]  # Remove empty lines


# Sort genes alphabetically for consistent plotting
gene_list_sorted <- sort(gene_list)
cat(sprintf("Genes to plot: %d genes\n", length(gene_list_sorted)))


# Function to process all experiments and create combined QV data for COSIGT ONLY
create_combined_qv_data <- function(merged_files, experiment_names, genes) {
  all_qv_data <- list()


  for (i in seq_along(merged_files)) {
    experiment_file <- merged_files[i]
    experiment_name <- experiment_names[i]


    cat(sprintf("\nProcessing %s...\n", experiment_name))
    cat(sprintf("  Reading: %s\n", experiment_file))


    if (!file.exists(experiment_file)) {
      warning(sprintf("File not found: %s. Skipping.", experiment_file))
      next
    }


    merged_data <- fread(experiment_file)
    cat(sprintf("  Total rows: %d, Unique genes: %d\n", 
                nrow(merged_data), length(unique(merged_data$gene_name))))


    # Filter to selected genes
    data_filtered <- merged_data %>%
      filter(gene_name %in% genes)


    if (nrow(data_filtered) == 0) {
      warning(sprintf("No data for selected genes in: %s", experiment_name))
      next
    }


    # Process COSIGT only
    qv_cols <- c("QV_1_cosigt", "QV_2_cosigt")


    qv_data <- data_filtered %>%
      select(sample, region, gene_name, all_of(qv_cols)) %>%
      pivot_longer(
        cols = all_of(qv_cols),
        names_to = "metric",
        values_to = "metric.values"
      ) %>%
      mutate(
        experiment = experiment_name
      )


    all_qv_data[[length(all_qv_data) + 1]] <- qv_data
  }


  # Combine all data
  combined <- bind_rows(all_qv_data)


  # Categorize QV values
  combined <- combined %>%
    mutate(
      quality = factor(
        case_when(
          is.infinite(metric.values) & metric.values < 0 ~ "failed",
          is.na(metric.values)       ~ "unknown",
          metric.values > 33         ~ "high: >33",
          metric.values > 23         ~ "mid: >23, <=33",
          metric.values > 17         ~ "low: >17, <= 23",
          metric.values <= 17        ~ "very low: <= 17",
          TRUE                       ~ "unknown"
        ),
        levels = c("failed", "unknown", "very low: <= 17", "low: >17, <= 23", "mid: >23, <=33", "high: >33")
      ),
      gene_name = factor(gene_name, levels = genes),
      experiment = factor(experiment, levels = experiment_names)
    )


  # Summarize
  qv_summary <- combined %>%
    count(experiment, gene_name, quality) %>%
    group_by(experiment, gene_name) %>%
    mutate(
      total = sum(n),
      percent = n / total * 100
    ) %>%
    ungroup()


  return(qv_summary)
}


# Function to create combined error rate data
create_combined_error_data <- function(merged_files, experiment_names, genes) {
  all_error_data <- list()


  for (i in seq_along(merged_files)) {
    experiment_file <- merged_files[i]
    experiment_name <- experiment_names[i]


    if (!file.exists(experiment_file)) next


    merged_data <- fread(experiment_file)


    data_filtered <- merged_data %>%
      filter(gene_name %in% genes)


    if (nrow(data_filtered) == 0) next


    diff_summary <- data_filtered %>%
      group_by(gene_name) %>%
      summarise(
        avg_error_rate_diff = mean(avg_error_rate_locityper - avg_error_rate_cosigt, na.rm = TRUE),
        n_samples = n(),
        .groups = "drop"
      ) %>%
      mutate(
        experiment = experiment_name,
        gene_name = factor(gene_name, levels = genes),
        abs_error_rate_diff = abs(avg_error_rate_diff),
        performance = ifelse(avg_error_rate_diff > 0, "cosigt_better", "locityper_better")
      )


    all_error_data[[i]] <- diff_summary
  }


  # Combine all data
  combined_error <- bind_rows(all_error_data) %>%
    mutate(experiment = factor(experiment, levels = experiment_names))


  return(combined_error)
}


# Create combined QV data (cosigt only)
cat("\n=== Creating QV comparison plot (cosigt only) ===\n")
qv_summary <- create_combined_qv_data(merged_files, experiment_names, gene_list_sorted)


# Create faceted QV plot - UPDATED with increased facet spacing
qv_plot <- ggplot(qv_summary, aes(x = gene_name, y = percent, fill = quality)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  facet_grid(experiment ~ ., scales = "fixed") +
  scale_fill_manual(
    values = qv_colors,
    name = expression(QV[pred])
  ) +
  labs(
    x = "Gene",
    y = expression("% of"~QV[pred])
  ) +
  theme_bw(base_family = "Helvetica") +
  theme(
    # Nature Methods: 5-7pt text, sans-serif (Helvetica) - but should be rescaled to the entire figure
    axis.title = element_text(size = 24, family = "Helvetica"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 24, family = "Helvetica"),
    axis.text.y = element_text(size = 24, family = "Helvetica"),
    strip.text.y = element_text(size = 24, family = "Helvetica", margin = margin(t = 8, b = 8, l = 5, r = 5, unit = "pt")),
    strip.background = element_rect(fill = "white"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 24, family = "Helvetica"),
    legend.title = element_text(size = 24, family = "Helvetica"),
    plot.margin = margin(t = 5, r = 20, b = 5, l = 5, unit = "pt"),
    panel.spacing.y = unit(1.2, "lines")
  ) +
  guides(fill = guide_legend(nrow = 1)) +
  scale_y_continuous(
    breaks = seq(0, 100, by = 25),
    expand = expansion(mult = c(0, 0.05))
  )


# Create combined error data
cat("\n=== Creating error rate plot ===\n")
error_data <- create_combined_error_data(merged_files, experiment_names, gene_list_sorted)


# Calculate GLOBAL max for background bars
global_max_diff <- max(error_data$abs_error_rate_diff, na.rm = TRUE)


error_data_with_max <- error_data %>%
  mutate(max_diff = global_max_diff)


# Create faceted error rate plot - UPDATED with increased facet spacing
error_plot <- ggplot(error_data_with_max, aes(x = gene_name)) +
  # Background bar at max height (light gray with dashed outline)
  geom_bar(aes(y = max_diff), stat = "identity", 
           fill = "gray90", color = "gray60", linetype = "dashed", linewidth = 0.1, width = 0.8) +
  # Colored bar at actual value
  geom_bar(aes(y = abs_error_rate_diff, fill = performance), stat = "identity", linewidth = 0.1, width = 0.8) +
  facet_grid(experiment ~ ., scales = "fixed") +
  scale_fill_manual(
    values = performance_colors,
    labels = c(
      "locityper_better" = "locityper better",
      "cosigt_better" = "cosigt better"
    ),
    name = "Performance"
  ) +
  labs(
    x = "Gene",
    y = "Average error rate difference\n|locityper - cosigt|"
  ) +
  theme_bw(base_family = "Helvetica") +
  theme(
    # Nature Methods: 5-7pt text, sans-serif (Helvetica)
    axis.title = element_text(size = 24, family = "Helvetica"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 24, family = "Helvetica"),
    axis.text.y = element_text(size = 24, family = "Helvetica"),
    strip.text.y = element_text(size = 24, family = "Helvetica", margin = margin(t = 8, b = 8, l = 5, r = 5, unit = "pt")),
    strip.background = element_rect(fill = "white"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 24, family = "Helvetica"),
    legend.title = element_text(size = 24, family = "Helvetica"),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 20, unit = "pt"),
    panel.spacing.y = unit(1.2, "lines") 
  ) +
  guides(fill = guide_legend(nrow = 1))


# Combine QV plot and error plot side by side
cat("\n=== Creating final combined plot ===\n")
final_plot <- plot_grid(
  qv_plot,
  error_plot,
  ncol = 2,
  align = 'h',
  axis = 'tb',
  rel_widths = c(1.25, 1),
  label_size = 22,
  label_fontface = "bold"
)


# Calculate plot dimensions
num_genes <- length(gene_list_sorted)
num_experiments <- length(merged_files)
plot_width <- max(22, num_genes * 0.55 + 7)
plot_height <- max(17, num_experiments * 2.2)


cat(sprintf("Plot dimensions: %.1f x %.1f inches\n", plot_width, plot_height))


# Save in multiple formats
output_png <- paste0(output_prefix, ".multi_comparison.png")
cat(sprintf("Saving PNG: %s\n", output_png))
ggsave(output_png, plot = final_plot, 
       width = plot_width, height = plot_height, 
       dpi = 600, limitsize = FALSE)


output_pdf <- paste0(output_prefix, ".multi_comparison.pdf")
cat(sprintf("Saving PDF: %s\n", output_pdf))
ggsave(output_pdf, plot = final_plot, 
       width = plot_width, height = plot_height, 
       limitsize = FALSE)


output_svg <- paste0(output_prefix, ".multi_comparison.svg")
cat(sprintf("Saving SVG: %s\n", output_svg))
ggsave(output_svg, plot = final_plot, 
       width = plot_width, height = plot_height, 
       device = "svg", limitsize = FALSE)



output_eps <- paste0(output_prefix, ".multi_comparison.eps")
cat(sprintf("Saving EPS: %s\n", output_eps))
ggsave(output_eps, plot = final_plot, 
       width = plot_width, height = plot_height, 
       device = "eps", limitsize = FALSE)


# === Calculate summary statistics for figure legend ===
cat("\n=== Summary Statistics ===\n")

# Enhanced QV statistics by experiment for BOTH cosigt and locityper
qv_stats_enhanced <- list()


for (i in seq_along(merged_files)) {
  experiment_file <- merged_files[i]
  experiment_name <- experiment_names[i]


  if (!file.exists(experiment_file)) next


  merged_data <- fread(experiment_file)


  data_filtered <- merged_data %>%
    filter(gene_name %in% gene_list_sorted)


  if (nrow(data_filtered) == 0) next


  # Process COSIGT
  qv_cols_cosigt <- c("QV_1_cosigt", "QV_2_cosigt")
  qv_data_cosigt <- data_filtered %>%
    select(sample, region, gene_name, all_of(qv_cols_cosigt)) %>%
    pivot_longer(
      cols = all_of(qv_cols_cosigt),
      names_to = "metric",
      values_to = "metric.values"
    ) %>%
    mutate(
      quality = case_when(
        is.infinite(metric.values) & metric.values < 0 ~ "failed",
        is.na(metric.values)       ~ "unknown",
        metric.values > 33         ~ "high: >33",
        metric.values > 23         ~ "mid: >23, <=33",
        metric.values > 17         ~ "low: >17, <= 23",
        metric.values <= 17        ~ "very low: <= 17",
        TRUE                       ~ "unknown"
      ),
      tool = "cosigt"
    )


  # Process LOCITYPER
  qv_cols_locityper <- c("QV_1_locityper", "QV_2_locityper")
  qv_data_locityper <- data_filtered %>%
    select(sample, region, gene_name, all_of(qv_cols_locityper)) %>%
    pivot_longer(
      cols = all_of(qv_cols_locityper),
      names_to = "metric",
      values_to = "metric.values"
    ) %>%
    mutate(
      quality = case_when(
        is.infinite(metric.values) & metric.values < 0 ~ "failed",
        is.na(metric.values)       ~ "unknown",
        metric.values > 33         ~ "high: >33",
        metric.values > 23         ~ "mid: >23, <=33",
        metric.values > 17         ~ "low: >17, <= 23",
        metric.values <= 17        ~ "very low: <= 17",
        TRUE                       ~ "unknown"
      ),
      tool = "locityper"
    )


  # Combine and summarize
  combined_qv <- bind_rows(qv_data_cosigt, qv_data_locityper)


  qv_stats_exp <- combined_qv %>%
    group_by(tool) %>%
    summarise(
      total_calls = n(),
      pct_high = 100 * sum(quality == "high: >33") / n(),
      pct_mid_or_higher = 100 * sum(quality %in% c("high: >33", "mid: >23, <=33")) / n(),
      pct_low_or_worse = 100 * sum(quality %in% c("low: >17, <= 23", "very low: <= 17", "failed", "unknown")) / n(),
      n_failed = sum(quality == "failed"),
      pct_failed = 100 * sum(quality == "failed") / n(),
      n_unknown = sum(quality == "unknown"),
      pct_unknown = 100 * sum(quality == "unknown") / n(),
      .groups = "drop"
    ) %>%
    mutate(experiment = experiment_name)


  qv_stats_enhanced[[i]] <- qv_stats_exp
}


qv_stats <- bind_rows(qv_stats_enhanced) %>%
  select(experiment, tool, total_calls, 
         n_failed, pct_failed, n_unknown, pct_unknown,
         pct_high, pct_mid_or_higher, pct_low_or_worse) %>%
  arrange(experiment, tool)

# Save QV statistics
qv_stats_file <- paste0(output_prefix, ".qv_statistics.tsv")
write_tsv(qv_stats, qv_stats_file)
cat(sprintf("Saved: %s\n", qv_stats_file))


# Error rate statistics by experiment
error_stats <- error_data %>%
  group_by(experiment) %>%
  summarise(
    n_genes = n(),
    cosigt_better = sum(performance == "cosigt_better"),
    locityper_better = sum(performance == "locityper_better"),
    pct_cosigt_better = 100 * sum(performance == "cosigt_better") / n(),
    pct_locityper_better = 100 * sum(performance == "locityper_better") / n(),
    median_abs_diff = median(abs_error_rate_diff, na.rm = TRUE),
    mean_abs_diff = mean(abs_error_rate_diff, na.rm = TRUE),
    .groups = "drop"
  )


# Save error statistics
error_stats_file <- paste0(output_prefix, ".error_statistics.tsv")
write_tsv(error_stats, error_stats_file)
cat(sprintf("Saved: %s\n", error_stats_file))


# Coverage comparison - by individual experiment
error_data_grouped <- error_data %>%
  group_by(experiment) %>%
  summarise(
    n_observations = n(),
    cosigt_better_count = sum(performance == "cosigt_better"),
    locityper_better_count = sum(performance == "locityper_better"),
    pct_cosigt_better = 100 * sum(performance == "cosigt_better") / n(),
    pct_locityper_better = 100 * sum(performance == "locityper_better") / n(),
    mean_abs_diff = mean(abs_error_rate_diff, na.rm = TRUE),
    median_abs_diff = median(abs_error_rate_diff, na.rm = TRUE),
    .groups = "drop"
  )


# Save coverage group comparison
coverage_stats_file <- paste0(output_prefix, ".coverage_comparison.tsv")
write_tsv(error_data_grouped, coverage_stats_file)
cat(sprintf("Saved: %s\n", coverage_stats_file))

cat("\n")
