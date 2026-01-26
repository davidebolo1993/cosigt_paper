#!/usr/bin/env Rscript


library(tidyverse)
library(cowplot)
library(data.table)


# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)


# Check if correct number of arguments provided
if (length(args) < 3) {
  cat("Usage: Rscript panel_b.r <output_prefix> <metric_type> <table1:exp1_name> [table2:exp2_name] ...\n")
  cat("  <output_prefix>: Prefix for output files (without extension)\n")
  cat("  <metric_type>: 'qv' or 'error_rate'\n")
  cat("  <tableN:expN_name>: Table path and experiment name separated by colon\n")
  quit(status = 1)
}


output_prefix <- args[1]
metric_type <- tolower(args[2])
exp_args <- args[3:length(args)]


# Validate metric type
if (!metric_type %in% c("qv", "error_rate")) {
  stop("metric_type must be 'qv' or 'error_rate'")
}


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


cat("=== Categorical Quality Comparison (HORIZONTAL BARS) ===\n")
cat(sprintf("Output prefix: %s\n", output_prefix))
cat(sprintf("Metric type: %s\n", metric_type))
cat(sprintf("Number of experiments: %d\n", length(merged_files)))
cat("\nExperiments:\n")
for (i in seq_along(merged_files)) {
  cat(sprintf("  %s: %s\n", experiment_names[i], merged_files[i]))
}


# Color-blind friendly palette (Tol palette) - same as panel D
qv_colors <- c(
  "failed" = "#000000",           # Black
  "unknown" = "#999999",          # Gray
  "very low: <= 17" = "#DDCC77",  # Sand/beige
  "low: >17, <= 23" = "#999933",  # Olive
  "mid: >23, <=33" = "#117733",   # Forest green
  "high: >33" = "#44AA99"         # Teal
)


# Function to categorize QV values
categorize_qv <- function(qv_values) {
  factor(
    case_when(
      is.na(qv_values) | (is.infinite(qv_values) & qv_values < 0) ~ "failed",
      qv_values > 33         ~ "high: >33",
      qv_values > 23         ~ "mid: >23, <=33",
      qv_values > 17         ~ "low: >17, <= 23",
      qv_values <= 17        ~ "very low: <= 17",
      TRUE                   ~ "unknown"
    ),
    levels = c("failed", "unknown", "very low: <= 17", "low: >17, <= 23", "mid: >23, <=33", "high: >33")
  )
}

# Function to categorize error rates (maintain QV correspondence)
# Error rate = 10^(-QV/10)
categorize_error_rate <- function(error_rates) {
  factor(
    case_when(
      is.infinite(error_rates) | is.na(error_rates) ~ "failed",
      error_rates > 0.02         ~ "very low: <= 17",   # QV <= 17
      error_rates > 0.005        ~ "low: >17, <= 23",   # QV 17-23
      error_rates > 0.0005       ~ "mid: >23, <=33",    # QV 23-33
      error_rates <= 0.0005      ~ "high: >33",         # QV > 33
      TRUE                       ~ "unknown"
    ),
    levels = c("failed", "unknown", "very low: <= 17", "low: >17, <= 23", "mid: >23, <=33", "high: >33")
  )
}


# Process all experiments
cat("\nProcessing experiments...\n")
all_category_data <- list()


for (i in seq_along(merged_files)) {
  experiment_file <- merged_files[i]
  experiment_name <- experiment_names[i]


  cat(sprintf("  Reading: %s\n", experiment_name))


  if (!file.exists(experiment_file)) {
    warning(sprintf("File not found: %s. Skipping.", experiment_file))
    next
  }


  merged_data <- fread(experiment_file)


  if (metric_type == "qv") {
    # Process QV values
    qv_cols_cosigt <- c("QV_1_cosigt", "QV_2_cosigt")
    qv_cols_locityper <- c("QV_1_locityper", "QV_2_locityper")


    # Cosigt QVs
    cosigt_data <- merged_data %>%
      select(sample, region, gene_name, all_of(qv_cols_cosigt)) %>%
      pivot_longer(
        cols = all_of(qv_cols_cosigt),
        names_to = "metric",
        values_to = "metric_value"
      ) %>%
      mutate(
        tool = "cosigt",
        experiment = experiment_name,
        quality = categorize_qv(metric_value)
      )


    # Locityper QVs
    locityper_data <- merged_data %>%
      select(sample, region, gene_name, all_of(qv_cols_locityper)) %>%
      pivot_longer(
        cols = all_of(qv_cols_locityper),
        names_to = "metric",
        values_to = "metric_value"
      ) %>%
      mutate(
        tool = "locityper",
        experiment = experiment_name,
        quality = categorize_qv(metric_value)
      )


  } else {  # error_rate
    # Process error rate values
    cosigt_data <- merged_data %>%
      select(sample, region, gene_name, avg_error_rate_cosigt) %>%
      mutate(
        tool = "cosigt",
        experiment = experiment_name,
        quality = categorize_error_rate(avg_error_rate_cosigt)
      )


    locityper_data <- merged_data %>%
      select(sample, region, gene_name, avg_error_rate_locityper) %>%
      mutate(
        tool = "locityper",
        experiment = experiment_name,
        quality = categorize_error_rate(avg_error_rate_locityper)
      )
  }


  # Combine both tools
  exp_data <- bind_rows(cosigt_data, locityper_data)
  all_category_data[[i]] <- exp_data
}


# Combine all experiments
combined_data <- bind_rows(all_category_data) %>%
  mutate(
    tool = factor(tool, levels = c("locityper", "cosigt")),  # Reversed order for better visual
    experiment = factor(experiment, levels = experiment_names)
  )


cat(sprintf("\nTotal data points: %d\n", nrow(combined_data)))


# Calculate summary statistics
cat("\nCalculating summary statistics...\n")


category_summary <- combined_data %>%
  group_by(experiment, tool, quality) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(experiment, tool) %>%
  mutate(
    total = sum(count),
    percentage = (count / total) * 100
  ) %>%
  ungroup()


# Detailed statistics by experiment and tool
detailed_stats <- combined_data %>%
  group_by(experiment, tool) %>%
  summarise(
    total_calls = n(),
    n_failed = sum(quality == "failed"),
    pct_failed = 100 * sum(quality == "failed") / n(),
    n_unknown = sum(quality == "unknown"),
    pct_unknown = 100 * sum(quality == "unknown") / n(),
    pct_high = 100 * sum(quality == "high: >33") / n(),
    pct_mid_or_higher = 100 * sum(quality %in% c("high: >33", "mid: >23, <=33")) / n(),
    pct_low_or_worse = 100 * sum(quality %in% c("low: >17, <= 23", "very low: <= 17", "failed", "unknown")) / n(),
    .groups = "drop"
  ) %>%
  arrange(experiment, tool)


# Create the plot - HORIZONTAL BARS VERSION (tools on y-axis, percentage on x-axis)
cat("\nCreating categorical barplot (horizontal bars)...\n")


# Determine axis labels based on metric type
if (metric_type == "qv") {
  x_label <- expression("% of"~QV[pred])
  legend_title <- expression(QV[pred])
} else {
  x_label <- "% of Error Rate Calls"
  legend_title <- "Quality"
}


# HORIZONTAL BARS: y = tool, x = percentage, coord_flip NOT needed
categorical_plot <- ggplot(subset(category_summary, quality!="failed"), aes(y = tool, x = percentage, fill = quality)) +
  geom_col(stat = "identity", position = "stack") + 
  facet_grid(experiment ~ ., scales = "free_y", space = "free_y") +  # Vertical faceting by experiment
  scale_fill_manual(
    values = qv_colors,
    name = legend_title
  ) +
  labs(
    x = x_label,  # x-axis shows percentage (horizontal)
    y = "Tool"
  ) +
  theme_bw(base_family = "Helvetica") +
  theme(
    axis.title = element_text(size = 20, family = "Helvetica"),
    axis.text.x = element_text(size = 18, family = "Helvetica"),
    axis.text.y = element_text(size = 18, family = "Helvetica"),
    strip.text.y = element_text(size = 20, family = "Helvetica", angle = 0,
                                margin = margin(t = 5, b = 5, l = 5, r = 5, unit = "pt")),
    strip.background = element_blank(),  # No gray box
    panel.border = element_blank(),  # Remove panel borders
    axis.line = element_line(color = "black", linewidth = 0.5),  # Add axis lines instead
    panel.spacing.y = unit(0.3, "lines"),  # Minimal spacing between facets
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 18, family = "Helvetica"),
    legend.title = element_text(size = 18, family = "Helvetica"),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")
  ) +
  guides(fill = guide_legend(nrow = 1)) +
  scale_x_continuous(
    breaks = seq(0, 100, by = 25),
    expand = expansion(mult = c(0, 0.05))
  )


# Calculate plot dimensions - HORIZONTAL bars (wide plot, shorter height)
num_experiments <- length(experiment_names)
plot_width <- 10  # Fixed width for horizontal bars
plot_height <- max(6, num_experiments * 1.2)  # Height scales with experiments


cat(sprintf("Plot dimensions: %.1f x %.1f inches\n", plot_width, plot_height))


# Save in multiple formats
output_png <- paste0(output_prefix, ".horizontal.", metric_type, ".png")
cat(sprintf("Saving PNG: %s\n", output_png))
ggsave(output_png, plot = categorical_plot, 
       width = plot_width, height = plot_height, 
       dpi = 600, limitsize = FALSE)


output_pdf <- paste0(output_prefix, ".horizontal.", metric_type, ".pdf")
cat(sprintf("Saving PDF: %s\n", output_pdf))
ggsave(output_pdf, plot = categorical_plot, 
       width = plot_width, height = plot_height, 
       limitsize = FALSE)


output_svg <- paste0(output_prefix, ".horizontal.", metric_type, ".svg")
cat(sprintf("Saving SVG: %s\n", output_svg))
ggsave(output_svg, plot = categorical_plot, 
       width = plot_width, height = plot_height, 
       device = "svg", limitsize = FALSE)


output_eps <- paste0(output_prefix, ".horizontal.", metric_type, ".eps")
cat(sprintf("Saving EPS: %s\n", output_eps))
ggsave(output_eps, plot = categorical_plot, 
       width = plot_width, height = plot_height, 
       device = "eps", limitsize = FALSE)


# Save summary statistics (same as other versions)
cat("\nSaving summary statistics...\n")


category_summary_file <- paste0(output_prefix, ".category_summary_", metric_type, ".tsv")
write_tsv(category_summary, category_summary_file)
cat(sprintf("Saved: %s\n", category_summary_file))


detailed_stats_file <- paste0(output_prefix, ".detailed_stats_", metric_type, ".tsv")
write_tsv(detailed_stats, detailed_stats_file)
cat(sprintf("Saved: %s\n", detailed_stats_file))
