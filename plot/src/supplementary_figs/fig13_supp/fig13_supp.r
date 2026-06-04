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



# Create combined QV data (cosigt only)
cat("\n=== Creating QV comparison plot (cosigt only) ===\n")
qv_summary <- create_combined_qv_data(merged_files, experiment_names, gene_list_sorted)
qv_summary$experiment<-as.character(qv_summary$experiment)
qv_summary$cov_cont<-gsub(".*_","",qv_summary$experiment)
qv_summary<-qv_summary %>% separate(cov_cont,sep="-",into=c("cov","cont"),remove=F)
qv_summary<-qv_summary %>% arrange(cov,cont)
qv_summary$experiment<-factor(qv_summary$experiment,levels=unique(qv_summary$experiment))
# Create faceted QV plot - UPDATED with increased facet spacing
qv_plot <- ggplot(qv_summary, aes(x = gene_name, y = percent, fill = quality)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8,show.legend = F) +
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
    strip.text.y = element_text(size = 24,angle=0, family = "Helvetica", hjust=0),
    strip.background = element_rect(fill = "white"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 24, family = "Helvetica"),
    legend.title = element_text(size = 24, family = "Helvetica"),
    plot.margin = margin(t = 5, r = 20, b = 5, l = 5, unit = "pt"),
    panel.spacing.y = unit(1, "lines")
  ) +
  guides(fill = guide_legend(nrow = 1)) +
  scale_y_continuous(
    breaks = seq(0, 100, by = 25),
    expand = expansion(mult = c(0, 0.05))
  )

category_summary<-qv_summary %>% group_by(experiment,quality) %>% summarise(n=sum(n))
category_summary<-category_summary %>% group_by(experiment) %>% mutate(totalExp=sum(n))
category_summary$percentage<-category_summary$n/category_summary$totalExp*100

category_summary$quality<-factor(as.character(category_summary$quality),
                                 levels= c("very low: <= 17", "low: >17, <= 23", "mid: >23, <=33", "high: >33"))

# HORIZONTAL BARS: y = tool, x = percentage, coord_flip NOT needed
categorical_plot <- ggplot(category_summary, aes(y = experiment, x = percentage, fill = quality)) +
  geom_col(position = "stack") + 
  #  facet_grid(experiment ~ ., scales = "free_y", space = "free_y") +  # Vertical faceting by experiment
  scale_fill_manual(
    values = qv_colors,
    name = "QV category"
  ) +
  labs(
    x = "Percentage",  # x-axis shows percentage (horizontal)
    y = "Experiment"
  ) +
  theme_bw(base_family = "Helvetica") +
  theme(
    axis.title = element_text(size = 24, family = "Helvetica"),
    axis.text.x = element_text(size = 24, family = "Helvetica"),
    axis.text.y = element_blank(),
    strip.text.y = element_blank(),
    strip.background = element_blank(),  # No gray box
    panel.border = element_blank(),  # Remove panel borders
    axis.line = element_line(color = "black", linewidth = 0.5),  # Add axis lines instead
    panel.spacing.y = unit(1, "lines"),  # Minimal spacing between facets
    legend.position = "bottom",
    legend.text = element_text(size = 18, family = "Helvetica"),
    legend.title = element_text(size = 18, family = "Helvetica"),
    plot.margin = margin(t = 5, r = 5, b = 50, l = 5, unit = "pt")
  ) +
  guides(fill = guide_legend(nrow = 2,byrow=T)) +
  scale_x_continuous(
    breaks = seq(0, 100, by = 25),
    expand = expansion(mult = c(0, 0.05))
  )+facet_grid(experiment ~ ., scales = "free_y")

#categorical_plot

final_plot<-plot_grid(qv_plot,categorical_plot,axis="b",greedy=F,rel_widths = c(0.7,0.3))
#final_plot

# Calculate plot dimensions
num_genes <- length(gene_list_sorted)
num_experiments <- length(merged_files)
plot_width <- max(22, num_genes * 0.55 + 7)
plot_height <- min(20, num_experiments * 2.4)


cat(sprintf("Plot dimensions: %.1f x %.1f inches\n", plot_width, plot_height))


# Save in multiple formats
output_png <- paste0(output_prefix, ".multi_comparison.png")
cat(sprintf("Saving PNG: %s\n", output_png))
ggsave(output_png, plot = final_plot, 
       width = plot_width, height = plot_height, 
       dpi = 600, limitsize = FALSE, create.dir = T)


output_pdf <- paste0(output_prefix, ".multi_comparison.pdf")
cat(sprintf("Saving PDF: %s\n", output_pdf))
ggsave(output_pdf, plot = final_plot, 
       width = plot_width, height = plot_height, 
       limitsize = FALSE, create.dir = T)


output_svg <- paste0(output_prefix, ".multi_comparison.svg")
cat(sprintf("Saving SVG: %s\n", output_svg))
ggsave(output_svg, plot = final_plot, 
       width = plot_width, height = plot_height, 
       device = "svg", limitsize = FALSE, create.dir = T)



output_eps <- paste0(output_prefix, ".multi_comparison.eps")
cat(sprintf("Saving EPS: %s\n", output_eps))
ggsave(output_eps, plot = final_plot, 
       width = plot_width, height = plot_height, 
       device = "eps", limitsize = FALSE, create.dir = T)


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
    filter(region %in% gene_list_sorted)


  if (nrow(data_filtered) == 0) next


  # Process COSIGT
  
  # Combine and summarize
  combined_qv <- data_filtered

  qv_stats_exp <- combined_qv %>% mutate(totalcalls=sum(n)) %>%
    group_by(quality,totalcalls) %>%
    summarise(
      calls = sum(n))
  qv_stats_summary<-data.frame(quality=c("high","mid or higher","low or worse"),n=c(0,0,0),percent=c(0,0,0))
  qv_stats_summary[1,"n"]<-sum(qv_stats_exp[qv_stats_exp$quality=="high: >33","calls"])
  qv_stats_summary[2,"n"]<-sum(qv_stats_exp[qv_stats_exp$quality %in% c("high: >33", "mid: >23, <=33"),"calls"])
  qv_stats_summary[3,"n"]<-sum(qv_stats_exp[qv_stats_exp$quality %in% c("low: >17, <= 23", "very low: <= 17", "failed", "unknown"),"calls"])
  qv_stats_summary$percent <- (qv_stats_summary$n/unique(qv_stats_exp$totalcalls))*100
  qv_stats_summary$experiment<-experiment_name
  qv_stats_enhanced[[i]]<-qv_stats_summary
  
}
  combined_summary<-bind_rows(qv_stats_enhanced)
  write.table(combined_summary,file=paste0(output_prefix, "summary_table.tsv"),quote=F,row.names =F, col.names = T,sep="\t")
  
