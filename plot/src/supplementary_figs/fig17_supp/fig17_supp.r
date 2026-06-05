#!/usr/bin/env Rscript
#
# Usage:
#   Rscript fig17_supp.r <output_prefix> <benchmark_tsv> <annotation_bed>
#
#   <output_prefix>   prefix for all output files (no extension needed)
#   <benchmark_tsv>   cosigt/locityper merged TSV (may be gzipped)
#   <annotation_bed>  BED-style TSV with columns:
#                     #chrom  start  end  gene  classification  ...
#
# Outputs:
#   <prefix>.svtype_qv.png / .pdf / .svg
#   <prefix>.svtype_qv_summary.tsv

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# ── Arguments ──────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  cat("Usage: Rscript cosigt_svtype_qv.r <output_prefix> <benchmark_tsv> <annotation_bed>\n")
  quit(status = 1)
}

output_prefix   <- args[1]
benchmark_file  <- args[2]
annotation_file <- args[3]

cat("=== COSIGT QV by SV Type ===\n")
cat(sprintf("Benchmark   : %s\n", benchmark_file))
cat(sprintf("Annotation  : %s\n", annotation_file))
cat(sprintf("Output      : %s.*\n\n", output_prefix))

# ── Palettes (colour-blind friendly, Tol) ─────────────────────────────────────
qv_levels <- c("very low: <= 17", "low: >17, <= 23", "mid: >23, <=33", "high: >33")

qv_colors <- c(
  "very low: <= 17" = "#DDCC77",
  "low: >17, <= 23" = "#999933",
  "mid: >23, <=33"  = "#117733",
  "high: >33"       = "#44AA99"
)

# ── Helper: QV → quality category ─────────────────────────────────────────────
categorize_qv <- function(x) {
  factor(
    case_when(
      x >  33 ~ "high: >33",
      x >  23 ~ "mid: >23, <=33",
      x >  17 ~ "low: >17, <= 23",
      x <= 17 ~ "very low: <= 17",
      TRUE    ~ NA_character_
    ),
    levels = qv_levels
  )
}

# ── Read files ─────────────────────────────────────────────────────────────────
cat("Reading files ...\n")

benchmark  <- fread(benchmark_file)
annotation <- fread(annotation_file)

# BED files often have '#chrom'; normalise first column name
setnames(annotation, 1L, "chrom")

cat(sprintf("  benchmark  : %d rows, %d unique genes\n",
            nrow(benchmark), uniqueN(benchmark$gene_name)))
cat(sprintf("  annotation : %d genes, SV types: %s\n",
            nrow(annotation),
            paste(sort(unique(annotation$classification)), collapse = ", ")))

# ── Join SV-type information ───────────────────────────────────────────────────
sv_info <- annotation[, .(gene_name = gene, sv_type = classification)]
sv_info <- sv_info[!duplicated(gene_name)]   # guard against duplicate rows

data <- merge(benchmark, sv_info, by = "gene_name", all.x = TRUE)

n_missing <- sum(is.na(data$sv_type))
if (n_missing > 0)
  warning(sprintf(
    "%d benchmark rows have no annotation match (sv_type = NA) – dropped",
    n_missing))

data <- data[!is.na(sv_type)]

cat(sprintf("  after join : %d rows | %d genes | SV types: %s\n\n",
            nrow(data),
            uniqueN(data$gene_name),
            paste(sort(unique(data$sv_type)), collapse = ", ")))

# ── Tidy: one row per haplotype QV (COSIGT only) ──────────────────────────────
qv_long <- as_tibble(data) %>%
  select(sample, gene_name, sv_type, QV_1_cosigt, QV_2_cosigt) %>%
  pivot_longer(
    cols      = c(QV_1_cosigt, QV_2_cosigt),
    names_to  = "haplotype",
    values_to = "qv"
  ) %>%
  mutate(quality = categorize_qv(qv)) %>%
  filter(!is.na(quality))

cat(sprintf("QV calls: %d\n", nrow(qv_long)))

# ── SV-type level summary ──────────────────────────────────────────────────────
sv_summary <- qv_long %>%
  count(sv_type, quality, .drop = FALSE) %>%
  group_by(sv_type) %>%
  mutate(
    total   = sum(n),
    percent = if_else(total == 0, 0, n / total * 100)
  ) %>%
  ungroup()

# Fixed SV type order: DEL, INS, COMPLEX (bottom to top on a horizontal chart)
sv_levels <- c("DEL", "INS", "COMPLEX")

sv_summary <- sv_summary %>%
  mutate(
    sv_type = factor(sv_type, levels = sv_levels),
    quality = factor(quality, levels = qv_levels)
  )

# Annotation: n samples / n genes per SV type
sv_annot <- qv_long %>%
  group_by(sv_type) %>%
  summarise(
    n_samples = n_distinct(sample),
    n_regions = n_distinct(gene_name),
    .groups   = "drop"
  ) %>%
  mutate(
    sv_type = factor(sv_type, levels = sv_levels),
    label   = as.character(n_regions)
  )

# ── Plot ───────────────────────────────────────────────────────────────────────
cat("Building plot ...\n")

p <- ggplot(sv_summary, aes(y = sv_type, x = percent, fill = quality)) +
  geom_bar(stat = "identity", position = "stack",
           width = 0.55, colour = "white", linewidth = 0.3) +
  # region count to the right of each bar
  geom_text(
    data        = sv_annot,
    aes(y       = sv_type, x = 101, label = label),
    inherit.aes = FALSE,
    hjust       = 0, vjust = 0.5,
    size        = 5.5, family = "Helvetica"
  ) +
  scale_fill_manual(
    values = qv_colors,
    name   = expression(QV[pred]),
    drop   = FALSE
  ) +
  scale_x_continuous(
    limits = c(0, 110),
    breaks = seq(0, 100, 25),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = expression("% of" ~ QV[pred]),
    y = "SV type"
  ) +
  theme_bw(base_family = "Helvetica") +
  theme(
    axis.title         = element_text(size = 20),
    axis.text.x        = element_text(size = 16),
    axis.text.y        = element_text(size = 18),
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    legend.text        = element_text(size = 14),
    legend.title       = element_text(size = 14),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.margin        = margin(10, 10, 5, 10)
  ) +
  guides(fill = guide_legend(nrow = 1))

# ── Save plot ──────────────────────────────────────────────────────────────────
cat("Saving plots ...\n")
for (ext in c("png", "pdf", "svg")) {
  f <- sprintf("%s.svtype_qv.%s", output_prefix, ext)
  ggsave(f, plot = p, width = 11, height = 4.5,
         dpi = 600, limitsize = FALSE)
  cat(sprintf("  saved: %s\n", f))
}

# ── Summary statistics ─────────────────────────────────────────────────────────
cat("Computing summary ...\n")

qv_summary <- qv_long %>%
  group_by(sv_type) %>%
  summarise(
    n_calls           = n(),
    n_samples         = n_distinct(sample),
    n_genes           = n_distinct(gene_name),
    pct_high          = round(100 * mean(quality == "high: >33"),            2),
    pct_mid_or_higher = round(100 * mean(quality %in% c("high: >33",
                                                         "mid: >23, <=33")), 2),
    pct_low_or_worse  = round(100 * mean(quality %in% c("low: >17, <= 23",
                                                         "very low: <= 17")),2),
    median_qv         = round(median(qv, na.rm = TRUE), 3),
    mean_qv           = round(mean(qv,   na.rm = TRUE), 3),
    .groups = "drop"
  )

print(as.data.frame(qv_summary))

summary_out <- paste0(output_prefix, ".svtype_qv_summary.tsv")
write_tsv(qv_summary, summary_out)
cat(sprintf("  saved: %s\n", summary_out))

cat("\n=== Done! ===\n")
