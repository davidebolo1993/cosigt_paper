library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(viridis)

# plot qv vs region average difference

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)


# Check if correct number of arguments provided
if (length(args) < 3) {
  cat("Usage: Rscript fig15_supp.r <distance_file> <output_prefix> <table1:exp1_name> [table2:exp2_name] ...\n")
  cat("  <distance_file>: .tsv file with region distances with columns chrom, start, end, name, region, distance\n")
  cat("  <output_prefix>: Prefix for output files (without extension)\n")
  cat("  <tableN:expN_name>: Table path and experiment name separated by colon\n")
  quit(status = 1)
}

distance_file<-args[1]
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
cat(sprintf("Output prefix: %s\n", output_prefix))
cat(sprintf("Number of experiments: %d\n", length(merged_files)))
cat("\nExperiments:\n")
for (i in seq_along(merged_files)) {
  cat(sprintf("  %s: %s\n", experiment_names[i], merged_files[i]))
}

qvs_list<- list()

for(i in seq_along(merged_files)){
  cat(sprintf("\nProcessing %s...\n", experiment_names[i]))
  cat(sprintf("  Reading: %s\n", experiment_names[i]))
  
  if (!file.exists(merged_files[i])) {
    warning(sprintf("File not found: %s. Skipping.", merged_files[i]))
    next
  }
  
  qvs<-fread(merged_files[i])
  qvs$experiment<-experiment_names[i]
  qvs_list[[i]]<-qvs
  
}
qvs<-do.call(rbind,qvs_list)

qvs$experiment<-factor(qvs$experiment,levels=unique(qvs$experiment))

qvs_melt<-melt(qvs,id.vars=c("sample","region","gene_name","experiment"),measure.vars=c("QV_1_cosigt","QV_2_cosigt"))
qvs_by_region<-qvs_melt %>% group_by(experiment,region) %>% summarise(meanqv=mean(value))

regions<-read.table(distance_file,sep="\t",header=T)
qvs_by_region<-merge(qvs_by_region,regions,by.x=c("region","experiment"),by.y=c("region","condition"))



final_plot<-ggplot(qvs_by_region,aes(x=meandist,y=meanqv))+geom_point(alpha=0.5)+theme_minimal()+
  scale_x_continuous(name="Mean normalized cluster hapdist")+facet_wrap(facets=vars(experiment),nrow=2)+
  theme(strip.text = element_text(hjust=0,size=12),legend.position="inside",
        legend.justification = "right",legend.position.inside = c(0.98,0.2))+
  scale_y_continuous(name="Mean QVpred")

plot_width=16
plot_height=8

cat(sprintf("Plot dimensions: %.1f x %.1f inches\n", plot_width, plot_height))


# Save in multiple formats
output_png <- paste0(output_prefix, ".png")
cat(sprintf("Saving PNG: %s\n", output_png))
ggsave(output_png, plot = final_plot, 
       width = plot_width, height = plot_height, 
       dpi = 600, limitsize = FALSE, create.dir = T)


output_pdf <- paste0(output_prefix, ".pdf")
cat(sprintf("Saving PDF: %s\n", output_pdf))
ggsave(output_pdf, plot = final_plot, 
       width = plot_width, height = plot_height, 
       limitsize = FALSE, create.dir = T)


output_svg <- paste0(output_prefix, ".svg")
cat(sprintf("Saving SVG: %s\n", output_svg))
ggsave(output_svg, plot = final_plot, 
       width = plot_width, height = plot_height, 
       device = "svg", limitsize = FALSE, create.dir = T)



output_eps <- paste0(output_prefix, ".eps")
cat(sprintf("Saving EPS: %s\n", output_eps))
ggsave(output_eps, plot = final_plot, 
       width = plot_width, height = plot_height, 
       device = "eps", limitsize = FALSE, create.dir = T)


# qvs_lao<-fread("../../data/review/paper_benchmark_cosine_vs_qvs_vs_edr_CMRG_LAO.tsv")
# regions_lao<-read.table("/group/soranzo/davide.bolognini/working/stable/cosigt_paper/resources/benchmarking/cmrgs/cmrgs_refined.hprcy1.bed")
# regions_lao$region<-paste0(regions_lao$V1,"_",regions_lao$V2,"_",regions_lao$V3)
# regions_lao$meandist<-lapply(regions_lao$region,read_dissimilarity,"/project/ham/paper_results/cosigt_lao_hprcy1_on_hgsvcv3_cmr_FLAGGER/") %>% unlist()
# qv_lao_by_region<-qvs_lao %>% group_by(region) %>% summarise(meanqvsum=mean(qv_sum),meanqvmax=mean(qvmax),meanqvfrac=mean(qv_fraction))
# qv_lao_by_region<-merge(qv_lao_by_region,regions_lao,by.x="region",by.y="region")
# ggplot(qv_lao_by_region,aes(x=meandist,y=meanqvfrac))+geom_point(alpha=0.5)+ggtitle("CMRGs 30X LAO")+theme_minimal()+scale_x_continuous(name="Normalized cluster hapdist")#+scale_x_log10()#
# 
# ggsave("qvfrac_vs_meandist_CMRGs_LAO.pdf",width=16,height=8)
# 
# ggplot(qv_lao_by_region,aes(x=meandist,y=meanqvsum))+geom_point(alpha=0.5)+ggtitle("CMRGs 30X LAO")+theme_minimal()+scale_x_continuous(name="Normalized cluster hapdist")#+scale_x_log10()#+geom_smooth(method="lm",formula=y~x)
# ggplot(qv_lao_by_region,aes(x=meandist,y=meanqvmax))+geom_point(alpha=0.5)+ggtitle("CMRGs 30X LAO")+theme_minimal()+scale_x_continuous(name="Normalized cluster hapdist")#+scale_x_log10()#+geom_smooth(method="lm",formula=y~x)
# 
# 
# model<-lm(data=qv_lao_by_region,formula=meanqvfrac~meandist)
# 
# 
# categorize_qv <- function(qv_values) {
#   factor(
#     case_when(
#       is.na(qv_values) | (is.infinite(qv_values) & qv_values < 0) ~ "failed",
#       qv_values > 33*2         ~ "high: >66",
#       qv_values > 23*2         ~ "mid: >46, <=66",
#       qv_values > 17*2         ~ "low: >34, <= 46",
#       qv_values <= 17*2        ~ "very low: <= 34",
#       TRUE                   ~ "unknown"
#     ),
#     levels = c("failed", "unknown", "very low: <= 34", "low: >34, <= 46", "mid: >46, <=66", "high: >66")
#   )
# }
# 
# qv_colors <- c(
#   "failed" = "#000000",           # Black
#   "unknown" = "#999999",          # Gray
#   "very low: <= 34" = "#DDCC77",  # Sand/beige
#   "low: >34, <= 46" = "#999933",  # Olive
#   "mid: >46, <=66" = "#117733",   # Forest green
#   "high: >66" = "#44AA99"         # Teal
# )
# 
# 
# qvs<-merge(qvs,regions,by.x="region",by.y="region")
# qvs$qv_category<-categorize_qv(qvs$qvsum)
# ggplot(qvs,aes(fill=qv_category,x=qv_category,y=meandist))+geom_violin(alpha=0.5,lwd=0.25)+
#   geom_boxplot(width=0.1,show.legend = F,lwd=0.25,outlier.size = 0.25)+scale_fill_manual(values=qv_colors)+
#   scale_y_continuous(name="Normalized cluster hapdist")
