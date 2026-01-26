## final comparison plot

library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(scales)
library(reshape2)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)


if (length(args) < 2) {
  cat("Usage: Rscript fig_e1.R <output_prefix> <hla_results>  \n")
  cat("  <output_prefix>: Prefix for output files\n")
  cat("  <hla_results>: Table of HLA typing results (output of get_hla_4dig.R) \n")
  quit(status = 1)
}


output_prefix <- args[1]
hla_table <- args[2]

cat("=== HLA typing results comparison ===\n")
cat(sprintf("Output prefix: %s\n", output_prefix))


hla_res <- read.table(hla_table,header=T)%>% select(gene,samples_with_data_poss,haplotype_accuracy_poss_cosigt,haplotype_accuracy_poss_t1k)
hla_res<-reshape2::melt(hla_res,measure.vars=c("haplotype_accuracy_poss_cosigt","haplotype_accuracy_poss_t1k"))
hla_res$Tool<-gsub("haplotype_accuracy_poss_","",hla_res$variable)

palette<-viridis(n = 2, option = "mako", direction = -1, begin = 0.2, end = 0.9)

samples_with_data<-hla_res %>% group_by(gene,samples_with_data_poss) %>% summarise(y=max(value))

hla_plot<-ggplot(hla_res)+geom_col(aes(x=gene,y=value,fill=Tool),position="dodge")+
  geom_text(data=samples_with_data,aes(x=gene,y=y+0.01,
                label=samples_with_data_poss),size=1.8,
            vjust = 0)+
  theme_minimal()+
  scale_x_discrete(name="Gene")+
  scale_fill_manual(values=palette)+
  scale_y_continuous(name="% of correct predictions",labels = scales::percent)+
  theme(text=element_text(size=5,family="Helvetica"),legend.key.size = unit(0.2, "cm"),
        axis.text = element_text(size=5,family="Helvetica"),legend.text = element_text(size=5,family="Helvetica"),
        legend.position="top")

plot_width <- 8.8
plot_height <- 5


cat(sprintf("Plot dimensions: %.1f x %.1f cm\n", plot_width, plot_height))


# Save in multiple formats
output_png <- paste0(output_prefix, ".hla_comparison.png")
cat(sprintf("\nSaving PNG: %s\n", output_png))
ggsave(output_png, plot = hla_plot, 
       width = plot_width, height = plot_height, 
       dpi = 600, limitsize = FALSE,units = "cm")


output_pdf <- paste0(output_prefix, ".hla_comparison.pdf")
cat(sprintf("Saving PDF: %s\n", output_pdf))
ggsave(output_pdf, plot = hla_plot, 
       width = plot_width, height = plot_height, 
       limitsize = FALSE,units = "cm")


output_svg <- paste0(output_prefix, ".hla_comparison.svg")
cat(sprintf("Saving SVG: %s\n", output_svg))
ggsave(output_svg, plot = hla_plot, 
       width = plot_width, height = plot_height, 
       device = "svg", limitsize = FALSE,units = "cm")


output_eps <- paste0(output_prefix, ".hla_comparison.eps")
cat(sprintf("Saving EPS: %s\n", output_eps))
ggsave(output_eps, plot = hla_plot, 
       width = plot_width, height = plot_height, 
       device = "eps", limitsize = FALSE,units = "cm")

