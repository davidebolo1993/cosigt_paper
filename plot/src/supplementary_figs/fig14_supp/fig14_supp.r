# plotting qv vs region length

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(viridis)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)


# Check if correct number of arguments provided
if (length(args) < 2) {
  cat("Usage: Rscript fig14_supp.r <output_prefix> <table1:exp1_name> [table2:exp2_name] ...\n")
  cat("  <output_prefix>: Prefix for output files (without extension)\n")
  cat("  <tableN:expN_name>: Table path and experiment name separated by colon\n")
  quit(status = 1)
}


output_prefix <- args[1]
exp_args <- args[2:length(args)]

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



## LEAVE-0-OUT

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
qvs_melt<-qvs_melt %>% separate(region,sep="_",into=c("chr","start","end"),remove=F,convert = T) %>% mutate(len=end-start)
qvs_by_region<-qvs_melt %>% group_by(experiment,region,len) %>% summarise(meanqv=mean(value))
qv_theoretical<-data.frame(len=seq(min(qvs_melt$len),max(qvs_melt$len),10))
qv_theoretical$qv<--10*log10(0.5/qv_theoretical$len)

colors_qv<-viridis(n = 3, option = "mako", direction = -1, begin = 0.2, end = 0.9)

names(colors_qv)<-c("mean qv per region","y ~ log(x)","theoretical QV max")

final_plot<-ggplot(qvs_by_region)+geom_point(aes(x=len/1000,y=meanqv,color="mean qv per region"),alpha=0.5)+
  geom_smooth(data=qvs_by_region,aes(x=len/1000,y=meanqv,color="y ~ log(x)"),method="lm",formula=y ~ log(x),fill="#3488A644")+
  geom_line(data=qv_theoretical,aes(x=len/1000,y=qv,color="theoretical QV max"))+
  facet_wrap(facets=vars(experiment),nrow=2,ncol=5,axes = "all")+
  theme_minimal()+#ggtitle("Mean QV by region vs region length - leave 0 out")+
  theme(strip.text = element_text(hjust=0,size=12),legend.position="inside",
        legend.justification = "right",legend.position.inside = c(0.95,0.25))+
  scale_x_continuous(name="Region length (kb)")+scale_color_manual(values=colors_qv,name=element_blank())+scale_y_continuous(name="Mean QVpred")
#ggsave("QV_by_SIZE.pdf",width=16,height=8)
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

# LEAVE-ALL-OUT
# 
# qvs_lao<-fread("../../fig2_main/panel_c_plot/panel_c.error_rates.tsv")
# 
# qvs_lao<-qvs_lao %>% separate(region,sep="_",into=c("chr","start","end"),remove=F,convert=T) %>% mutate(len=end-start)
# 
# qvs_lao_byregion<-qvs_lao %>% group_by(condition,region,gene,len) %>% summarise(mean_qv_frac=mean(qv_fraction))
# 
# ggplot(qvs_lao_byregion,aes(x=len/1000,y=mean_qv_frac))+geom_point(aes(color="mean qv per region"),alpha=0.5)+
#   facet_wrap(facets=vars(condition),nrow=2)+geom_smooth(aes(color="y ~ log(x)"),method="lm",formula=y~log(x),fill="#3488A644")+
#   theme_minimal()+
#   scale_x_continuous(name="length (kb)")+scale_color_manual(values=colors_qv)+scale_y_continuous(name="Mean QVfrac")+
#   ggtitle("Mean QV by region vs region length - leave all out")+
#   theme(strip.text = element_text(hjust=0,size=12))
# ggsave("qvs_vs_len_lao.pdf",width=6,height=6)

