#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(tidyr)
library(rjson)
library(RColorBrewer)
library(scales)
library(data.table)
library(gggenes)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)


# Check if correct number of arguments provided
if (length(args) < 3) {
  cat("Usage: Rscript panel_d.r <output_folder> <input_folder> <region>\n")
  cat("  <output_folder>: Folder for output files\n")
  cat("  <input folder>: Folder with input files. Requires Cosigt intermediate outputs in the folder: 
      clustering data  (.clusters.medoids.tsv and .clusters.json),
      fai index from bedtools/getfasta (.fasta.gz.fai),
      Pangene results from pangene/assemblies (.plot.bed.gz)")
  cat("  <region>: Region coordinates in format chr_start_end\n")
  quit(status = 1)
}

output_folder <- args[1]
input_folder <- args[2]
region <- args[3]


viz_files<-c(paste0(input_folder,"/",region,".plot.bed.gz"),
        paste0(input_folder,"/",region,".clusters.json"),
        paste0(input_folder,"/",region,".fasta.gz.fai"),
        paste0(input_folder,"/",region,".clusters.medoids.tsv"))

gene_data <- fread(viz_files[1])
cluster_mapping <- fromJSON(file=viz_files[2])    
fai_index<-fread(viz_files[3])
medoids<-fread(viz_files[4],header=F)

# add hapl start-end from fai
fai_index$start<-1
fai_index<-fai_index[,c(1,6,2)]
colnames(fai_index)<-c("molecule","mstart","mend")
gene_data<-merge(gene_data,fai_index,all.y=T)

# map clusters
clusters <- sapply(gene_data$molecule, function(mol) {
  cluster_mapping[[mol]]
})
gene_data$c <- clusters  

# order
gene_data <- gene_data[order(c, molecule)]
gene_data$start<-ifelse(is.na(gene_data$start),1,gene_data$start)
gene_data$end<-ifelse(is.na(gene_data$end),1,gene_data$end)
gene_data$strand<-ifelse(is.na(gene_data$strand),1,gene_data$strand)
#plot

mini_gene<-gene_data[gene_data$molecule %in% medoids$V2,]
unique_mol<-gene_data[!duplicated(gene_data$molecule),]
abundance<-data.frame(table(unique_mol$c))
mini_gene<-merge(mini_gene,abundance,by.x="c",by.y="Var1")
mini_gene$label<-paste0(mini_gene$c,"\n",mini_gene$Freq," haplotype(s)")
mini_gene<-mini_gene[order(mini_gene$Freq,decreasing = T),]
mini_gene$molecule<-gsub("#hapl.*","",mini_gene$molecule)
mini_gene$label<-factor(mini_gene$label,levels=unique(mini_gene$label))

p<-ggplot(mini_gene, aes(xmin=start, xmax=end, y=molecule,
                         fill=gene, forward=strand, label=gene )) +
  geom_gene_arrow(aes(xmin=mstart, xmax=mend, y=molecule,forward = T),fill="#eeeeee",
                  color="#777777",show.legend = F,
                  arrowhead_height = unit(0, "mm"),
                  arrowhead_width = unit(0, "mm"),arrow_body_height = unit(3,"mm"))+
  geom_gene_arrow(arrowhead_height = unit(3, "mm"),
                  arrowhead_width = unit(1, "mm"),
                  show.legend=TRUE,arrow_body_height = unit(3,"mm")) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() +
  labs(y="Haplotypes",x="Genomic Position") +
  facet_grid(rows=vars(label), scales="free_y", space="free_y") + 
  theme(text=element_text(size=18,family="Helvetica"),strip.text.y = element_text(angle = 0, hjust = 0),
        strip.background = element_rect(fill= "white"),legend.position="top",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) + 
  guides(fill = guide_legend(nrow = 1))+
  scale_y_discrete(
    expand = expansion(add=c(0.02,0.02))
  )

plot_height=14
plot_width=33
ggsave(paste0(output_folder,"/figure1d.pdf"), height=plot_height, width=plot_width, limitsize = FALSE,units="cm",create.dir = T)
