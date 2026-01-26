#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(tidyr)
library(rjson)
library(RColorBrewer)
library(scales)
library(data.table)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)


# Check if correct number of arguments provided
if (length(args) < 3) {
  cat("Usage: Rscript panel_c.r <output_folder> <input_folder> <region>\n")
  cat("  <output_folder>: Folder for output files\n")
  cat("  <input folder>: Folder with input files. Requires Cosigt intermediate outputs in the folder: 
      clustering data  (.clusters.medoids.tsv and .clusters.json),
      node length file from odgi/view (.node.length.tsv), 
      haplotype coverage file from odgi/paths (.tsv.gz)\n")
  cat("  <region>: Region coordinates in format chr_start_end\n")
  quit(status = 1)
}

output_folder <- args[1]
input_folder <- args[2]
region <- args[3]

setDTthreads(1)
viz_files<-c(paste0(input_folder,"/",region,".tsv.gz"),
        paste0(input_folder,"/",region,".node.length.tsv"),
        paste0(input_folder,"/",region,".clusters.json"),
        paste0(input_folder,"/",region,".clusters.medoids.tsv"))

coverage_file<-viz_files[1]
node_length_file<-viz_files[2]
cluster_file<-viz_files[3]
medoids<-fread(viz_files[4],header=F)

#mod nodes
pad_node_ids <- function(ids) {
  node_nums <- as.integer(sub("node\\.", "", ids))
  sprintf("node.%06d", node_nums)
}

#read inputs
#pangenome coverage over nodes 
coverage_data <- fread(coverage_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
coverage_long <- coverage_data %>%
  pivot_longer(cols = -path.name, names_to = "node_id", values_to = "coverage") %>%
  rename(path_name = path.name)
coverage_long$node_id <- pad_node_ids(coverage_long$node_id)

#node lengths
node_length_data <- fread(node_length_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
length_df<-data.frame(
  node_id = node_length_data$V1,
  length = node_length_data$V2,
  stringsAsFactors = FALSE
)
length_df$node_id<-pad_node_ids(length_df$node_id)

#clustering info
clustering_data<-fromJSON(file = cluster_file)
clustering_df <- data.frame(
  path_name = names(clustering_data),
  cluster = as.character(clustering_data),
  stringsAsFactors = FALSE
)

#merge
viz_data <- merge(coverage_long, length_df, by = "node_id")
viz_data <- merge(viz_data, clustering_df, by = "path_name", all.x = TRUE)

viz_data <- viz_data %>%
  arrange(node_id, path_name) %>% 
  group_by(path_name) %>%
  mutate(
    cumulative_length = cumsum(length),
    start_pos = lag(cumulative_length, default = 0),
    end_pos = cumulative_length
  ) %>%
  ungroup()


spectral_colors <- brewer.pal(11, "Spectral")  # 11-color spectral
max_cov_palette <- 2 + length(spectral_colors) - 1  # max coverage before clamping

color_map <- c(
  "0" = "white",
  "1" = "grey60"
)

for (i in 2:max_cov_palette) {
  color_map[as.character(i)] <- spectral_colors[i - 1]  # offset by 1
}

#mod viz data
viz_data$coverage_clamped <- as.character(
  ifelse(viz_data$coverage >= max_cov_palette, max_cov_palette, viz_data$coverage)
)


viz_data$cluster_num<-as.numeric(gsub("HaploGroup","",viz_data$cluster))
clusters<-unique(viz_data[,c("cluster_num","cluster")])
viz_data$cluster<-factor(viz_data$cluster,levels=clusters$cluster)

# select medoids, sort, get frequency info and change labels

viz_medoid<-viz_data[viz_data$path_name %in% medoids$V2,]
unique_mol<-viz_data[!duplicated(viz_data$path_name),]
abundance<-data.frame(table(unique_mol$cluster))
viz_medoid<-merge(viz_medoid,abundance,by.x="cluster",by.y="Var1")
viz_medoid$label<-paste0("HaploGroup",viz_medoid$cluster_num,"\n",viz_medoid$Freq," haplotype(s)")
viz_medoid<-viz_medoid[order(viz_medoid$Freq,decreasing = T),]
viz_medoid$path_name<-gsub("#hapl.*","",viz_medoid$path_name)
viz_medoid$label<-factor(viz_medoid$label,levels=unique(viz_medoid$label))

#plot

p <- ggplot(viz_medoid, 
            aes(x = start_pos, xend = end_pos, 
                y = path_name, 
                color = coverage_clamped)) +
  geom_segment(linewidth = 5) +
  scale_color_manual(
    values = color_map,
    na.value = "transparent"
  ) +
  scale_y_discrete(
    expand = expansion(add=c(0.2,0.2))
  ) +
  facet_grid(label ~ ., scales = "free_y", space = "free_y") +
  labs(
    x = "Genomic Position",
    y = "Haplotypes",
    color = "Coverage"
  ) +
  theme_bw() +
  theme(
    text=element_text(size=18,family="Helvetica"),
    strip.text.y = element_text(angle = 0, hjust = 0),
    strip.background = element_rect(fill= "white"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position="top")+ guides(colour = guide_legend(nrow = 1))

#save
plot_height=14
plot_width=33
ggsave(paste0(output_folder,"/figure1c.pdf"), height=plot_height, width=plot_width, limitsize = FALSE,units="cm",create.dir = T)
