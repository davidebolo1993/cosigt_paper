library(ggplot2)

files<-list.files(getwd(), pattern=".sorted.bed", recursive=T)
all<-list()


for (f in files) {

    df<-data.table::fread(f)
    sample<-unlist(strsplit(unlist(strsplit(f, "/"))[3], ".", fixed=T))[c(1,2)]
    sample_id<-paste(sample, collapse="#")
    all[[sample_id]]<-data.frame(value=df$V2, sample=sample_id)

}

all_df<-data.table::data.table(do.call(rbind,all))
p<-ggplot(all_df, aes(x=value))+
  geom_area(stat ="bin")+facet_wrap(~sample, ncol=10, nrow=10)

ggsave("freq.pdf", width=25, height=25)

range<-seq(min(all_df$value), max(all_df$value)-50000, by=50000)
found<-0

for (r in range) {

    if (found >= 5) break

    if (length(unique(all_df[value>=r & value<=r+50000]$sample)) == 100) {

        found<-found+1
        df<-data.frame(chrom="grch38#chr20", start=r, end=r+50000)
        data.table::fwrite(df, "roi.bed", append=TRUE, col.names=F, row.names=F, sep="\t", quote=F)
        freq<-table(all_df[value>=r & value<=r+50000]$sample)
        data.table::fwrite(data.frame(freq), paste0(r, ".tsv"), col.names=T, row.names=F, sep="\t", quote=F)

    }

}
