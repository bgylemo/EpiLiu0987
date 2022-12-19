setwd("work_wd")
library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(rtracklayer)
library(GenomicRanges)
library(ggplotify)
library(karyoploteR)
library(UpSetR)
library(grid)
library(gridExtra)
library(VennDiagram)
library(raster)
library(rstatix)
library(ggpubr)

# read in ggplot theemes (made by Antonio Lentini) and plot parameters.
source("r_sources/AL_ggplot_themes.R")
source("r_sources/plot_parameters.R")

##### Read in data ##### 

# Annotation data
# tuki df
#tukiainen <- subset(fread("annotation_data/landscape.Suppl.Table.13.csv", header = T), !is.na(XCI_across_tissues) & XCI_across_tissues != "")
tukiainen <- fread("annotation_data/Suppl.Table.5.csv", header = T)
names(tukiainen)[names(tukiainen) == 'Gene name'] <- 'Gene_name'


# Supplementary table 4
suppl_table4 <- fread("supplementary_tables/Supplementary_Table4.tsv")

## load gene annotations
gtf_raw <- import("annotation_data/gencode.v41.annotation.gtf")
gtf_raw <- gtf_raw[seqnames(gtf_raw) == "chrX",]



#################### Supplementary Figure 1b ###########################

# Include only gene types observed in our data set
types_to_include <- unique(data.frame(gtf_raw[gtf_raw$type == "gene" & gtf_raw$gene_name %in% suppl_table4$gene ,])$gene_type)

markerz <- data.frame(gtf_raw[gtf_raw$type == "gene",])
markerz$hits <- "Gencode v.41"
markerz[markerz$gene_name %in% suppl_table4$gene & !markerz$gene_name %in% tukiainen$Gene_name,]$hits <- "Ours"
markerz[markerz$gene_name %in% tukiainen$Gene_name & markerz$gene_name %in% suppl_table4$gene,]$hits <- "Also_in_Tuki"

ggsave2(filename = paste0(plot_dir, "Supplementary_Figure_1b.pdf"), height = 4, width=3,
ggplot(data = markerz[markerz$gene_type %in% types_to_include,], aes(x=gene_type, y=..count.., fill = hits)) + 
  geom_bar(stat="count", position = "dodge") + 
  geom_text(stat="count", aes(label=..count..), position = position_dodge(width=0.85), size = 2)+ 
  theme_AL_box_rotX()+
  theme(axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6),
        axis.title = element_text(size=6))+
  labs(x="count", y="")
)
