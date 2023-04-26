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
source("project/Genome_biology_submission/r_sources/plot_parameters.R")

##### Read in data ##### 

# Supplementary table 4
suppl_table4 <- fread("project/Genome_biology_submission/supplementary_tables/Supplementary_Table4.tsv")

suppl_table7 <- fread("project/Genome_biology_submission/supplementary_tables/Supplementary_Table5.tsv")

# keep only highest covered hSNP.
suppl_table4_filtered <- data.table(suppl_table4[suppl_table4$keep == T,])

####### Suppl figure 1c #######
# get overlapping genes in overlapping tissues
# find overlapping genes
overlap1 <- intersect(unique(suppl_table4_filtered[suppl_table4_filtered$individual == "UPIC",]$gene), unique(suppl_table4_filtered[suppl_table4_filtered$individual == "ZZPU",]$gene))
overlap2 <- intersect(unique(suppl_table4_filtered[suppl_table4_filtered$individual == "UPIC",]$gene), unique(suppl_table4_filtered[suppl_table4_filtered$individual == "13PLJ",]$gene))
overlap3 <- intersect(unique(suppl_table4_filtered[suppl_table4_filtered$individual == "ZZPU",]$gene), unique(suppl_table4_filtered[suppl_table4_filtered$individual == "13PLJ",]$gene))

genes_overlap <- unique(c(overlap1, overlap2, overlap3))

# find overlapping tissues
overlap4 <- intersect(unique(suppl_table4_filtered[suppl_table4_filtered$individual == "UPIC",]$tissue_id), unique(suppl_table4_filtered[suppl_table4_filtered$individual == "ZZPU",]$tissue_id))
overlap5 <- intersect(unique(suppl_table4_filtered[suppl_table4_filtered$individual == "UPIC",]$tissue_id), unique(suppl_table4_filtered[suppl_table4_filtered$individual == "13PLJ",]$tissue_id))
overlap6 <- intersect(unique(suppl_table4_filtered[suppl_table4_filtered$individual == "ZZPU",]$tissue_id), unique(suppl_table4_filtered[suppl_table4_filtered$individual == "13PLJ",]$tissue_id))

UPIC_ZZPU <- na.omit(dcast(suppl_table4_filtered[suppl_table4_filtered$gene %in% overlap1 & suppl_table4_filtered$tissue_id %in% overlap4,], gene + tissue_id ~ individual, value.var = "allelic_expression"))
UPIC_13PLJ <- na.omit(dcast(suppl_table4_filtered[suppl_table4_filtered$gene %in% overlap2 & suppl_table4_filtered$tissue_id %in% overlap5,], gene + tissue_id ~ individual, value.var = "allelic_expression"))
ZZPU_13PLJ <- na.omit(dcast(suppl_table4_filtered[suppl_table4_filtered$gene %in% overlap3 & suppl_table4_filtered$tissue_id %in% overlap6,], gene + tissue_id ~ individual, value.var = "allelic_expression"))

UPIC_ZZPU$overlaping_between <- "UPIC vs ZZPU"
UPIC_13PLJ$overlaping_between <- "UPIC vs 13PLJ"
ZZPU_13PLJ$overlaping_between <- "ZZPU vs 13PLJ"

colnames(UPIC_ZZPU) <- c("gene", "tissue_id", "participant1", "participant2", "overlaping_between")
colnames(UPIC_13PLJ) <- c("gene", "tissue_id", "participant1", "participant2", "overlaping_between")
colnames(ZZPU_13PLJ) <- c("gene", "tissue_id", "participant1", "participant2", "overlaping_between")

df_overlaps <- rbind(UPIC_ZZPU, UPIC_13PLJ, ZZPU_13PLJ)

cor_test(data = df_overlaps, participant1, participant2, method = "spearman")

df_overlaps_with_class <- merge(df_overlaps, suppl_table7, by = "gene")


ggsave2(filename = paste(plot_dir, "Supplementary_Fig_1C.pdf"),
ggplot(df_overlaps_with_class, aes(x=participant1, y=participant2)) + 
  geom_point(aes(shape = category,col = overlaping_between)) + 
  geom_smooth(method = "lm", level=0.95) +
  stat_cor(method = "spearman") + 
  theme_AL_box()
)

####### Suppl figure 1d #######
# Are the hetSNPs the same?
UPIC_UPIC_ZZPU_hetSNPs <- suppl_table4_filtered[suppl_table4_filtered$gene %in% overlap1 & suppl_table4_filtered$tissue_id %in% overlap4 & suppl_table4_filtered$keep == T & suppl_table4_filtered$individual == "UPIC",] %>% dplyr::select(gene, position, tissue_id, individual, allelic_expression)
ZZPU_UPIC_ZZPU_hetSNPs <- suppl_table4_filtered[suppl_table4_filtered$gene %in% overlap1 & suppl_table4_filtered$tissue_id %in% overlap4 & suppl_table4_filtered$keep == T & suppl_table4_filtered$individual == "ZZPU",] %>% dplyr::select(gene, position, tissue_id, individual, allelic_expression)
smash_UPIC_ZZPU_hetSNPs <- merge(UPIC_UPIC_ZZPU_hetSNPs, ZZPU_UPIC_ZZPU_hetSNPs, by = c("gene", "position", "tissue_id"), all = T)
smash_UPIC_ZZPU_hetSNPs$tessst <- ifelse(is.na(smash_UPIC_ZZPU_hetSNPs$individual.x) | is.na(smash_UPIC_ZZPU_hetSNPs$allelic_expression.x) | is.na(smash_UPIC_ZZPU_hetSNPs$individual.y) | is.na(smash_UPIC_ZZPU_hetSNPs$allelic_expression.y), yes = "unique", no = "shared")


UPIC_UPIC_13PLJ_hetSNPs <- suppl_table4_filtered[suppl_table4_filtered$gene %in% overlap2 & suppl_table4_filtered$tissue_id %in% overlap5 & suppl_table4_filtered$keep == T & suppl_table4_filtered$individual == "UPIC",] %>% dplyr::select(gene, position, tissue_id, individual, allelic_expression)
ZZPU_UPIC_13PLJ_hetSNPs <- suppl_table4_filtered[suppl_table4_filtered$gene %in% overlap2 & suppl_table4_filtered$tissue_id %in% overlap5 & suppl_table4_filtered$keep == T & suppl_table4_filtered$individual == "13PLJ",] %>% dplyr::select(gene, position, tissue_id, individual, allelic_expression)
smash_UPIC_13PLJ_hetSNPs <- merge(UPIC_UPIC_13PLJ_hetSNPs, ZZPU_UPIC_13PLJ_hetSNPs, by = c("gene", "position", "tissue_id"), all = T)
smash_UPIC_13PLJ_hetSNPs$tessst <- ifelse(is.na(smash_UPIC_13PLJ_hetSNPs$individual.x) | is.na(smash_UPIC_13PLJ_hetSNPs$allelic_expression.x) | is.na(smash_UPIC_13PLJ_hetSNPs$individual.y) | is.na(smash_UPIC_13PLJ_hetSNPs$allelic_expression.y), yes = "unique", no = "shared")


UPIC_ZZPU_13PLJ_hetSNPs <- suppl_table4_filtered[suppl_table4_filtered$gene %in% overlap3 & suppl_table4_filtered$tissue_id %in% overlap6 & suppl_table4_filtered$keep == T & suppl_table4_filtered$individual == "ZZPU",] %>% dplyr::select(gene, position, tissue_id, individual, allelic_expression)
ZZPU_ZZPU_13PLJ_hetSNPs <- suppl_table4_filtered[suppl_table4_filtered$gene %in% overlap3 & suppl_table4_filtered$tissue_id %in% overlap6 & suppl_table4_filtered$keep == T & suppl_table4_filtered$individual == "13PLJ",] %>% dplyr::select(gene, position, tissue_id, individual, allelic_expression)
smash_ZZPU_13PLJ_hetSNPs <- merge(UPIC_ZZPU_13PLJ_hetSNPs, ZZPU_ZZPU_13PLJ_hetSNPs, by = c("gene", "position", "tissue_id"), all = T)
smash_ZZPU_13PLJ_hetSNPs$tessst <- ifelse(is.na(smash_ZZPU_13PLJ_hetSNPs$individual.x) | is.na(smash_ZZPU_13PLJ_hetSNPs$allelic_expression.x) | is.na(smash_ZZPU_13PLJ_hetSNPs$individual.y) | is.na(smash_ZZPU_13PLJ_hetSNPs$allelic_expression.y), yes = "unique", no = "shared")

countz_UPIC_ZZPU <- smash_UPIC_ZZPU_hetSNPs %>% dplyr::group_by(tessst) %>% dplyr::summarise(count=n())
countz_UPIC_ZZPU$overlaping_between <- "UPIC vs ZZPU"
countz_UPIC_ZZPU$perc <- countz_UPIC_ZZPU[countz_UPIC_ZZPU$tessst == "shared",]$count / (countz_UPIC_ZZPU[countz_UPIC_ZZPU$tessst == "shared",]$count + countz_UPIC_ZZPU[countz_UPIC_ZZPU$tessst != "shared",]$count) * 100

countz_UPIC_13PLJ <- smash_UPIC_13PLJ_hetSNPs %>% dplyr::group_by(tessst) %>% dplyr::summarise(count=n())
countz_UPIC_13PLJ$overlaping_between <- "UPIC vs 13PLJ"
countz_UPIC_13PLJ$perc <- countz_UPIC_13PLJ[countz_UPIC_13PLJ$tessst == "shared",]$count / (countz_UPIC_13PLJ[countz_UPIC_13PLJ$tessst == "shared",]$count + countz_UPIC_13PLJ[countz_UPIC_13PLJ$tessst != "shared",]$count) * 100

countz_ZZPU_13PLJ <- smash_ZZPU_13PLJ_hetSNPs %>% dplyr::group_by(tessst) %>% dplyr::summarise(count=n())
countz_ZZPU_13PLJ$overlaping_between <- "ZZPU vs 13PLJ"
countz_ZZPU_13PLJ$perc <- countz_ZZPU_13PLJ[countz_ZZPU_13PLJ$tessst == "shared",]$count / (countz_ZZPU_13PLJ[countz_ZZPU_13PLJ$tessst == "shared",]$count + countz_ZZPU_13PLJ[countz_ZZPU_13PLJ$tessst != "shared",]$count) * 100

ggsave2(filename = paste0(plot_dir, "Supplementary_Fig_1D.pdf"),
ggplot(data=rbind(countz_UPIC_ZZPU, countz_UPIC_13PLJ, countz_ZZPU_13PLJ), aes(x=overlaping_between, y=count, fill=tessst)) + 
  geom_bar(stat="identity", position = "dodge") + 
  geom_text(y=900, aes(label = round(perc, 2))) + 
  geom_text(aes(label=count)) + 
  theme_AL_box()
)