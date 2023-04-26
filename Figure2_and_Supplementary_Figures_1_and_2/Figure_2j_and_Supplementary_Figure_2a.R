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

# Annotation data
# PAR genes
PAR <- rbind(fread("project/Genome_biology_submission/annotation_data/group-715_PAR1.csv"),
             fread("project/Genome_biology_submission/annotation_data/group-716_PAR2.csv"))

# tuki df
tukiainen <- fread("project/Genome_biology_submission/annotation_data/landscape.Suppl.Table.13.csv", header = T)


# Supplementary table 4
suppl_table4 <- fread("project/Genome_biology_submission/supplementary_tables/Supplementary_Table4.tsv")

# Supplementary table 6
suppl_table5 <- fread("project/Genome_biology_submission/supplementary_tables/Supplementary_Table5.tsv")


# keep only highest covered hSNP.
suppl_table4_filtered <- data.table(suppl_table4[suppl_table4$keep == T,])

oooo <- suppl_table4_filtered[suppl_table4_filtered$gene %in% suppl_table5[suppl_table5$category %in% c("variable"),]$gene,]

oooo_short <- dplyr::select(oooo, individual, gene, refCount, altCount, tissue_id, position)

oooo_short_long = melt(oooo_short, id.vars = c("individual", "gene", "tissue_id", "position"),
                       measure.vars = c("refCount", "altCount"))

oooo_short_long_stats <- oooo_short_long %>% dplyr::group_by(gene, tissue_id, variable) %>% dplyr::summarise(meanz = mean(value))

urderz <- c("ANOS1", "ARSD", "CA5B", "DIP2KB", "MED14", "NAA10", "PRKX", "SMC1A",  "TMEM187", "TRAPPC2", "UBA1")

order_fin <- c(urderz, unique(oooo_short_long_stats[!oooo_short_long_stats$gene %in% urderz,]$gene))

# CURRENTLY USED
ggsave2(filename = paste0(plot_dir, "Supplementary_fig_2A.pdf"), height = 240, width = 180, limitsize = F, units = "mm", 
        plot_grid(
ggplot(oooo_short_long_stats, aes(x=tissue_id, y=meanz, fill = factor(variable, levels = c("refCount", "altCount")))) + 
  geom_bar(position = "stack", stat = "identity") + 
  facet_wrap(factor(gene, levels = order_fin)~., scales = "free_y", ncol = 3) +
  theme_AL_box_rotX() + 
  theme(axis.text.x = element_text(size=5), 
        axis.text.y = element_text(size=5), 
        legend.position = "none", 
        legend.text = element_text(size=5), 
        strip.text = element_text(size=5, margin=margin(t=0, b=0)),
        strip.background = element_blank()) +
  labs(x="", y="")+
  scale_fill_manual(values=c("refCount" = "#F15A22", "altCount" = "#6DC8BF")) +
  scale_y_continuous(limits = c(0,NA))+
  geom_hline(yintercept = 8, lty=2)
))



suppl_table4_filtered <- data.table(suppl_table4[suppl_table4$keep == T,])

novel <- c("ARSD","MED14")
known_true <- c("ANOS1","NAA10")
known_false <- c("USP11", "MAP7D3")

oooo2 <- suppl_table4_filtered[suppl_table4_filtered$gene %in% c(novel, known_true, known_false),]

oooo_short2 <- dplyr::select(oooo2, individual, gene, refCount, altCount, tissue_id, position)

oooo_short_long2 = melt(oooo_short2, id.vars = c("individual", "gene", "tissue_id", "position"),
                        measure.vars = c("refCount", "altCount"))

oooo_short_long_stats2 <- oooo_short_long2 %>% dplyr::group_by(gene, tissue_id, variable) %>% dplyr::summarise(meanz = mean(value))


ggsave2(filename = paste0(plot_dir, "Figure_2J.pdf"), height = 85, width = 75, limitsize = F, units = "mm", 
        plot_grid(
          ggplot(oooo_short_long_stats2[oooo_short_long_stats2$gene %in% c(novel, known_true, known_false),], aes(x=tissue_id, y=meanz, fill = factor(variable, levels = c("refCount", "altCount", "effectSize")))) + 
            geom_bar(position = "stack", stat = "identity") + 
            facet_wrap(factor(gene, levels = c(novel, known_true, known_false))~., scales = "free_y", ncol = 1) +
            theme_AL_box_rotX() + 
            theme(axis.text.x = element_text(size=5), 
                  axis.text.y = element_text(size=5), 
                  legend.position = "none", 
                  legend.text = element_text(size=5), 
                  strip.text = element_text(size=5, margin=margin(t=0, b=0)),
                  strip.background = element_blank()) +
            labs(x="", y="")+
            scale_fill_manual(values=c("refCount" = "#F15A22", "altCount" = "#6DC8BF")) +
            scale_y_continuous(limits = c(0,NA))+
            geom_hline(yintercept = 8, lty=2)
          )
        )



