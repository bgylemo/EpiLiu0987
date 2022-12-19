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

# Supplementary table 4
suppl_table4 <- fread("supplementary_tables/Supplementary_Table4.tsv")
# keep only highest covered hSNP.
suppl_table4_filtered <- data.table(suppl_table4[suppl_table4$keep == T,])

########## fig 2d,2e,2f and 2g (ploting AE for individual genes across tissues) ###############

goiz <- c("AKAP17A", "GTPBP6", "DDX3X", "PUDP", "ACOT9", "IDS", "ARSD", "ANOS1")

ggsave2(filename = paste0(plot_dir, "figure_2d_e_f_g.pdf"), height = 60, width = 180, limitsize = F, units = "mm", 
        ggplot(data=suppl_table4_filtered[suppl_table4_filtered$gene %in% goiz,], aes(x=tissue_id, y=allelic_expression)) + 
          geom_point(aes(col=individual)) +
          facet_wrap(factor(gene, levels = goiz)~., ncol = 2) + 
          theme_AL_box_rotX() + 
          theme(axis.text.x = element_text(size=6), 
                axis.text.y = element_text(size=6), 
                legend.position = "none", 
                legend.text = element_text(size=6), 
                strip.text.x = element_text(size=6, margin=margin(t=0, b=0))) +
          geom_hline(yintercept = 0.4, lty=2, size = 0.5)+
          scale_y_continuous(limits = c(0,0.5), breaks = c(0,0.1,0.2,0.3,0.4,0.5))+
          labs(x="", y="")
)
