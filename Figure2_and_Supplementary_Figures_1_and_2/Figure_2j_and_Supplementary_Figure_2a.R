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
# PAR genes
PAR <- rbind(fread("annotation_data/group-715_PAR1.csv"),
             fread("annotation_data/group-716_PAR2.csv"))

# tuki df
tukiainen <- fread("annotation_data/landscape.Suppl.Table.13.csv", header = T)


# Supplementary table 4
suppl_table4 <- fread("supplementary_tables/Supplementary_Table4.tsv")

# Supplementary table 6
suppl_table6 <- fread("supplementary_tables/Supplementary_Table6.tsv")


# keep only highest covered hSNP.
suppl_table4_filtered <- data.table(suppl_table4[suppl_table4$keep == T,])

oooo <- suppl_table4_filtered[suppl_table4_filtered$gene %in% suppl_table6[suppl_table6$category == "variable",]$gene,]

oooo_short <- dplyr::select(oooo, individual, gene, refCount, altCount, tissue_id, position, allelic_expression)

oooo_short_long = melt(oooo_short, id.vars = c("individual", "gene", "tissue_id", "position"),
                       measure.vars = c("refCount", "altCount", "allelic_expression"))

oooo_short_long$splitter <- ifelse(oooo_short_long$variable == "allelic_expression", yes = "AE", no = "Counts")




oooo_short_long_stats <- oooo_short_long %>% dplyr::group_by(gene, tissue_id, variable, splitter) %>% dplyr::summarise(meanz = mean(value))

urderz <- c("ANOS1", "ARSD", "CA5B", "DIP2KB", "MED14", "NAA10", "PRKX", "SMC1A",  "TMEM187", "TRAPPC2", "UBA1")

order_fin <- c(urderz, unique(oooo_short_long_stats[!oooo_short_long_stats$gene %in% urderz,]$gene))

ggsave2(filename = paste0(plot_dir, "Supplementary_figure_2a_part1.pdf"), height = 400, width = 300, limitsize = F, units = "mm", 
        plot_grid(
    ggplot(oooo_short_long_stats, aes(x=tissue_id, y=meanz, fill = factor(variable, levels = c("refCount", "altCount", "effectSize")))) + 
      geom_col(position = "dodge", width = .75) + 
      geom_point(col="black", size = 0.5)+ 
      facet_wrap(factor(gene, levels = order_fin)~splitter, scales = "free_y", ncol = 6) +
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
      geom_hline(yintercept = 0.5, lty=2, color = "black")+
      geom_hline(yintercept = 0.4, lty=2, color = "red")
    
))

ggsave2(filename = paste0(plot_dir, "Supplementary_figure_2a_part2.pdf"), height = 400, width = 300, limitsize = F, units = "mm", 
        plot_grid(
          ggplot(oooo_short_long_stats, aes(x=tissue_id, y=meanz, fill = factor(variable, levels = c("refCount", "altCount", "effectSize")))) + 
            geom_col(position = "dodge", width = .75) + 
            geom_point(col="black", size = 0.5)+ 
            facet_wrap(factor(gene, levels = order_fin)~splitter, scales = "free_y", ncol = 6) +
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
            geom_hline(yintercept = 0.5, lty=2, color = "black")+
            geom_hline(yintercept = 8, lty=2, color = "blue")
          
        ))







suppl_table4_filtered <- data.table(suppl_table4[suppl_table4$keep == T,])

novel <- c("ARSD","MED14")
known_true <- c("ANOS1","NAA10")
known_false <- c("USP11", "MAP7D3")

oooo2 <- suppl_table4_filtered[suppl_table4_filtered$gene %in% c(novel, known_true, known_false),]

oooo_short2 <- dplyr::select(oooo2, individual, gene, refCount, altCount, tissue_id, position, allelic_expression)

oooo_short_long2 = melt(oooo_short2, id.vars = c("individual", "gene", "tissue_id", "position"),
                       measure.vars = c("refCount", "altCount", "allelic_expression"))

oooo_short_long2$splitter <- ifelse(oooo_short_long2$variable == "allelic_expression", yes = "AE", no = "Counts")




oooo_short_long_stats2 <- oooo_short_long2 %>% dplyr::group_by(gene, tissue_id, variable, splitter) %>% dplyr::summarise(meanz = mean(value))




ggsave2(filename = paste0(plot_dir, "Figure_2j_part1.pdf"), height = 100, width = 150, limitsize = F, units = "mm", 
        plot_grid(
          ggplot(oooo_short_long_stats2[oooo_short_long_stats2$gene %in% c(novel, known_true, known_false),], aes(x=tissue_id, y=meanz, fill = factor(variable, levels = c("refCount", "altCount", "effectSize")))) + 
            geom_col(position = "dodge", width = .75) + 
            geom_point(col="black", size = 0.5)+ 
            facet_wrap(factor(gene, levels = c(novel, known_true, known_false))~splitter, scales = "free_y", ncol = 2) +
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
            geom_hline(yintercept = 0.5, lty=2, color = "black")+
            geom_hline(yintercept = 0.4, lty=2, color = "red")
          
        ))

ggsave2(filename = paste0(plot_dir, "Figure_2j_part2.pdf"), height = 100, width = 150, limitsize = F, units = "mm", 
        plot_grid(
          ggplot(oooo_short_long_stats2[oooo_short_long_stats2$gene %in% c(novel, known_true, known_false),], aes(x=tissue_id, y=meanz, fill = factor(variable, levels = c("refCount", "altCount", "effectSize")))) + 
            geom_col(position = "dodge", width = .75) + 
            geom_point(col="black", size = 0.5)+ 
            facet_wrap(factor(gene, levels = c(novel, known_true, known_false))~splitter, scales = "free_y", ncol = 2) +
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
            geom_hline(yintercept = 0.5, lty=2, color = "black")+
            geom_hline(yintercept = 8, lty=2, color = "blue")
          
        ))
