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
library(ggalluvial)

# read in ggplot theemes (made by Antonio Lentini) and plot parameters.
source("r_sources/AL_ggplot_themes.R")
source("r_sources/plot_parameters.R")

##### Read in data ##### 

# Annotation data

# Annotation data
# PAR genes
PAR <- rbind(fread("annotation_data/group-715_PAR1.csv"),
             fread("annotation_data/group-716_PAR2.csv"))

# tuki df
tukiainen <- fread("annotation_data/landscape.Suppl.Table.13.csv", header = T)


# Supplementary table 6
suppl_table6 <- fread("supplementary_tables/Supplementary_Table6.tsv")


##### Alluvial ######
df_classification <- dplyr::select(tukiainen[tukiainen$Gene_name %in% suppl_table6$gene,], Gene_name, Reported_XCI_status)
colnames(df_classification) <- c("gene", "category")
df_classification$axes <- "Tukiainen"
df_classification <- df_classification[!(df_classification$gene == "IDS" & df_classification$category == "Unknown"),]

df_gylemo <- dplyr::select(suppl_table6, gene, category)
df_gylemo$axes <- "Gylemo"
df_gylemo[df_gylemo$category == "monoallelic",]$category <- "Inactive"
df_gylemo[df_gylemo$category == "biallelic",]$category <- "Escape"
df_gylemo[df_gylemo$category == "variable",]$category <- "Variable"


novelz <- data.frame(gene = setdiff(df_gylemo$gene, df_classification$gene),
                     category = "not_investigated")
colnames(novelz) <- c("gene", "category")
novelz$axes <- "Tukiainen"

df_classification <- rbind(df_classification, novelz)

smash_fin <- rbind(df_classification, df_gylemo)


#smash_fin <- unique(smash_fin[smash_fin$gene %in% PAR$`Approved symbol`,])

is_lodes_form(smash_fin, key = "axes", value = "category", id = "gene")

smash_fin[smash_fin$gene %in% c("ANOS1", "PRKX", "ARSD", "MED14", "NAA10", "CA5B", "SMC1A", "DIP2KB", "TMEM187", "TRAPPC2", "UBA1") & smash_fin$axes == "Gylemo",]$category <- "curated variable"
smash_fin[smash_fin$category == "Variable" & smash_fin$axes == "Gylemo",]$category <- "Undeterminable"
smash_fin[is.na(smash_fin$category)& smash_fin$axes == "Gylemo",]$category <- "PAR"
smash_fin[smash_fin$axes == "Tukiainen" & smash_fin$gene %in% tukiainen[tukiainen$Region == "PAR"]$Gene_name,]$category <- "PAR"

smash_fin$category <- factor(smash_fin$category, levels = c("PAR","not_investigated", "Escape", "Variable", "Unknown", "curated variable", "Undeterminable", "Inactive"))



smash_fin %>% dplyr::group_by(category, axes) %>% dplyr::count()

smash_fin[smash_fin$axes == "Tukiainen",]$axes <- "Previous"
smash_fin$axes <- factor(smash_fin$axes, levels = c("Previous", "Gylemo"))

ggsave2(filename = paste0(plot_dir, "Figure_2k.pdf"), width=3, height = 4,
        plot_grid(ggplot(smash_fin, aes(alluvium = gene, x = axes, stratum = category)) + 
                    geom_flow(aes(fill=category)) +
                    geom_stratum(aes(fill=category))  + 
                    theme_AL_simple(legend.title = element_blank())+
                    theme(axis.text.x = element_text(size=6),
                          axis.text.y = element_text(size=6),
                          legend.text = element_text(size=6))+
                    labs(x="",y="gene count")+
                    scale_fill_manual(values = c("PAR" = "yellow", "Escape" = "red", "Inactive" = "#869a9a", "Variable" = "purple", "Unknown" = "blue", "Undeterminable" = "green", "curated variable" = "pink", "not_investigated" = "black"))
        )
)

