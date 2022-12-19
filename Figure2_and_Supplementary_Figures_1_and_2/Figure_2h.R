setwd("work_wd")
library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(VennDiagram)

source("r_sources/AL_ggplot_themes.R")
source("r_sources/plot_parameters.R")

##### Read in data ##### 

PAR <- rbind(fread("annotation_data/group-715_PAR1.csv"),
             fread("annotation_data/group-716_PAR2.csv"))


# Supplementary table 4
suppl_table4 <- fread("supplementary_tables/Supplementary_Table4.tsv")

# keep only highest covered hSNP.
suppl_table4_filtered <- data.table(suppl_table4[suppl_table4$keep == T,])

# read in our classifications
suppl_table6 <- fread("supplementary_tables/Supplementary_Table6.tsv")

# tuki df
tukiainen <- fread("annotation_data/landscape.Suppl.Table.13.csv", header = T)



PAR_our = unique(suppl_table4[suppl_table4$gene %in% suppl_table6[suppl_table6$category == "PAR",]$gene,]$gene)
Esc_our = unique(suppl_table4[suppl_table4$gene %in% suppl_table6[suppl_table6$category == "biallelic",]$gene,]$gene)
Variable_our = unique(suppl_table4[suppl_table4$gene %in% suppl_table6[suppl_table6$category == "variable",]$gene,]$gene)
Inactive_our = unique(suppl_table4[suppl_table4$gene %in% suppl_table6[suppl_table6$category == "monoallelic",]$gene,]$gene)
PAR_tukiainen = unique(tukiainen[tukiainen$Gene_name %in% PAR$`Approved symbol` & tukiainen$Gene_name %in% suppl_table4$gene,]$Gene_name)
Esc_tukiainen = unique(tukiainen[tukiainen$Reported_XCI_status == "Escape" & tukiainen$Gene_name %in% suppl_table4$gene,]$Gene_name)
Variable_tukiainen = unique(tukiainen[tukiainen$Reported_XCI_status == "Variable" & tukiainen$Gene_name %in% suppl_table4$gene,]$Gene_name)
Inactive_tukiainen = unique(tukiainen[tukiainen$Reported_XCI_status == "Inactive" & tukiainen$Gene_name %in% suppl_table4$gene,]$Gene_name)

# Chart
fig2X_PAR <- venn.diagram(
  x = list(PAR_our, PAR_tukiainen),
  category.names = c("Our" , "Tukiainen"),
  filename = NULL,
  imagetype="tiff",
  main = "PAR",
  fill = c("red", "blue"), alpha = c(0.33, 0.33), cex = 1,cat.fontface = 2,
  lty =2,
  disable.logging = T
)

fig2X_escape <- venn.diagram(
  x = list(Esc_our, Esc_tukiainen),
  category.names = c("Our" , "Tukiainen"),
  filename = NULL,
  imagetype="tiff",
  main = "Escape",
  fill = c("red", "blue"), alpha = c(0.33, 0.33), cex = 1,cat.fontface = 2,
  lty =2,
  disable.logging = T
)

fig2X_inacs <- venn.diagram(
  x = list(Inactive_our, Inactive_tukiainen),
  category.names = c("Our" , "Tukiainen"),
  filename = NULL,
  imagetype="tiff",
  main = "Inactive",
  fill = c("red", "blue"), alpha = c(0.33, 0.33), cex = 1,cat.fontface = 2,
  lty =2,
  disable.logging = T
)

fig2X_variable <- venn.diagram(
  x = list(Variable_our, Variable_tukiainen),
  category.names = c("Our" , "Tukiainen"),
  filename = NULL,
  imagetype="tiff",
  main = "Variable",
  fill = c("red", "blue"), alpha = c(0.33, 0.33), cex = 1,cat.fontface = 2,
  lty =2,
  disable.logging = T
)


ggsave2(filename = paste0(plot_dir, "Fig2h.pdf"), height = 8,  width=12,
        plot_grid(fig2X_PAR, fig2X_escape, fig2X_inacs, fig2X_variable)
)
