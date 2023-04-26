library(data.table)
library(cowplot)
library(VennDiagram)
source("project/Genome_biology_submission/r_sources/plot_parameters.R")

##### Read in data ##### 
# Supplementary table 4
suppl_table4 <- fread("project/Genome_biology_submission/supplementary_tables/Supplementary_Table4.tsv")

# keep only highest covered hSNP.
suppl_table4_filtered <- data.table(suppl_table4[suppl_table4$keep == T,])


####### figure 2a and figure 2b ######
# Generate 3 sets
set1_snps <- unique(suppl_table4_filtered[suppl_table4_filtered$individual == "UPIC",]$position)
set2_snps <- unique(suppl_table4_filtered[suppl_table4_filtered$individual == "ZZPU",]$position)
set3_snps <- unique(suppl_table4_filtered[suppl_table4_filtered$individual == "13PLJ",]$position)

# Chart
fig2a <- venn.diagram(
  x = list(set1_snps, set2_snps, set3_snps),
  category.names = c("UPIC" , "ZZPU " , "13PLJ"),
  filename = NULL,
  imagetype="tiff",
  main = paste("All hSNPs =", length(unique(suppl_table4_filtered$position))),
  fill = c("red", "blue", "yellow"), alpha = c(0.33, 0.33, 0.33), cex = 1, cat.fontface = 2,
  lty =2
)


set1_genes <- unique(suppl_table4_filtered[suppl_table4_filtered$individual == "UPIC",]$gene)
set2_genes <- unique(suppl_table4_filtered[suppl_table4_filtered$individual == "ZZPU",]$gene)
set3_genes <- unique(suppl_table4_filtered[suppl_table4_filtered$individual == "13PLJ",]$gene)

# Chart
fig2b <- venn.diagram(
  x = list(set1_genes, set2_genes, set3_genes),
  category.names = c("UPIC" , "ZZPU " , "13PLJ"),
  filename = NULL,
  imagetype="tiff",
  main = paste("All genes =", length(unique(suppl_table4_filtered$gene))),
  fill = c("red", "blue", "yellow"), alpha = c(0.33, 0.33, 0.33), cex = 1,cat.fontface = 2,
  lty =2
)



ggsave2(filename = paste0(plot_dir, "Fig2a_2b.pdf"), height = 4,  width=8,
        plot_grid(fig2a, fig2b)
)