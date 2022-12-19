setwd("work_wd")
library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(rtracklayer)
library(GenomicRanges)
library(ggplotify)
library(karyoploteR)
# read in ggplot theemes (made by Antonio Lentini) and plot parameters.
source("r_sources/AL_ggplot_themes.R")
source("r_sources/plot_parameters.R")

##### Read in data ##### 

# Supplementary table 4
suppl_table4 <- fread("supplementary_tables/Supplementary_Table4.tsv")

suppl_table6 <- fread("supplementary_tables/Supplementary_Table6.tsv")

# read in gencode v41
gtf_raw <- import("annotation_data/gencode.v41.annotation.gtf")
gtf_raw <- gtf_raw[seqnames(gtf_raw) == "chrX",]




####### figure 2c ######
# make zoom
entire_X_zoom.region <- toGRanges(data.frame("chrX", 1, 160000000))

# keep only gene markers we are interested in
markerz_our <- gtf_raw[gtf_raw$gene_name %in% unique(suppl_table4$gene) & gtf_raw$type == "gene",]

# Add our gene categories
markerz_our_cat <- merge(markerz_our, suppl_table6, by.x = "gene_name", by.y = "gene")

markerz_our_cat_PAR <- markerz_our_cat[markerz_our_cat$start < 2781479 | markerz_our_cat$category == "PAR",]
markerz_our_cat_always_biallelic <- markerz_our_cat[markerz_our_cat$category == "biallelic" & markerz_our_cat$start > 2781479,]
markerz_our_cat_always_monoallelic <- markerz_our_cat[markerz_our_cat$category == "monoallelic" & markerz_our_cat$start > 2781479,]
markerz_our_cat_variable <- markerz_our_cat[markerz_our_cat$category == "variable" & markerz_our_cat$start > 2781479,]

plot_fin_alles <- as.ggplot(expression(
  kp <- plotKaryotype(plot.type=3, chromosomes = "chrX",  main = "gene overlap", genome = "hg38", zoom = entire_X_zoom.region),
  kpAddBaseNumbers(kp),
  kpPoints(kp, 
           data = markerz_our_cat[markerz_our_cat$gene_name == "XIST",],
           chr = markerz_our_cat[markerz_our_cat$gene_name == "XIST",]$seqnames,
           x=markerz_our_cat[markerz_our_cat$gene_name == "XIST",]$start,
           y=0.3,
           data.panel = 1,
           col="red"
  ),
  kpPoints(kp, 
           data = markerz_our_cat[markerz_our_cat$gene_name == "DDX3X",],
           chr = markerz_our_cat[markerz_our_cat$gene_name == "DDX3X",]$seqnames,
           x=markerz_our_cat[markerz_our_cat$gene_name == "DDX3X",]$start,
           y=0.3,
           data.panel = 1,
           col="blue"
  ),
  kpPoints(kp, 
           data = markerz_our_cat_PAR,
           chr = markerz_our_cat_PAR$seqnames,
           x=markerz_our_cat_PAR$start,
           y=0.1,
           data.panel = 1
  ),
  kpAxis(kp, data.panel = 1, r0 = 0.1, r1 = 0.1, cex=0.5,
         labels = "PAR"), 
  kpPoints(kp, 
           data = markerz_our_cat_always_biallelic,
           chr = markerz_our_cat_always_biallelic$seqnames,
           x=markerz_our_cat_always_biallelic$start,
           y=0.15,
           data.panel = 1
  ),
  kpAxis(kp, data.panel = 1, r0 = 0.15, r1 = 0.15, cex=0.5, 
         labels = "biallelic"),
  kpPoints(kp, 
           data = markerz_our_cat_always_monoallelic,
           chr = markerz_our_cat_always_monoallelic$seqnames,
           x=markerz_our_cat_always_monoallelic$start,
           y=0.2,
           data.panel = 1
  ),
  kpAxis(kp, data.panel = 1, r0 = 0.2, r1 = 0.2, cex=0.5, 
         labels = "monoallelic"),
  kpPoints(kp, 
           data = markerz_our_cat_variable,
           chr = markerz_our_cat_variable$seqnames,
           x=markerz_our_cat_variable$start,
           y=0.25,
           data.panel = 1
  ),
  kpAxis(kp, data.panel = 1, r0 = 0.25, r1 = 0.25, cex=0.5, 
         labels = "variable")
  )
  )

ggsave2(filename = paste0(plot_dir, "Fig_2c.pdf"), height = 14, width = 24, units = "cm",
        plot_grid(plot_fin_alles)
)


######## Suppl_Fig 2a ######

markerz <- gtf_raw[gtf_raw$gene_name %in% unique(suppl_table4$gene) & gtf_raw$type == "gene",]

# I want to mark the Tukiainen UPIC genes to highlight we have much higher gene coverage.
Tuki_UPIC <- fread("annotation_data/Suppl.Table.5.csv")

markerz$TUKI <- ifelse(markerz$gene_name %in% Tuki_UPIC$`Gene name`, yes = "in_TUKI", no = "unique")
markerz$TUKI <- factor(markerz$TUKI, levels = c("in_TUKI", "unique"))

entire_X_zoom.region <- toGRanges(data.frame("chrX", 1, 160000000))

plot_fin_allez_gene_labels <- as.ggplot(expression(
  kp <- plotKaryotype(plot.type=3, chromosomes = "chrX",  main = "gene overlap", genome = "hg38", zoom = entire_X_zoom.region),
  kpAddBaseNumbers(kp),
  kpPlotMarkers(kp, 
                data = markerz[1:31],
                labels=markerz[1:31]$gene_name, 
                label.color = colByCategory(markerz[1:31]$TUKI, colors=c("#FDB913", "#524FA1")),
                data.panel = 1,
                cex = 0.35,
                adjust.label.position = F,
                ignore.chromosome.ends = T
  ),
  kpPlotMarkers(kp, 
                data = markerz[32:63],
                labels=markerz[32:63]$gene_name, 
                label.color = colByCategory(markerz[32:63]$TUKI, colors=c("#FDB913", "#524FA1")),
                data.panel = 2, 
                cex = 0.35,
                adjust.label.position = F,
                ignore.chromosome.ends = T
  ),
  kpPlotMarkers(kp, 
                data = markerz[64:95],
                labels=markerz[64:95]$gene_name, 
                label.color = colByCategory(markerz[1:95]$TUKI, colors=c("#FDB913", "#524FA1")),
                data.panel = 1, 
                cex = 0.35,
                adjust.label.position = F,
                ignore.chromosome.ends = T
  ),
  kpPlotMarkers(kp, 
                data = markerz[96:127],
                labels=markerz[96:127]$gene_name, 
                label.color = colByCategory(markerz[96:127]$TUKI, colors=c("#FDB913", "#524FA1")),
                data.panel = 2, 
                cex = 0.35,
                adjust.label.position = F,
                ignore.chromosome.ends = T
  ),
  kpPlotMarkers(kp, 
                data = markerz[128:159],
                labels=markerz[128:159]$gene_name, 
                label.color = colByCategory(markerz[128:159]$TUKI, colors=c("#FDB913", "#524FA1")),
                data.panel = 1, 
                cex = 0.35,
                adjust.label.position = F,
                ignore.chromosome.ends = T
  ),
  kpPlotMarkers(kp, 
                data = markerz[160:191],
                labels=markerz[160:191]$gene_name, 
                label.color = colByCategory(markerz[169:191]$TUKI, colors=c("#FDB913", "#524FA1")),
                data.panel = 2, 
                cex = 0.35,
                adjust.label.position = F,
                ignore.chromosome.ends = T
  ),
  kpPlotMarkers(kp, 
                data = markerz[192:223],
                labels=markerz[192:223]$gene_name, 
                label.color = colByCategory(markerz[192:223]$TUKI, colors=c("#FDB913", "#524FA1")),
                data.panel = 1, 
                cex = 0.35,
                adjust.label.position = F,
                ignore.chromosome.ends = T
  ),
  kpPlotMarkers(kp, 
                data = markerz[224:249],
                labels=markerz[224:249]$gene_name, 
                label.color = colByCategory(markerz[224:249]$TUKI, colors=c("#FDB913", "#524FA1")),
                data.panel = 2, 
                cex = 0.35,
                adjust.label.position = F,
                ignore.chromosome.ends = T
  ),
  kpPlotMarkers(kp, 
                data = markerz[250:280],
                labels=markerz[250:280]$gene_name, 
                label.color = colByCategory(markerz[250:280]$TUKI, colors=c("#FDB913", "#524FA1")),
                data.panel = 1, 
                cex = 0.35,
                adjust.label.position = F,
                ignore.chromosome.ends = T
  ),
  kpPlotMarkers(kp, 
                data = markerz[281:311],
                labels=markerz[281:311]$gene_name, 
                label.color = colByCategory(markerz[281:311]$TUKI, colors=c("#FDB913", "#524FA1")),
                data.panel = 2, 
                cex = 0.35,
                adjust.label.position = F,
                ignore.chromosome.ends = T
  ),
  kpPlotMarkers(kp, 
                data = markerz[312:342],
                labels=markerz[312:342]$gene_name, 
                label.color = colByCategory(markerz[312:342]$TUKI, colors=c("#FDB913", "#524FA1")),
                data.panel = 1, 
                cex = 0.35,
                adjust.label.position = F,
                ignore.chromosome.ends = T
  ),
  kpPlotMarkers(kp, 
                data = markerz[343:373],
                labels=markerz[343:373]$gene_name, 
                label.color = colByCategory(markerz[343:373]$TUKI, colors=c("#FDB913", "#524FA1")),
                data.panel = 2, 
                cex = 0.35,
                adjust.label.position = F,
                ignore.chromosome.ends = T
  ),
  kpPlotMarkers(kp, 
                data = markerz[374:392],
                labels=markerz[374:392]$gene_name, 
                label.color = colByCategory(markerz[374:392]$TUKI, colors=c("#FDB913", "#524FA1")),
                data.panel = 1, 
                cex = 0.35,
                adjust.label.position = F,
                ignore.chromosome.ends = T
  ),
  kpPlotMarkers(kp, 
                data = markerz[393:401],
                labels=markerz[393:401]$gene_name, 
                label.color = colByCategory(markerz[393:401]$TUKI, colors=c("#FDB913", "#524FA1")),
                data.panel = 2, 
                cex = 0.35,
                adjust.label.position = F,
                ignore.chromosome.ends = T
  )
  )
  )


ggsave2(filename = paste0(plot_dir, "Suppl_Fig_2A.pdf"), height = 14, width = 24, units = "cm",
        plot_grid(plot_fin_allez_gene_labels)
)

