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


# Supplementary table 5
suppl_table5 <- fread("project/Genome_biology_submission/supplementary_tables/Supplementary_Table5.tsv")

suppl_table6 <- fread("project/Genome_biology_submission/supplementary_tables/Supplementary_Table6.tsv")


# TPM figure, our classifications
df_TPM_ours <- merge(suppl_table6, dplyr::select(suppl_table5, gene, category), by = "gene")

df_TPM_ours[df_TPM_ours$category == "monoallelic",]$category <- "inactive"
df_TPM_ours[df_TPM_ours$category == "biallelic",]$category <- "escape"

df_TPM_ours_keep <- df_TPM_ours %>% dplyr::group_by(tissue_id, gene) %>% dplyr::summarise(meanz=mean(value))
df_TPM_ours_filt <- merge(df_TPM_ours, df_TPM_ours_keep, by = c("gene", "tissue_id"))
df_TPM_ours_filt <- df_TPM_ours[df_TPM_ours$value > 1,]

df_TPM_ours_stats2 <- df_TPM_ours_filt %>% dplyr::group_by(category, tissue_id) %>% rstatix::get_summary_stats(value, type = "common")

stat.test_our <- compare_means(ref.group = "variable", value ~ category,  data = df_TPM_ours_filt, method = "wilcox.test",p.adjust.method = "BH")
stat.test_our$y.position <- c(200,250,275)


# TPM figure, Tukiainens classifications
df_TPM_tuki <- merge(suppl_table6, dplyr::select(tukiainen[tukiainen$Reported_XCI_status != "Unknown",], Gene_name, Reported_XCI_status), by.x = "gene", by.y = "Gene_name")

df_TPM_tuki[df_TPM_tuki$Reported_XCI_status == "Variable",]$Reported_XCI_status <- "variable"
df_TPM_tuki[df_TPM_tuki$Reported_XCI_status == "Escape",]$Reported_XCI_status <- "escape"
df_TPM_tuki[df_TPM_tuki$Reported_XCI_status == "Inactive",]$Reported_XCI_status <- "inactive"
df_TPM_tuki[df_TPM_tuki$gene %in% PAR$`Approved symbol`,]$Reported_XCI_status <- "PAR"

df_TPM_tuki_keep <- df_TPM_tuki %>% dplyr::group_by(tissue_id, gene) %>% dplyr::summarise(meanz=mean(value))
df_TPM_tuki_filt <- merge(df_TPM_tuki, df_TPM_tuki_keep, by = c("gene", "tissue_id"))

df_TPM_tuki_filt <- df_TPM_tuki_filt[df_TPM_tuki_filt$meanz > 1,]

df_TPM_tuki_stats2 <- df_TPM_tuki_filt %>% dplyr::group_by(Reported_XCI_status, tissue_id) %>% rstatix::get_summary_stats(value, type = "common")

stat.test_tuki <- compare_means(ref.group = "variable", value ~ Reported_XCI_status,  data = df_TPM_tuki_filt, method = "wilcox.test",p.adjust.method = "BH")
stat.test_tuki$y.position <- c(70,85,100)



ggsave2(filename = paste0(plot_dir, "Supplementary_Fig1e.pdf"), width=45, height = 100, unit = "mm",
        plot_grid(ncol=1,
          plot_grid(ggplot(df_TPM_ours_stats2, 
                           aes(x=factor(category, levels = c("PAR","escape","inactive","variable" )), y=mean)) + 
                      geom_violin() +
                      geom_jitter(aes(col=tissue_id, shape = tissue_id), width = 0.25, size = 0.5) + 
                      scale_color_manual(values=color_valzz, limits = force) + 
                      scale_shape_manual(values=shape_valzz, limits = force) + 
                      labs(y="", x="") + 
                      theme_AL_box_rotX() + 
                      #coord_cartesian(ylim=c(0,91)) +
                      #scale_y_continuous(breaks = c(0,20,40,60,80))+
                      stat_pvalue_manual(stat.test_our, label = "p.signif", size = 3.25)+
                      theme(legend.position = "none",
                            axis.text.x = element_text(size=6),
                            axis.text.y = element_text(size=6))),
          
          plot_grid(ggplot(df_TPM_tuki_stats2, 
                           aes(x=factor(Reported_XCI_status, levels = c("PAR","escape","inactive", "variable" )), y=mean)) + 
                      geom_violin() +
                      geom_jitter(aes(col=tissue_id, shape = tissue_id), width = 0.25, size = 0.5) + 
                      scale_color_manual(values=color_valzz, limits = force) + 
                      scale_shape_manual(values=shape_valzz, limits = force) + 
                      labs(y="", x="") + 
                      theme_AL_box_rotX() + 
                      coord_cartesian(ylim=c(0,NA)) +
                      #scale_y_continuous(breaks = c(0,20,40,60,80))+
                      stat_pvalue_manual(stat.test_tuki, label = "p.signif", size = 3.25)+
                      theme(legend.position = "none",
                            axis.text.x = element_text(size=6),
                            axis.text.y = element_text(size=6))))
)


new_tissue_cats <- data.frame(c("brain-amyg" = "brain",
                                "brain-ante" ="brain" ,
                                "brain-caud" ="brain",
                                "brain-cehe" = "brain", 
                                "brain-cere" = "brain" ,
                                "brain-cort" = "brain" ,
                                "brain-frco" = "brain" ,
                                "brain-hipp" = "brain" ,
                                "brain-hypo" = "brain",
                                "brain-nucl" = "brain",
                                "brain-puta" = "brain",
                                "brain-spin" = "brain",
                                "brain-subs" ="brain",
                                "esophagus-gaju" = "esophagus" ,
                                "esophagus-muco" ="esophagus",
                                "esophagus-musc" = "esophagus",
                                "skin-lleg" ="skin",
                                "skin-supr" = "skin",
                                "heart-atri" = "heart",
                                "heart-vent" = "heart",
                                "adipose-subc" = "adipose",
                                "adipose-visc" ="adipose",
                                "artery-aort" = "artery",
                                "artery-coro" = "artery" ,
                                "artery-tibi" = "artery",
                                "adrenal gland" = "adrenal gland",
                                "kidney-cort" = "kidney",
                                "kidney-medu" = "kidney",
                                "bladder" = "bladder",
                                "nerve" = "nerve",
                                "salivary gland" = "salivary gland",
                                "thyroid" = "thyroid",
                                "pituitary" = "pituitary",
                                "pancreas" = "pancreas",
                                "spleen" = "spleen",
                                "liver" = "liver",
                                "lung" = "lung",
                                "breast" ="breast",
                                "muscle" = "muscle",
                                "small intestine" = "small intestine",
                                "colon-sigm" = "colon",
                                "colon-tran" = "colon",
                                "stomach" ="stomach",
                                "lymphocytes" = "lymphocytes",
                                "fibroblasts" = "fibroblasts",
                                "whole blood" = "whole blood",
                                "cervix-ecto" = "cervix",
                                "cervix-endo" = "cervix",
                                "uterus" = "uterus",
                                "vagina" = "vagina",
                                "ovary" = "ovary",
                                "fallopian tube" = "fallopian tube"))

colnames(new_tissue_cats) <- c("new_cat")
new_tissue_cats$tissue_id <- rownames(new_tissue_cats)

df_TPM_ours_filt_new_cat <- merge(df_TPM_ours_filt, new_tissue_cats, by = "tissue_id")
df_TPM_tuki_filt_new_cat <- merge(df_TPM_tuki_filt, new_tissue_cats, by = "tissue_id")

df_TPM_ours_filt_new_cat_stats2 <- df_TPM_ours_filt_new_cat %>% dplyr::group_by(category, new_cat) %>% rstatix::get_summary_stats(value, type = "common")
df_TPM_tuki_filt_new_cat_stats2 <- df_TPM_tuki_filt_new_cat %>% dplyr::group_by(Reported_XCI_status, new_cat) %>% rstatix::get_summary_stats(value, type = "common")



ggsave2(filename = paste0(plot_dir, "Fig2h.pdf"), width=4, height = 3,
        plot_grid(
          plot_grid(ggplot(df_TPM_ours_filt_new_cat_stats2, 
                           aes(x=factor(category, levels = c("PAR","escape","inactive","variable" )), y=mean)) + 
                      geom_violin() +
                      geom_jitter(aes(col=new_cat), width = 0.25, size = 1) + 
                      scale_color_manual(values=new_tissue_cats_color_valzz, limits = force)+
                      labs(y="", x="") + 
                      theme_AL_box_rotX() + 
                      #coord_cartesian(ylim=c(0,91)) +
                      #scale_y_continuous(breaks = c(0,20,40,60,80))+
                      stat_pvalue_manual(stat.test_our, label = "p.signif", size = 3.25)+
                      theme(legend.position = "none",
                            axis.text.x = element_text(size=6),
                            axis.text.y = element_text(size=6))),
          
          plot_grid(ggplot(df_TPM_tuki_filt_new_cat_stats2, 
                           aes(x=factor(Reported_XCI_status, levels = c("PAR","escape","inactive", "variable" )), y=mean)) + 
                      geom_violin() +
                      geom_jitter(aes(col=new_cat), width = 0.25, size = 1) + 
                      scale_color_manual(values=new_tissue_cats_color_valzz, limits = force)+
                      labs(y="", x="") + 
                      theme_AL_box_rotX() + 
                      coord_cartesian(ylim=c(0,NA)) +
                      #scale_y_continuous(breaks = c(0,20,40,60,80))+
                      stat_pvalue_manual(stat.test_tuki, label = "p.signif", size = 3.25)+
                      theme(legend.position = "none",
                            axis.text.x = element_text(size=6),
                            axis.text.y = element_text(size=6))))
)

countz <- unique(dplyr::select(df_TPM_ours_filt, participant, tissue_id)) %>% dplyr::group_by(tissue_id) %>% dplyr::count()

ggsave2(filename = paste0(plot_dir, "Supplementary_Figure_1f.pdf"), width=4, height = 3,
        ggplot(data=countz, aes(x=reorder(tissue_id, n), y=n)) + geom_col() + theme_AL_box_rotX() + geom_text(aes(label=n), size = 3)
)
