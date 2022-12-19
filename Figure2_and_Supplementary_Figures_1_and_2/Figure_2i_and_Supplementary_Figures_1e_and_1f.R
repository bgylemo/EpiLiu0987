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
suppl_table6 <- fread("supplementary_tables/Supplementary_Table6.tsv")


######### Figure 2i #############
df_TPM_tissues <- melt(fread("supplementary_data/22_11_07_expression_matrix_unfiltered.tsv"))
df_TPM_tissues$variable <- gsub(df_TPM_tissues$variable, pattern = "\\.", replacement = "-")
meta_full <- fread("annotation_data/metadata_RNAseq_full_for_tissues.tsv")

meta_full[meta_full$tissue_id == "Adipose_Subcutaneous",]$tissue_id <- "adipose-subc"
meta_full[meta_full$tissue_id == "Adipose_Visceral_Omentum",]$tissue_id <- "adipose-visc"
meta_full[meta_full$tissue_id == "Artery_Aorta",]$tissue_id <- "artery-aort"
meta_full[meta_full$tissue_id == "Artery - Aorta",]$tissue_id <- "artery-aort"
meta_full[meta_full$tissue_id == "Artery_Coronary",]$tissue_id <- "artery-coro"
meta_full[meta_full$tissue_id == "Artery - Coronary",]$tissue_id <- "artery-coro"
meta_full[meta_full$tissue_id == "Artery_Tibial",]$tissue_id <- "artery-tibi"
meta_full[meta_full$tissue_id == "Brain_Amygdala",]$tissue_id <- "brain-amyg"
meta_full[meta_full$tissue_id == "Brain_Anterior_cingulate_cortex_BA24",]$tissue_id <- "brain-ante"
meta_full[meta_full$tissue_id == "Brain_Caudate_basal_ganglia",]$tissue_id <- "brain-caud"
meta_full[meta_full$tissue_id == "Brain_Cerebellar_Hemisphere",]$tissue_id <- "brain-cehe"
meta_full[meta_full$tissue_id == "Brain_Cerebellum",]$tissue_id <- "brain-cere"
meta_full[meta_full$tissue_id == "Brain_Cortex",]$tissue_id <- "brain-cort"
meta_full[meta_full$tissue_id == "Brain_Frontal_Cortex_BA9",]$tissue_id <- "brain-frco"
meta_full[meta_full$tissue_id == "Brain_Hippocampus",]$tissue_id <- "brain-hipp"
meta_full[meta_full$tissue_id == "Brain_Hypothalamus",]$tissue_id <- "brain-hypo"
meta_full[meta_full$tissue_id == "Brain_Nucleus_accumbens_basal_ganglia",]$tissue_id <- "brain-nucl"
meta_full[meta_full$tissue_id == "Brain_Putamen_basal_ganglia",]$tissue_id <- "brain-puta"
meta_full[meta_full$tissue_id == "Brain_Spinal_cord_cervical_c-1",]$tissue_id <- "brain-spin"
meta_full[meta_full$tissue_id == "Brain_Substantia_nigra",]$tissue_id <- "brain-subs"
meta_full[meta_full$tissue_id == "Breast_Mammary_Tissue",]$tissue_id <- "breast"
meta_full[meta_full$tissue_id == "Minor_Salivary_Gland",]$tissue_id <- "salivary gland"
meta_full[meta_full$tissue_id == "Cervix_Ectocervix",]$tissue_id <- "cervix-ecto"
meta_full[meta_full$tissue_id == "Cervix_Endocervix",]$tissue_id <- "cervix-endo"
meta_full[meta_full$tissue_id == "Colon_Sigmoid",]$tissue_id <- "colon-sigm"
meta_full[meta_full$tissue_id == "Colon_Transverse",]$tissue_id <- "colon-tran"
meta_full[meta_full$tissue_id == "Colon - Transverse",]$tissue_id <- "colon-tran"
meta_full[meta_full$tissue_id == "Esophagus_Gastroesophageal_Junction",]$tissue_id <- "esophagus-gaju"
meta_full[meta_full$tissue_id == "Esophagus_Mucosa",]$tissue_id <- "esophagus-muco"
meta_full[meta_full$tissue_id == "Esophagus - Mucosa",]$tissue_id <- "esophagus-muco"
meta_full[meta_full$tissue_id == "Esophagus_Muscularis",]$tissue_id <- "esophagus-musc"
meta_full[meta_full$tissue_id == "Esophagus - Muscularis",]$tissue_id <- "esophagus-musc"
meta_full[meta_full$tissue_id == "Cells_Cultured_fibroblasts",]$tissue_id <- "fibroblasts"
meta_full[meta_full$tissue_id == "Heart_Atrial_Appendage",]$tissue_id <- "heart-atri"
meta_full[meta_full$tissue_id == "Heart_Left_Ventricle",]$tissue_id <- "heart-vent"
meta_full[meta_full$tissue_id == "Kidney_Cortex",]$tissue_id <- "kidney-cort"
meta_full[meta_full$tissue_id == "Kidney - Cortex",]$tissue_id <- "kidney-cort"
meta_full[meta_full$tissue_id == "Kidney_Medulla",]$tissue_id <- "kidney-medu"
meta_full[meta_full$tissue_id == "Cells_EBV-transformed_lymphocytes",]$tissue_id <- "lymphocytes"
meta_full[meta_full$tissue_id == "Cells - EBV-transformed lymphocytes",]$tissue_id <- "lymphocytes"
meta_full[meta_full$tissue_id == "Muscle_Skeletal",]$tissue_id <- "muscle"
meta_full[meta_full$tissue_id == "Nerve_Tibial",]$tissue_id <- "nerve"
meta_full[meta_full$tissue_id == "Skin_Not_Sun_Exposed_Suprapubic",]$tissue_id <- "skin-supr"
meta_full[meta_full$tissue_id == "Skin_Sun_Exposed_Lower_leg",]$tissue_id <- "skin-lleg"
meta_full[meta_full$tissue_id == "Skin - Sun Exposed (Lower leg)",]$tissue_id <- "skin-lleg"
meta_full[meta_full$tissue_id == "Small_Intestine_Terminal_Ileum",]$tissue_id <- "small intestine"
meta_full[meta_full$tissue_id == "Adrenal_Gland",]$tissue_id <- "adrenal gland"
meta_full[meta_full$tissue_id == "Bladder",]$tissue_id <- "bladder"
meta_full[meta_full$tissue_id == "Liver",]$tissue_id <- "liver"
meta_full[meta_full$tissue_id == "Lung",]$tissue_id <- "lung"
meta_full[meta_full$tissue_id == "Pancreas",]$tissue_id <- "pancreas"
meta_full[meta_full$tissue_id == "Pituitary",]$tissue_id <- "pituitary"
meta_full[meta_full$tissue_id == "Spleen",]$tissue_id <- "spleen"
meta_full[meta_full$tissue_id == "Stomach",]$tissue_id <- "stomach"
meta_full[meta_full$tissue_id == "Thymocytes",]$tissue_id <- "thymocytes"
meta_full[meta_full$tissue_id == "Thyroid",]$tissue_id <- "thyroid"
meta_full[meta_full$tissue_id == "Whole_Blood",]$tissue_id <- "whole blood"
meta_full[meta_full$tissue_id == "Whole Blood",]$tissue_id <- "whole blood"
meta_full[meta_full$tissue_id == "Uterus",]$tissue_id <- "uterus"
meta_full[meta_full$tissue_id == "Vagina",]$tissue_id <- "vagina"
meta_full[meta_full$tissue_id == "Ovary",]$tissue_id <- "ovary"
meta_full[meta_full$tissue_id == "Fallopian_Tube",]$tissue_id <- "fallopian tube"

meta_tissues <- subset(meta_full, meta_full$`entity:sample_id` %in% gsub(unique(df_TPM_tissues$variable), pattern = "\\.", replacement = "-"))
meta_tissues[meta_tissues$tissue_id %in% c("Cells_EBV-transformed_lymphocytes"),]$tissue_id <- "Cells_EBV_transformed_lymphocytes"

df_TPM_tissues_anno <- merge(df_TPM_tissues, meta_tissues, by.x = "variable", by.y = "entity:sample_id")

df_TPM_femmes <- df_TPM_tissues_anno[df_TPM_tissues_anno$sex == "Female",]





# TPM figure, our classifications
df_TPM_ours <- merge(df_TPM_femmes, dplyr::select(suppl_table6, gene, category), by = "gene")

df_TPM_ours[df_TPM_ours$category == "monoallelic",]$category <- "inactive"
df_TPM_ours[df_TPM_ours$category == "biallelic",]$category <- "escape"

df_TPM_ours_keep <- df_TPM_ours %>% dplyr::group_by(tissue_id, gene) %>% dplyr::summarise(meanz=mean(value))
df_TPM_ours_filt <- merge(df_TPM_ours, df_TPM_ours_keep, by = c("gene", "tissue_id"))
df_TPM_ours_filt <- df_TPM_ours[df_TPM_ours$value > 1,]

df_TPM_ours_stats2 <- df_TPM_ours_filt %>% dplyr::group_by(category, tissue_id) %>% rstatix::get_summary_stats(value, type = "common")

stat.test_our <- compare_means(ref.group = "variable", value ~ category,  data = df_TPM_ours_filt, method = "wilcox.test",p.adjust.method = "BH")
stat.test_our$y.position <- c(200,250,275)


# TPM figure, Tukiainens classifications
df_TPM_tuki <- merge(df_TPM_femmes, dplyr::select(tukiainen[tukiainen$Reported_XCI_status != "Unknown",], Gene_name, Reported_XCI_status), by.x = "gene", by.y = "Gene_name")

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



ggsave2(filename = paste0(plot_dir, "Supplementary_Fig1e_new.pdf"), width=45, height = 100, unit = "mm",
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



ggsave2(filename = paste0(plot_dir, "Fig2h_new.pdf"), width=4, height = 3,
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
