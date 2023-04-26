table_dir <- "project/Genome_biology_submission/supplementary_tables/"
plot_dir <- "project/Genome_biology_submission/plots/"


color_valzz <- c("brain-amyg" = "#0072B5FF",
                 "brain-ante" ="#0072B5FF" ,
                 "brain-caud" ="#0072B5FF",
                 "brain-cehe" = "#0072B5FF", 
                 "brain-cere" = "#0072B5FF" ,
                 "brain-cort" = "#0072B5FF" ,
                 "brain-frco" = "#0072B5FF" ,
                 "brain-hipp" = "#0072B5FF" ,
                 "brain-hypo" = "#0072B5FF",
                 "brain-nucl" = "#0072B5FF",
                 "brain-puta" = "#0072B5FF",
                 "brain-spin" = "#0072B5FF",
                 "brain-subs" ="#0072B5FF",
                 "esophagus-gaju" = "#4DBBD5FF" ,
                 "esophagus-muco" ="#4DBBD5FF",
                 "esophagus-musc" = "#4DBBD5FF",
                 "skin-lleg" ="#00A087FF",
                 "skin-supr" = "#00A087FF",
                 "heart-atri" = "#3C5488FF",
                 "heart-vent" = "#3C5488FF",
                 "adipose-subc" = "#F39B7FFF",
                 "adipose-visc" ="#F39B7FFF",
                 "artery-aort" = "#91D1C2FF",
                 "artery-coro" = "#91D1C2FF" ,
                 "artery-tibi" = "#91D1C2FF",
                 "adrenal gland" = "#8491B4FF",
                 "kidney-cort" = "#7E6148FF",
                 "kidney-medu" = "#7E6148FF",
                 "bladder" = "#DC0000FF",
                 "nerve" = "#B09C85FF",
                 "salivary gland" = "#BC3C29FF",
                 "thyroid" = "#E64B35FF",
                 "pituitary" = "#E18727FF",
                 "pancreas" = "#20854EFF",
                 "spleen" = "#7876B1FF",
                 "liver" = "#6F99ADFF",
                 "lung" = "#FFDC91FF",
                 "breast" ="#EE4C97FF",
                 "muscle" = "#878787",
                 "small intestine" = "#DF8F44FF",
                 "colon-sigm" = "#DF8F44FF",
                 "colon-tran" = "#DF8F44FF",
                 "stomach" ="#DF8F44FF",
                 "lymphocytes" = "grey",
                 "fibroblasts" = "#B24745FF",
                 "whole blood" = "#79AF97FF",
                 "cervix-ecto" = "#FFA500",
                 "cervix-endo" = "#FFA500",
                 "uterus" = "#FFA500",
                 "vagina" = "#FFA500",
                 "ovary" = "#FFA500",
                 "fallopian tube" = "#FFA500")


shape_valzz <- c("brain-amyg" = 1,
                 "brain-ante" = 2 ,
                 "brain-caud" = 3,
                 "brain-cehe" = 4, 
                 "brain-cere" = 5,
                 "brain-cort" = 6,
                 "brain-frco" = 7,
                 "brain-hipp" = 8,
                 "brain-hypo" = 9,
                 "brain-nucl" = 10,
                 "brain-puta" = 11,
                 "brain-spin" = 12,
                 "brain-subs" = 13,
                 "esophagus-gaju" = 1,
                 "esophagus-muco" = 2,
                 "esophagus-musc" = 3,
                 "skin-lleg" = 1,
                 "skin-supr" = 2,
                 "heart-atri" = 1,
                 "heart-vent" = 2,
                 "adipose-subc" = 1,
                 "adipose-visc" = 2,
                 "artery-aort" = 1,
                 "artery-coro" = 2 ,
                 "artery-tibi" = 3,
                 "kidney-cort" = 1,
                 "kidney-medu" = 2,
                 "colon-sigm" = 1,
                 "colon-tran" = 2,
                 "adrenal gland" = 1,
                 "bladder" = 1,
                 "nerve" = 1,
                 "salivary gland" = 1,
                 "thyroid" = 1,
                 "pituitary" = 1,
                 "pancreas" = 1,
                 "spleen" = 1,
                 "liver" = 1,
                 "lung" = 1,
                 "breast" = 1,
                 "muscle" = 1,
                 "small intestine" = 1,
                 "stomach" = 1,
                 "lymphocytes" = 1,
                 "fibroblasts" = 1,
                 "whole blood" = 1,
                 "cervix-ecto" = 1,
                 "cervix-endo" = 2,
                 "uterus" = 1,
                 "vagina" = 1,
                 "ovary" = 1,
                 "fallopian tube" = 1)

new_tissue_cats_color_valzz <- c(
  "brain" = "#0072B5FF",
  "esophagus" = "#4DBBD5FF" ,
  "skin" ="#00A087FF",
  "heart" = "#3C5488FF",
  "adipose" = "#F39B7FFF",
  "artery" = "#91D1C2FF",
  "adrenal gland" = "#8491B4FF",
  "kidney" = "#7E6148FF",
  "bladder" = "#DC0000FF",
  "nerve" = "#B09C85FF",
  "salivary gland" = "#BC3C29FF",
  "thyroid" = "#E64B35FF",
  "pituitary" = "#E18727FF",
  "pancreas" = "#20854EFF",
  "spleen" = "#7876B1FF",
  "liver" = "#6F99ADFF",
  "lung" = "#FFDC91FF",
  "breast" ="#EE4C97FF",
  "muscle" = "#878787",
  "small intestine" = "#DF8F44FF",
  "colon" = "#DF8F44FF",
  "stomach" ="#DF8F44FF",
  "lymphocytes" = "grey",
  "fibroblasts" = "#B24745FF",
  "whole blood" = "#79AF97FF",
  "cervix" = "#FFA500",
  "uterus" = "#FFA500",
  "vagina" = "#FFA500",
  "ovary" = "#FFA500",
  "fallopian tube" = "#FFA500")


## Functions to make plots look similar across plottings. Written by Antonio Lentini.
require(ggplot2)
theme_AL_simple <- function (...) { 
  theme_bw(base_size=12, base_family="") %+replace% 
    theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(colour="black"),axis.line=element_line(colour="black"), ... )
}
theme_AL_simple_rotX <- function (...) { 
  theme_bw(base_size=12, base_family="") %+replace% 
    theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line=element_line(colour="black"),axis.text = element_text(colour="black"),axis.text.x=element_text(angle = 90,hjust=1), ... )
}
theme_AL_box <- function (...) { 
  theme_bw(base_size=12, base_family="") %+replace% 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(colour="black"),axis.line=element_line(colour="black"), ... )
}
theme_AL_box_rotX <- function (...) { 
  theme_bw(base_size=12, base_family="") %+replace% 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(colour="black"),axis.line=element_line(colour="black"), axis.text.x=element_text(angle = 90,hjust=1), ... )
}
theme_AL_box_rotX_45 <- function (...) { 
  theme_bw(base_size=12, base_family="") %+replace% 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(colour="black"),axis.line=element_line(colour="black"), axis.text.x=element_text(angle = 45), ... )
}
theme_AL_simple_rotX_45 <- function (...) { 
  theme_bw(base_size=12, base_family="") %+replace% 
    theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line=element_line(colour="black"),axis.text = element_text(colour="black"),axis.text.x=element_text(angle = 45), ... )
}