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

#############################################################################################
###################################### DATA PROCESSING ######################################
#############################################################################################

##### Read in data ##### 

# Annotation data
# PAR genes
PAR <- rbind(fread("annotation_data/group-715_PAR1.csv"),
             fread("annotation_data/group-716_PAR2.csv"))

# tuki df
tukiainen <- fread("annotation_data/landscape.Suppl.Table.13.csv", header = T)

# our data
# ASE tables dir
WES_file_dir_ASE_screen <- dir("supplementary_data/ASEReadCounter_output/WES", full.names = T, recursive = F)

WGS_file_dir_ASE_screen <- dir("supplementary_data/ASEReadCounter_output/WGS", full.names = T) 


# read in and add file name as column
WES_df_ase_screen <- do.call(rbind, lapply(WES_file_dir_ASE_screen, function(x) cbind(read.csv(x, header = T, sep = "\t"), name=strsplit(x,'\\.')[[1]][1])))
WGS_df_ase_screen <- do.call(rbind, lapply(WGS_file_dir_ASE_screen, function(x) cbind(read.csv(x, header = T, sep = "\t"), name=strsplit(x,'\\.')[[1]][1])))

# remove dir from file name column
WES_df_ase_screen$name <- gsub(WES_df_ase_screen$name, pattern = "supplementary_data/ASEReadCounter_output/WES/", replacement = "")
WGS_df_ase_screen$name <- gsub(WGS_df_ase_screen$name, pattern = "supplementary_data/ASEReadCounter_output/WGS/", replacement = "")

# add participant column
WES_df_ase_screen$participant <- sub(WES_df_ase_screen$name, pattern = "^([^-]*-[^-]*).*", replacement = "\\1")
WGS_df_ase_screen$participant <- sub(WGS_df_ase_screen$name, pattern = "^([^-]*-[^-]*).*", replacement = "\\1")

# WES data dir
WES_file_dir_WES_screen <- dir("supplementary_data/ASEReadCounter_output/WES_readcounts", full.names = T)
WGS_file_dir_WES_screen <- dir("supplementary_data/ASEReadCounter_output/WGS_readcounts", full.names = T) 


# read in WES data and add file name as column
WES_df_wes_screen <- do.call(rbind, lapply(WES_file_dir_WES_screen, function(x) cbind(read.csv(x, header = F, sep = "\t"), name=strsplit(x,'\\.')[[1]][1])))
WGS_df_wes_screen <- do.call(rbind, lapply(WGS_file_dir_WES_screen, function(x) cbind(read.csv(x, header = F, sep = "\t"), name=strsplit(x,'\\.')[[1]][1])))

# remove dir from file name column
WES_df_wes_screen$name <- gsub(WES_df_wes_screen$name, pattern = "supplementary_data/ASEReadCounter_output/WES_readcounts/", replacement = "")
WGS_df_wes_screen$name <- gsub(WGS_df_wes_screen$name, pattern = "supplementary_data/ASEReadCounter_output/WGS_readcounts/", replacement = "")

WES_df_ase_screen$position <- as.numeric(WES_df_ase_screen$position)
WES_df_wes_screen$V2 <- as.numeric(WES_df_wes_screen$V2)

WGS_df_ase_screen$position <- as.numeric(WGS_df_ase_screen$position)
WGS_df_wes_screen$V2 <- as.numeric(WGS_df_wes_screen$V2)



# merge WES and ASE dfs
WES_df_merged <- merge(WES_df_ase_screen, WES_df_wes_screen, by.x = c("participant", "position"), by.y = c("name", "V2"), all.x = T)
WGS_df_merged <- merge(WGS_df_ase_screen, WGS_df_wes_screen, by.x = c("participant", "position"), by.y = c("name", "V2"), all.x = T)

# change column order and rename columns
WES_df_merged <- dplyr::select(WES_df_merged, contig, position, refAllele, altAllele, V5, V6, V7, refCount, altCount, totalCount, name, participant)
WGS_df_merged <- dplyr::select(WGS_df_merged, contig, position, refAllele, altAllele, V5, V6, V7, refCount, altCount, totalCount, name, participant)


colnames(WES_df_merged) <- c("contig","position","refAllele","altAllele",
                             "altCount_WES", "refCount_WES","totalCount_WES",
                             "refCount","altCount","totalCount", 
                             "name" ,"participant")
colnames(WGS_df_merged) <- c("contig","position","refAllele","altAllele",
                             "altCount_WES", "refCount_WES","totalCount_WES",
                             "refCount","altCount","totalCount", 
                             "name" ,"participant")

WES_df_merged <- as.data.table(WES_df_merged)
WGS_df_merged <- as.data.table(WGS_df_merged)

# add gene names
chrom_order <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22", "chrX")

## load gene annotations
gtf_raw <- import("annotation_data/gencode.v41.annotation.gtf")
gtf_raw <- gtf_raw[seqnames(gtf_raw) == "chrX",]

# exclude all but standard contigs
WES_df_merged <- WES_df_merged[WES_df_merged$contig %in% chrom_order]
WGS_df_merged <- WGS_df_merged[WGS_df_merged$contig %in% chrom_order]

## add gene annotations to WES/ASE df
gr.WES <- with(WES_df_merged, GRanges(seqnames=contig, IRanges(position,width = 1)))
ol.WES <- findOverlaps(gr.WES, gtf_raw[gtf_raw$type == "gene"])
nm.WES <- tapply(gtf_raw[gtf_raw$type == "gene"]$gene_name[subjectHits(ol.WES)],queryHits(ol.WES),function(x) paste0(unique(x),collapse=";") )

WES_df_merged[, gene := NA]
WES_df_merged$gene[as.numeric(names(nm.WES))] <- nm.WES


## add gene annotations to WGS/ASE df
gr.WGS <- with(WGS_df_merged, GRanges(seqnames=contig, IRanges(position,width = 1)))
ol.WGS <- findOverlaps(gr.WGS, gtf_raw[gtf_raw$type == "exon"])
nm.WGS <- tapply(gtf_raw[gtf_raw$type == "exon"]$gene_name[subjectHits(ol.WGS)],queryHits(ol.WGS),function(x) paste0(unique(x),collapse=";") )

WGS_df_merged[, gene := NA]
WGS_df_merged$gene[as.numeric(names(nm.WGS))] <- nm.WGS

# make sure all counts are numeric
WES_df_merged$altCount <- as.numeric(WES_df_merged$altCount)
WES_df_merged$refCount <- as.numeric(WES_df_merged$refCount)
WES_df_merged$totalCount <- as.numeric(WES_df_merged$totalCount)
WES_df_merged$altCount_WES <- as.numeric(WES_df_merged$altCount_WES)
WES_df_merged$refCount_WES <- as.numeric(WES_df_merged$refCount_WES)
WES_df_merged$totalCount_WES <- as.numeric(WES_df_merged$totalCount_WES)

WGS_df_merged$altCount <- as.numeric(WGS_df_merged$altCount)
WGS_df_merged$refCount <- as.numeric(WGS_df_merged$refCount)
WGS_df_merged$totalCount <- as.numeric(WGS_df_merged$totalCount)
WGS_df_merged$altCount_WES <- as.numeric(WGS_df_merged$altCount_WES)
WGS_df_merged$refCount_WES <- as.numeric(WGS_df_merged$refCount_WES)
WGS_df_merged$totalCount_WES <- as.numeric(WGS_df_merged$totalCount_WES)


## keep only chrX
chrX_WES_df_merged <- WES_df_merged[WES_df_merged$contig == "chrX",]
chrX_WGS_df_merged <- WGS_df_merged[WGS_df_merged$contig == "chrX",]

## merge WES/ASE and WGS/ASE dfs, keep only overlapping
overlap_WGS_WES <- merge(chrX_WES_df_merged,
                         chrX_WGS_df_merged, 
                         by = c("position", "name", "gene", "contig", "refAllele", "altAllele", "participant"),
                         all=F)

# merge counts
overlap_WGS_WES$altCount <- overlap_WGS_WES$altCount.x
overlap_WGS_WES$refCount <- overlap_WGS_WES$refCount.x
overlap_WGS_WES$totalCount <- overlap_WGS_WES$totalCount.x

overlap_WGS_WES$altVariantDepth <- overlap_WGS_WES$altCount_WES.x + overlap_WGS_WES$altCount_WES.y
overlap_WGS_WES$refVariantDepth <- overlap_WGS_WES$refCount_WES.x + overlap_WGS_WES$refCount_WES.y
overlap_WGS_WES$totalVariantDepth <- overlap_WGS_WES$totalCount_WES.x + overlap_WGS_WES$totalCount_WES.x

# reorder columns
overlap_WGS_WES_fin <- dplyr::select(overlap_WGS_WES, contig, position, participant, name, gene, refAllele, altAllele, refVariantDepth, altVariantDepth, totalVariantDepth, refCount, altCount, totalCount)

# add label stating that WES and WGS overlap at that position in that sample
overlap_WGS_WES_fin$overlap <- "WGS_WES"

# find WES/ASE unique, add label that states its unique to WES
WES_unique <- chrX_WES_df_merged[!chrX_WES_df_merged$position %in% overlap_WGS_WES_fin$position,]
WES_unique <- dplyr::select(WES_unique, contig, position, participant, name, gene, refAllele, altAllele, refCount_WES, altCount_WES, totalCount_WES, refCount, altCount, totalCount)
colnames(WES_unique) <- c("contig", "position", "participant", "name", "gene", "refAllele", "altAllele","refVariantDepth", "altVariantDepth", "totalVariantDepth", "refCount", "altCount", "totalCount")
WES_unique$overlap <- "WES_only"

# find WGS/ASE unique, add label that states its unique to WGS
WGS_unique <- chrX_WGS_df_merged[!chrX_WGS_df_merged$position %in% overlap_WGS_WES_fin$position,]
WGS_unique <- dplyr::select(WGS_unique, contig, position, participant, name, gene, refAllele, altAllele, refCount_WES, altCount_WES, totalCount_WES, refCount, altCount, totalCount)
colnames(WGS_unique) <- c("contig", "position", "participant", "name", "gene", "refAllele", "altAllele","refVariantDepth", "altVariantDepth", "totalVariantDepth", "refCount", "altCount", "totalCount")
WGS_unique$overlap <- "WGS_only"

# merge the three dfs
super_smash <- rbind(overlap_WGS_WES_fin, WES_unique, WGS_unique)

# filter
super_smash$read_count_filter <- ifelse((super_smash$totalVariantDepth >= 16 & super_smash$totalCount >= 8), yes = TRUE, no = FALSE)

super_smash$alt_ref_filter <- ifelse((super_smash$altVariantDepth >= 8 & super_smash$refVariantDepth >= 8), yes = TRUE, no = FALSE) 

super_smash$minor_allele_filter <- ifelse((super_smash$refVariantDepth/super_smash$totalVariantDepth >= 0.1 & super_smash$altVariantDepth/super_smash$totalVariantDepth >= 0.1) | (super_smash$refCount/super_smash$totalCount >= 0.1 & super_smash$altCount/super_smash$totalCount >= 0.1) , yes = TRUE, no = FALSE)

super_smash$pass_all_filters <- ifelse(super_smash$minor_allele_filter == TRUE & 
                                         super_smash$read_count_filter == TRUE & 
                                         super_smash$alt_ref_filter == TRUE, 
                                       yes = TRUE, no = FALSE)

# add metadata (tissue)
sample_meta <- fread("annotation_data/sample.tsv")

# UPIC is not included in standard annotation file, use older annotation version and merge with recent
UPIC_meta <- dplyr::select(subset(fread("annotation_data/fixed_GTEx_Data_V6_Annotations_SampleAttributesDS.txt"), PARTICIPANT == "GTEX-UPIC" & SMGEBTCHT == "TrueSeq.v1"), SAMPID, SMTSD)

# change colnames to match between the annotations
colnames(UPIC_meta) <- c("entity:sample_id", "tissue_id")

# merge annotations
meta_fin <- rbind(dplyr::select(sample_meta, "entity:sample_id", tissue_id), UPIC_meta)

# shorten tissue names
meta_fin[meta_fin$tissue_id == "Adipose_Subcutaneous",]$tissue_id <- "adipose-subc"
meta_fin[meta_fin$tissue_id == "Adipose_Visceral_Omentum",]$tissue_id <- "adipose-visc"
meta_fin[meta_fin$tissue_id == "Artery_Aorta",]$tissue_id <- "artery-aort"
meta_fin[meta_fin$tissue_id == "Artery - Aorta",]$tissue_id <- "artery-aort"
meta_fin[meta_fin$tissue_id == "Artery_Coronary",]$tissue_id <- "artery-coro"
meta_fin[meta_fin$tissue_id == "Artery - Coronary",]$tissue_id <- "artery-coro"
meta_fin[meta_fin$tissue_id == "Artery_Tibial",]$tissue_id <- "artery-tibi"
meta_fin[meta_fin$tissue_id == "Brain_Amygdala",]$tissue_id <- "brain-amyg"
meta_fin[meta_fin$tissue_id == "Brain_Anterior_cingulate_cortex_BA24",]$tissue_id <- "brain-ante"
meta_fin[meta_fin$tissue_id == "Brain_Caudate_basal_ganglia",]$tissue_id <- "brain-caud"
meta_fin[meta_fin$tissue_id == "Brain_Cerebellar_Hemisphere",]$tissue_id <- "brain-cehe"
meta_fin[meta_fin$tissue_id == "Brain_Cerebellum",]$tissue_id <- "brain-cere"
meta_fin[meta_fin$tissue_id == "Brain_Cortex",]$tissue_id <- "brain-cort"
meta_fin[meta_fin$tissue_id == "Brain_Frontal_Cortex_BA9",]$tissue_id <- "brain-frco"
meta_fin[meta_fin$tissue_id == "Brain_Hippocampus",]$tissue_id <- "brain-hipp"
meta_fin[meta_fin$tissue_id == "Brain_Hypothalamus",]$tissue_id <- "brain-hypo"
meta_fin[meta_fin$tissue_id == "Brain_Nucleus_accumbens_basal_ganglia",]$tissue_id <- "brain-nucl"
meta_fin[meta_fin$tissue_id == "Brain_Putamen_basal_ganglia",]$tissue_id <- "brain-puta"
meta_fin[meta_fin$tissue_id == "Brain_Spinal_cord_cervical_c-1",]$tissue_id <- "brain-spin"
meta_fin[meta_fin$tissue_id == "Brain_Substantia_nigra",]$tissue_id <- "brain-subs"
meta_fin[meta_fin$tissue_id == "Breast_Mammary_Tissue",]$tissue_id <- "breast"
meta_fin[meta_fin$tissue_id == "Minor_Salivary_Gland",]$tissue_id <- "salivary gland"
meta_fin[meta_fin$tissue_id == "Cervix_Ectocervix",]$tissue_id <- "cervix-ecto"
meta_fin[meta_fin$tissue_id == "Cervix_Endocervix",]$tissue_id <- "cervix-endo"
meta_fin[meta_fin$tissue_id == "Colon_Sigmoid",]$tissue_id <- "colon-sigm"
meta_fin[meta_fin$tissue_id == "Colon_Transverse",]$tissue_id <- "colon-tran"
meta_fin[meta_fin$tissue_id == "Colon - Transverse",]$tissue_id <- "colon-tran"
meta_fin[meta_fin$tissue_id == "Esophagus_Gastroesophageal_Junction",]$tissue_id <- "esophagus-gaju"
meta_fin[meta_fin$tissue_id == "Esophagus_Mucosa",]$tissue_id <- "esophagus-muco"
meta_fin[meta_fin$tissue_id == "Esophagus - Mucosa",]$tissue_id <- "esophagus-muco"
meta_fin[meta_fin$tissue_id == "Esophagus_Muscularis",]$tissue_id <- "esophagus-musc"
meta_fin[meta_fin$tissue_id == "Esophagus - Muscularis",]$tissue_id <- "esophagus-musc"
meta_fin[meta_fin$tissue_id == "Cells_Cultured_fibroblasts",]$tissue_id <- "fibroblasts"
meta_fin[meta_fin$tissue_id == "Heart_Atrial_Appendage",]$tissue_id <- "heart-atri"
meta_fin[meta_fin$tissue_id == "Heart_Left_Ventricle",]$tissue_id <- "heart-vent"
meta_fin[meta_fin$tissue_id == "Kidney_Cortex",]$tissue_id <- "kidney-cort"
meta_fin[meta_fin$tissue_id == "Kidney - Cortex",]$tissue_id <- "kidney-cort"
meta_fin[meta_fin$tissue_id == "Kidney_Medulla",]$tissue_id <- "kidney-medu"
meta_fin[meta_fin$tissue_id == "Cells_EBV-transformed_lymphocytes",]$tissue_id <- "lymphocytes"
meta_fin[meta_fin$tissue_id == "Cells - EBV-transformed lymphocytes",]$tissue_id <- "lymphocytes"
meta_fin[meta_fin$tissue_id == "Muscle_Skeletal",]$tissue_id <- "muscle"
meta_fin[meta_fin$tissue_id == "Nerve_Tibial",]$tissue_id <- "nerve"
meta_fin[meta_fin$tissue_id == "Skin_Not_Sun_Exposed_Suprapubic",]$tissue_id <- "skin-supr"
meta_fin[meta_fin$tissue_id == "Skin_Sun_Exposed_Lower_leg",]$tissue_id <- "skin-lleg"
meta_fin[meta_fin$tissue_id == "Skin - Sun Exposed (Lower leg)",]$tissue_id <- "skin-lleg"
meta_fin[meta_fin$tissue_id == "Small_Intestine_Terminal_Ileum",]$tissue_id <- "small intestine"
meta_fin[meta_fin$tissue_id == "Adrenal_Gland",]$tissue_id <- "adrenal gland"
meta_fin[meta_fin$tissue_id == "Bladder",]$tissue_id <- "bladder"
meta_fin[meta_fin$tissue_id == "Liver",]$tissue_id <- "liver"
meta_fin[meta_fin$tissue_id == "Lung",]$tissue_id <- "lung"
meta_fin[meta_fin$tissue_id == "Pancreas",]$tissue_id <- "pancreas"
meta_fin[meta_fin$tissue_id == "Pituitary",]$tissue_id <- "pituitary"
meta_fin[meta_fin$tissue_id == "Spleen",]$tissue_id <- "spleen"
meta_fin[meta_fin$tissue_id == "Stomach",]$tissue_id <- "stomach"
meta_fin[meta_fin$tissue_id == "Thymocytes",]$tissue_id <- "thymocytes"
meta_fin[meta_fin$tissue_id == "Thyroid",]$tissue_id <- "thyroid"
meta_fin[meta_fin$tissue_id == "Whole_Blood",]$tissue_id <- "whole blood"
meta_fin[meta_fin$tissue_id == "Whole Blood",]$tissue_id <- "whole blood"
meta_fin[meta_fin$tissue_id == "Uterus",]$tissue_id <- "uterus"
meta_fin[meta_fin$tissue_id == "Vagina",]$tissue_id <- "vagina"
meta_fin[meta_fin$tissue_id == "Ovary",]$tissue_id <- "ovary"
meta_fin[meta_fin$tissue_id == "Fallopian_Tube",]$tissue_id <- "fallopian tube"

df_merged_meta <- merge(super_smash, meta_fin, by.x = "name", by.y = "entity:sample_id")

# calculate effectSize
df_merged_meta$effectSize <- abs(0.5 - df_merged_meta$refCount / (df_merged_meta$refCount + df_merged_meta$altCount))

# filter, keep only passing stringent filtering
df_merged_meta <- data.table(df_merged_meta[df_merged_meta$pass_all_filters == T])

# find ENSG000 genes
df_merged_meta_ensg <- df_merged_meta[grepl(df_merged_meta$gene, pattern = "ENSG0000")]

# remove rows with only ENSG0000 in the name
df_merged_meta_ensg <- df_merged_meta_ensg[grepl(df_merged_meta_ensg$gene, pattern = ";"),]
df_merged_meta_ensg$gene <- gsub(df_merged_meta_ensg$gene, pattern = "ENSG...........", replacement = "")
df_merged_meta_ensg$gene <- gsub(df_merged_meta_ensg$gene, pattern = ";", replacement = "")

# remove if empty gene name
df_merged_meta_ensg <- df_merged_meta_ensg[!df_merged_meta_ensg$gene == "",]

# remove ensg genes from original data frame
df_merged_meta <- df_merged_meta[!grepl(df_merged_meta$gene, pattern = "ENSG0000") & !grepl(df_merged_meta$gene, pattern = ";"),]

# which genes are added from the ENSG shennanigans?
setdiff(df_merged_meta_ensg$gene, df_merged_meta$gene)

# finally, merge
df_merged_meta <- rbind(df_merged_meta, df_merged_meta_ensg)

# remove overlapping genes
df_merged_meta <- df_merged_meta[!df_merged_meta$gene %in% c("CA5BP1CA5B"  ,   "GKGK-AS1"   ,    "L1CAML1CAM-AS1" ,   "SYN1TIMP1"  ,    "TSIXXIST"),]

# add keep filter to only keep hSNP per gene with highest coverage
topz <- df_merged_meta %>% dplyr::group_by(gene, tissue_id, participant) %>% top_n(1, totalCount) %>% filter(row_number()==1)
topz$keep <- T

# merge
df_merged_meta <- merge(df_merged_meta, topz, by = colnames(df_merged_meta), all = T)
df_merged_meta[is.na(df_merged_meta$keep)]$keep <- F

# shorten name
df_merged_meta$participant <- gsub(df_merged_meta$participant, pattern = "GTEX-", replacement = "")

# make unfiltered data frame
df_merged_meta_unfiltered <- df_merged_meta[!is.na(df_merged_meta$gene),]



########## Supplementary table 4 ##############

#  remove superfluous columns.
suppl_table4 <- dplyr::select(df_merged_meta_unfiltered, -name, -overlap, -read_count_filter, -alt_ref_filter, -minor_allele_filter, -pass_all_filters)

# order by hetSNP position.
suppl_table4 <- suppl_table4[order(suppl_table4$position),]

# change colnames
names(suppl_table4)[names(suppl_table4) == 'participant'] <- 'individual'
names(suppl_table4)[names(suppl_table4) == 'effectSize'] <- 'allelic_expression'

# restructure column order to make it more readable.
suppl_table4 <- dplyr::select(suppl_table4, contig, position, gene, refAllele, altAllele, refVariantDepth, altVariantDepth, totalVariantDepth, refCount, altCount, totalCount, individual, tissue_id, allelic_expression, keep)

# Write table.
write.table(suppl_table4, file = paste0(table_dir, "Supplementary_Table4.tsv"), quote = F, row.names = F, sep = "\t")




########## Supplementary table 5 ##############
# keep only highest covered hSNP.
suppl_table4_filtered <- data.table(suppl_table4[suppl_table4$keep == T,])

# Remove superfluous columns
suppl_table5 <- dplyr::select(suppl_table4_filtered, gene, individual, tissue_id, allelic_expression)

# add AE tag
suppl_table5$AE <- ifelse(suppl_table5$allelic_expression < 0.4, yes = "biallelic", no = "monoallelic")

# merge individual and tissues
suppl_table5$individual_tissueID <- paste(suppl_table5$individual, suppl_table5$tissue_id, sep = "_")

suppl_table5_counts <- suppl_table5 %>% dplyr::group_by(gene, AE) %>% count()
colnames(suppl_table5_counts) <- c("gene", "AE", "count_per_category")

dcast_suppl_table5_counts <- dcast(suppl_table5_counts, formula = gene~AE, value.var = "count_per_category")

dcast_suppl_table5_counts[is.na(dcast_suppl_table5_counts$biallelic),]$biallelic <- 0
dcast_suppl_table5_counts[is.na(dcast_suppl_table5_counts$monoallelic),]$monoallelic <- 0

dcast_suppl_table5_counts$perc_bial <- dcast_suppl_table5_counts$biallelic/(dcast_suppl_table5_counts$monoallelic+dcast_suppl_table5_counts$biallelic)*100

dcast_suppl_table5_counts$category <- "placeholder"
dcast_suppl_table5_counts[dcast_suppl_table5_counts$perc_bial >= 85,]$category <- "biallelic"
dcast_suppl_table5_counts[dcast_suppl_table5_counts$perc_bial <= 15,]$category <- "monoallelic"
dcast_suppl_table5_counts[dcast_suppl_table5_counts$perc_bial < 85 & dcast_suppl_table5_counts$perc_bial > 15,]$category <- "variable"

dcast_suppl_table5_counts[dcast_suppl_table5_counts$gene %in% PAR$`Approved symbol`,]$category <- "PAR"

# write table
write.table(dcast_suppl_table5_counts, file = paste0(table_dir, "Supplementary_Table5.tsv"), quote = F, row.names = F, sep = "\t")


########## Supplementary table 6 ##############
suppl_table6 <- dplyr::select(dcast_suppl_table5_counts, gene, category)

suppl_table6[suppl_table6$gene %in% c("ANOS1", "PRKX", "ARSD", "MED14", "NAA10", "CA5B", "SMC1A", "DIP2KB", "TMEM187", "TRAPPC2", "UBA1"),]$category <- "curated variable"
suppl_table6[suppl_table6$category == "variable",]$category <- "Undeterminable"

# write table
write.table(suppl_table6, file = paste0(table_dir, "Supplementary_Table6.tsv"), quote = F, row.names = F, sep = "\t")



########## Supplementary table 7 ##############
suppl_table7 <- unique(dplyr::select(df_merged_meta_unfiltered, name, tissue_id))

names(suppl_table7)[names(suppl_table7) == 'name'] <- 'sample_id'

# Write table.
write.table(suppl_table7, file = paste0(table_dir, "Supplementary_Table7.tsv"), quote = F, row.names = F, sep = "\t")