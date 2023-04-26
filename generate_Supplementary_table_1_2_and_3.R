library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(rtracklayer)
library(GenomicRanges)

# plot parameters
source("project/Genome_biology_submission/r_sources/plot_parameters.R")

# RNA-seq samples selected
selected_samples <- fread("project/Genome_biology_submission/supplementary_tables/selected_samples.tsv")

# ASE tables dir
file_dir_ASE_screen <- dir("project/Genome_biology_submission/supplementary_data/screen/ase", full.names = T) 

# read in and add file name as column
df_ase_screen <- do.call(rbind, lapply(file_dir_ASE_screen, function(x) cbind(read.csv(x, header = T, sep = "\t"), name=strsplit(x,'\\.')[[1]][1])))

# remove dir from file name column
df_ase_screen$name <- gsub(df_ase_screen$name, pattern = "project/Genome_biology_submission/supplementary_data/screen/ase/", replacement = "")

# add participant column
df_ase_screen$participant <- sub(df_ase_screen$name, pattern = "^([^-]*-[^-]*).*", replacement = "\\1")

# WES data dir
file_dir_WES_screen <- dir("project/Genome_biology_submission/supplementary_data/screen/wes", full.names = T) 

# read in WES data and add file name as column
df_wes_screen <- do.call(rbind, lapply(file_dir_WES_screen, function(x) cbind(read.csv(x, header = F, sep = "\t"), name=strsplit(x,'\\.')[[1]][1])))

# remove dir from file name column
df_wes_screen$name <- gsub(df_wes_screen$name, pattern = "project/Genome_biology_submission/supplementary_data/screen/wes/", replacement = "")

# make sure all numerical columns are numerical.
df_ase_screen$position <- as.numeric(df_ase_screen$position)
df_wes_screen$V2 <- as.numeric(df_wes_screen$V2)

# merge WES and ASE dfs
df_merged <- merge(df_ase_screen, df_wes_screen, by.x = c("participant", "position"), by.y = c("name", "V2"), all.x = T)

# change column order and rename columns
df_merged <- dplyr::select(df_merged, contig, position, refAllele, altAllele, V5, V6, V7, refCount, altCount, totalCount, name, participant)
colnames(df_merged) <- c("contig","position","refAllele","altAllele",
                         "altCount_WES", "refCount_WES","totalCount_WES",
                         "refCount","altCount","totalCount", 
                         "name" ,"participant")

# make into data table
df_merged <- as.data.table(df_merged)

# add gene names.
chrom_order <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22", "chrX")

## load gene annotations.
gtf_raw <- import("project/Genome_biology_submission/annotation_data/gencode.v41.annotation.gtf")

# exclude all but standard contigs.
df_merged <- df_merged[df_merged$contig %in% chrom_order]

## add gene annotations to merged WES/ASE df.
gr.ase <- with(df_merged, GRanges(seqnames=contig, IRanges(position,width = 1)))
ol.ase <- findOverlaps(gr.ase, gtf_raw[gtf_raw$type == "gene"])
nm.ase <- tapply(gtf_raw[gtf_raw$type == "gene"]$gene_name[subjectHits(ol.ase)],queryHits(ol.ase),function(x) paste0(unique(x),collapse=";") )
df_merged[, gene := NA]
df_merged$gene[as.numeric(names(nm.ase))] <- nm.ase

# make sure all numerical columns are numerical.
df_merged$altCount <- as.numeric(df_merged$altCount)
df_merged$refCount <- as.numeric(df_merged$refCount)
df_merged$totalCount <- as.numeric(df_merged$totalCount)

df_merged$altCount_WES <- as.numeric(df_merged$altCount_WES)
df_merged$refCount_WES <- as.numeric(df_merged$refCount_WES)
df_merged$totalCount_WES <- as.numeric(df_merged$totalCount_WES)


# Add a raw read count filter
df_merged$read_count_filter <- ifelse((df_merged$totalCount_WES >= 40 & df_merged$totalCount > 10), yes = TRUE, no = FALSE)

# Add raw alt/ref read count filter
df_merged$alt_ref_filter <- ifelse((df_merged$altCount_WES >= 20 & df_merged$refCount_WES >= 20), yes = TRUE, no = FALSE) 

# Add a $minor_allele_filter tag, where both and minor allele reads need to be equal to or more than 10 % of total reads, either in the WES or RNA-seq data. Essentially we make sure we observe both alleles somewhere in the data.
df_merged$minor_allele_filter <- ifelse((df_merged$refCount_WES/df_merged$totalCount_WES >= 0.1 & df_merged$altCount_WES/df_merged$totalCount_WES >= 0.1) | (df_merged$refCount/df_merged$totalCount >= 0.1 & df_merged$altCount/df_merged$totalCount >= 0.1) , yes = TRUE, no = FALSE)

# Add pass all filter tag
df_merged$pass_all_filters <- ifelse(df_merged$minor_allele_filter == TRUE & df_merged$read_count_filter == TRUE & df_merged$alt_ref_filter == TRUE, yes = TRUE, no = FALSE)

# remove not passing
df_merged <- df_merged[df_merged$pass_all_filters == T]

# calculate effectSize
df_merged$effectSize <- abs(0.5 - df_merged$refCount / (df_merged$refCount + df_merged$altCount))

# select inactive genes
tukiainen <- fread("project/Genome_biology_submission/annotation_data/landscape.Suppl.Table.13.csv", header = T)

# keep only inacs
df_merged_inacs <- df_merged[df_merged$gene %in% tukiainen[tukiainen$Reported_XCI_status == "Inactive",]$Gene_name,]

# add keep filter to only keep hSNP per gene with highest coverage
df_merged_inacs[, keep := order(totalCount,decreasing = T) == 1, by=c("name","gene")]

# add metadata (tissue)
sample_meta <- fread("project/Genome_biology_submission/annotation_data/sample.tsv")

# add UPIC meta (it is not included in the above metadata)
UPIC_meta <- dplyr::select(subset(fread("project/Genome_biology_submission/annotation_data/fixed_GTEx_Data_V6_Annotations_SampleAttributesDS.txt"), PARTICIPANT == "GTEX-UPIC" & SMGEBTCHT == "TrueSeq.v1"), SAMPID, SMTSD)

#change colnames to match
colnames(UPIC_meta) <- c("entity:sample_id", "tissue_id")

# bind together
meta_fin <- rbind(dplyr::select(sample_meta, "entity:sample_id", tissue_id), UPIC_meta)

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

# merge by sample name
df_merged_inacs_meta <- merge(df_merged_inacs, meta_fin, by.x = "name", by.y = "entity:sample_id")

# make the df into a data table
df_merged_inacs_meta <- data.table(df_merged_inacs_meta)




#### Make Supplementary Table 1 ####
suppl_table1 <- unique(dplyr::select(df_merged_inacs_meta[keep == T & name %in% selected_samples$sample_id],gene))

# change colnames
colnames(suppl_table1) <- c("gene")

# write table
write.table(suppl_table1, file = paste0(table_dir, "Supplementary_Table1.tsv"), quote = F, row.names = F, sep = "\t")




#### Make Supplementary Table 2 ####
# select genes which we have coverage of, and the selected tissues. Remove superfluous columns.
suppl_table2 <- dplyr::select(df_merged_inacs_meta[keep == T & name %in% selected_samples$sample_id], -name, -read_count_filter, -alt_ref_filter, -minor_allele_filter, -pass_all_filters, -keep)
# order by hetSNP position.
suppl_table2 <- suppl_table2[order(suppl_table2$position),]

# change colnames.
colnames(suppl_table2) <- c("contig","position","refAllele","altAllele","refVariantDepth","altVariantDepth","totalVariantDepth", "refCount", "altCount", "totalCount","individual","gene","allelic_expression","tissue_id")

# restructure column order to make it more readable.
suppl_table2 <- dplyr::select(suppl_table2, contig, position, gene, refAllele, altAllele, refVariantDepth, altVariantDepth, totalVariantDepth, refCount, altCount, totalCount, individual, tissue_id, allelic_expression)

# Shorten individual names (remove 'GTEX-').
suppl_table2$individual <- gsub(suppl_table2$individual, pattern = "GTEX-", replacement = "")

# Write table.
write.table(suppl_table2, file = paste0(table_dir, "Supplementary_Table2.tsv"), quote = F, row.names = F, sep = "\t")





#### Make Supplementary Table 3 ####
#Select the individuals that were interesting in Figure 1B, i.e. the ones with a mean allelic of inactive genes > 0.4.
df_merged_inacs_meta_statz <- df_merged_inacs_meta[keep == T & name %in% selected_samples$sample_id,] %>% dplyr::group_by(participant, tissue_id) %>% rstatix::get_summary_stats(effectSize, type = "common")

# select females of interest.
fois <- df_merged_inacs_meta_statz[df_merged_inacs_meta_statz$mean > 0.4,]$participant

# include 6 controls.
set.seed(190524)
selected <- c(sample(df_merged_inacs_meta_statz[!df_merged_inacs_meta_statz$participant %in% fois,]$participant, 6), fois)

# make new df with only selected samples.
df_merged_inacs_meta_selected_participants <- df_merged_inacs_meta[df_merged_inacs_meta$participant %in% selected,]

# keep only filtered and the hetSNP with the highest RNA-seq coverage.
df_merged_inacs_meta_selected_participants_fin <- df_merged_inacs_meta_selected_participants[keep == T,]

#  remove superfluous columns.
suppl_table3 <- dplyr::select(df_merged_inacs_meta_selected_participants_fin, -name, -read_count_filter, -alt_ref_filter, -minor_allele_filter, -pass_all_filters, -keep)

# order by hetSNP position.
suppl_table3 <- suppl_table3[order(suppl_table3$position),]

# change colnames
colnames(suppl_table3) <- c("contig","position","refAllele","altAllele","refVariantDepth","altVariantDepth","totalVariantDepth", "refCount", "altCount", "totalCount","individual","gene","allelic_expression","tissue_id")

# restructure column order to make it more readable.
suppl_table3 <- dplyr::select(suppl_table3, contig, position, gene, refAllele, altAllele, refVariantDepth, altVariantDepth, totalVariantDepth, refCount, altCount, totalCount, individual, tissue_id, allelic_expression)

# Shorten individual names (remove 'GTEX-').
suppl_table3$individual <- gsub(suppl_table3$individual, pattern = "GTEX-", replacement = "")

# Write table.
write.table(suppl_table3, file = paste0(table_dir, "Supplementary_Table3.tsv"), quote = F, row.names = F, sep = "\t")
