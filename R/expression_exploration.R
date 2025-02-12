library(dplyr)
library(tidyverse)
library(ggplot2)


low_SE <- readRDS("data_new/Pre_Exercise/annotated_low_SE_introns.RDS")

high_Se <- readRDS("data_new/Pre_Exercise/annotated_perfect_SE_introns.RDS")
length(unique(high_Se$transcript_ID))

normalized_counts <- readRDS("data_new/gene_counts/all_normalised_gene_counts.RDS")

# get the lowly SE genes in the normalised counts

low_SE_counts <- normalized_counts[normalized_counts$gene_ID %in% low_SE$ensembl_gene_id_version,]

high_SE_counts <- normalized_counts[normalized_counts$gene_ID %in% high_Se$ensembl_gene_id_version,]

copd_norm <- readRDS("data_new/gene_counts/cpm_normalised_COPD_counts.RDS")

ct_norm <- readRDS("data_new/gene_counts/cpm_normalised_contratrain_counts.RDS")
ao_norm <- readRDS("data_new/gene_counts/cpm_normalised_AO_counts.RDS")

low_se_copd <- copd_norm[copd_norm$transcript_ID %in% low_SE$transcript_ID,]


low_se_ct <- ct_norm[ct_norm$transcript_ID %in% low_SE$transcript_ID,]

low_se_ao <- ao_norm[ao_norm$transcript_ID %in% low_SE$transcript_ID,]

high_se_copd <- copd_norm[copd_norm$transcript_ID %in% high_Se$transcript_ID,]
