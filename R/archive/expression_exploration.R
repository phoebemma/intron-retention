library(dplyr)
library(tidyverse)
library(ggplot2)


low_SE <- readRDS("data_new/Pre_Exercise/annotated_low_SE_introns.RDS")

high_Se <- readRDS("data_new/Pre_Exercise/annotated_perfect_SE_introns.RDS")
length(unique(high_Se$transcript_ID))

normalized_counts <- readRDS("data_new/gene_counts/all_normalised_genename_counts.RDS")

# get the lowly SE genes in the normalised counts

low_SE_counts <- normalized_counts[normalized_counts$transcript_ID %in% low_SE$transcript_ID,]

high_SE_counts <- normalized_counts[normalized_counts$transcript_ID %in% high_Se$transcript_ID,]


RT_affected_genes <- readRDS("data_new/processed_data/annotation_RT_effects.RDS")

RT <- normalized_counts[normalized_counts$gene_name %in% RT_affected_genes$external_gene_name,]
