library(dplyr)
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(ggridges)
library(UpSetR)
library(ggplot2)
library(ggplotify)
library(clusterProfiler)
library(org.Hs.eg.db)
library(magick)
library(seqwrap)


# Load original splicing data 
all_splicing <- readRDS("data/all_splice.RDS") %>%
  drop_na()
all_metadata <- readRDS("data/all_full_metadata.RDS")


# Load the batch-corrected gene expression data
gene_exp_df <- readRDS("data_new/gene_counts/batch_corrected_genecounts.RDS")


# Load the informed binomial model
informed_binomial <- readRDS("data/informed_binom_model.RDS")
informed_binom_sum <- seqwrap_summarise(informed_binomial)


binom_model_outputs <- informed_binom_sum$summaries %>% 
  dplyr::select(-group) %>%
  dplyr::filter(term != "(Intercept)") %>%
  drop_na() %>%
  group_by(term) %>%                                  
  mutate(adj.p = p.adjust(p.value, method = "fdr")) %>% 
  ungroup() %>%
  filter(adj.p <= 0.05)




# Load the gene annotation file
gene_annotation <- readRDS("data/ensembl_gene_annotation.RDS")

# Load one file from which we will extract intron length
# This is valid as only introns quantified in all samples were included in the analyses
intron_length <- readr::read_tsv("data_new/Alpha_Omega_SpliceQ_outputs/A_102.tsv") %>%
  
  distinct(across(6:ncol(.)), .keep_all = T) %>% # Removes duplicates based on columns 6 to end
  mutate(transcript_ID = paste0(transcript_ID, "_", intron_ID, "_", chr),
         intron_length = abs((sj3start - sj5end) + 1) ) %>% # Ensures positive length regardless of strand
  dplyr::select(transcript_ID, intron_length)


binom_model_outputs1 <- binom_model_outputs %>%
  inner_join(intron_length, by = c("target" = "transcript_ID")) %>%
  mutate(effect = case_when(estimate > 0 & adj.p <= 0.05 ~ "Improved SE", 
                            estimate < 0 & adj.p <= 0.05 ~ "Reduced SE" ,
                            estimate < 0 & adj.p > 0.05 ~ "No effect",
                            estimate > 0 & adj.p > 0.05 ~ "No effect")) %>%
  mutate(transcript_ID = str_split(target, "_",simplify= T) [,1]) %>%
  inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version"))




Beta_binom_Model <-  readRDS("data/full_model.RDS")

beta_binom_sum <- seqwrap_summarise(Beta_binom_Model)


beta_binom_model_outputs <- beta_binom_sum$summaries %>% 
  dplyr::select(-group) %>%
  dplyr::filter(term != "(Intercept)") %>%
  drop_na() %>%
  group_by(term) %>%                                  
  mutate(adj.p = p.adjust(p.value, method = "fdr")) %>% 
  ungroup() %>%
  filter(adj.p <= 0.05)



beta_binom_model_outputs1 <- beta_binom_model_outputs %>%
  inner_join(intron_length, by = c("target" = "transcript_ID")) %>%
  mutate(effect = case_when(estimate > 0 & adj.p <= 0.05 ~ "Improved SE", 
                            estimate < 0 & adj.p <= 0.05 ~ "Reduced SE" ,
                            estimate < 0 & adj.p > 0.05 ~ "No effect",
                            estimate > 0 & adj.p > 0.05 ~ "No effect")) %>%
  mutate(transcript_ID = str_split(target, "_",simplify= T) [,1]) %>%
  inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version"))
