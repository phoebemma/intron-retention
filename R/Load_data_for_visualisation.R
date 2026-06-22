library(dplyr)
library(tidyverse)
library(seqwrap)

# define color scale for images
## Color scale ##
colors <-  c("#d7191c", "#fdae61", "#ffffbf", "#abd9e9", "#2c7bb6", "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")

all_splice_df <- readRDS("data/Trainome_all_splice_df.RDS")



metadata <- readRDS("data/Trainome_metadata.RDS")



# Load the gene annotation file
gene_annotation <- readRDS("data/ensembl_gene_annotation.RDS")

# Load one file from which we will extract intron length
# This is valid as only introns quantified in all samples were included in the analyses
intron_length <- readr::read_tsv("data_new/Alpha_Omega_SpliceQ_outputs/A_102.tsv") %>%
  distinct(across(6:ncol(.)), .keep_all = T) %>%
  mutate(transcript_ID = paste0(transcript_ID, "_", intron_ID, "_", chr),
         intron_length = abs((sj3start - sj5end) + 1)) %>%
  group_by(gene_ID) %>%
  mutate(number_introns = n()) %>%
  ungroup()
#  load the gene expression dataset
gene_exp_df <- readRDS("data_new/gene_counts/batch_corrected_genecounts.RDS") 



# ZERO-INFLATED BETABINOMIAL MODEL WITHOUT INTERACTION 

# The prediction with type = "response
zi_pred_resp <- readRDS("data/splined_pred_resp.RDS")


# The prediction with type = "zprob"

zi_pred_zprob <- readRDS("data/splined_pred_zero_inf.RDS")


# The prediction with type = "conditional

zi_pred_cond <- readRDS("data/splined_pred_cond.RDS")


# age slopes of the zero inflated model

zi_age_slopes <- readRDS("data/splined_age_slopes.RDS")




## THE ZI MODEL THAT INVESTIGATES INTERACTION BETWEEN AGING AND EXERCISE

# Only type = "response" was used in this

zi_interaction_pred <- readRDS("data/zi_interaction_predictions.RDS")

# showing the interaction between aging and exercise which is best investigated using this model 

zi_interaction_time <- readRDS( "data/zi_time_effects.RDS")
