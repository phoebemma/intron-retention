library(dplyr)
library(tidyverse)
library(seqwrap)

# define color scale for images
## Color scale ##
colors <-  c("#d7191c", "#fdae61", "#ffffbf", "#abd9e9", "#2c7bb6", "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")

all_splice_df <- readRDS("data/Trainome_all_splice_df.RDS")



metadata <- readRDS("data/Trainome_metadata.RDS")




# Load the binary model
binom_model <- readRDS("data/Trainome_binom_model.RDS") %>%
  seqwrap_summarise()


# Load the beta-binomial model and extract its summary
beta_binom_model<- readRDS("data/Trainome_full_model.RDS") %>%
  seqwrap_summarise()



# Load the gene annotation file
gene_annotation <- readRDS("data/ensembl_gene_annotation.RDS")

# Load one file from which we will extract intron length
# This is valid as only introns quantified in all samples were included in the analyses
intron_length <- readr::read_tsv("data_new/Alpha_Omega_SpliceQ_outputs/A_102.tsv") %>%
  
  distinct(across(6:ncol(.)), .keep_all = T) %>% # Removes duplicates based on columns 6 to end
  mutate(transcript_ID = paste0(transcript_ID, "_", intron_ID, "_", chr),
         intron_length = abs((sj3start - sj5end) + 1) ) %>% # Ensures positive length regardless of strand
  dplyr::select(transcript_ID, intron_length)

#  load the gene expression dataset
gene_exp_df <- readRDS("data_new/gene_counts/batch_corrected_genecounts.RDS") 
