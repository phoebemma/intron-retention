library(dplyr)
library(tidyverse)



source("R/Trainome_functions.R")
 # load the metadata 


metadata <- readRDS("data_new/processed_data/all_full_metadata.RDS")

# Load the different metadata


# Load the full splicing data
copd_splice_df <- extract_splice_q_updated("data_new/COPD_spliceq_outputs/")

vol_splice_df <- extract_splice_q_updated("data_new/Volume_spliceq_outputs/")
ct_splice_df <- extract_splice_q_updated("data_new/Contratrain_SpliceQ_output/")
SRP102542_splice_df <- extract_splice_q_updated("data_new/SRP102542_SpliceQ_outputs/")
AOD_splice_df <- extract_splice_q_updated("data_new/Alpha_Omega_SpliceQ_outputs/")
Relief_full_splice <- extract_splice_q_updated("data_new/Relief_SpliceQ_outputs/")


#  copd
intsect <- intersect(colnames(copd_splice_df), metadata$seq_sample_id)

# select only sample_ids in the metadata
copd_splice <- copd_splice_df %>%
  subset(select = c("transcript_ID", intsect))

saveRDS(copd_splice, "data/copd_splice.RDS")
#Volume
# make samplename match that in metadata, by removing everything before "."
colnames(vol_splice_df) <- gsub(".*?\\.", "", colnames(vol_splice_df) )

intsect <- intersect(colnames(vol_splice_df), metadata$seq_sample_id)



vol_splice <- vol_splice_df %>%
  subset(select = c("transcript_ID", intsect))

saveRDS(vol_splice, "data/volume_splice.RDS")

# Contratrain
# remove everything after the underscore
colnames(ct_splice_df)[-1] <- gsub("\\_S\\d+", "", colnames(ct_splice_df)[-1] )
intsect <- intersect(colnames(ct_splice_df), metadata$seq_sample_id)


ct_splice <- ct_splice_df %>%
  subset(select = c("transcript_ID", intsect))

saveRDS(ct_splice, "data/contratrain_splice.RDS")
# SRP102542

intsect <- intersect(colnames(SRP102542_splice_df), metadata$seq_sample_id)


SRP102542_splice <- SRP102542_splice_df %>%
  subset(select = c("transcript_ID", intsect))


saveRDS(SRP102542_splice, "data/srp102542_splice.RDS")
# Alpha_omega

intsect <- intersect(colnames(AOD_splice_df), metadata$seq_sample_id)


AOD_splice <- AOD_splice_df %>%
  subset(select = c("transcript_ID", intsect))


saveRDS(AOD_splice, "data/Alpha_Omege_splice.RDS")
# Relief
intsect <- intersect(colnames(Relief_full_splice), metadata$seq_sample_id)

Relief_splice <- Relief_full_splice %>%
  subset(select = c("transcript_ID", intsect))

saveRDS(Relief_splice, "data/relief_splice.RDS")


all_splice_df <-copd_splice  %>%
  inner_join(vol_splice, by = "transcript_ID") %>%
  inner_join(ct_splice, by = "transcript_ID") %>%
  inner_join(AOD_splice, by = "transcript_ID") %>%
  inner_join(SRP102542_splice, by = "transcript_ID") %>%
  inner_join(Relief_splice, by = "transcript_ID") 


saveRDS(all_splice_df, "data/all_splice.RDS")
