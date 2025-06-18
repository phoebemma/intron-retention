library(dplyr)
library(tidyverse)



source("R/Trainome_functions.R")
 # load the individual metadata 
#COPD metadata
copd_metadata <- readRDS("data_new/processed_data/copd_metadata.RDS") %>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)


# Volume_data
Vol_metadata <- readRDS("data_new/processed_data/volume_metadata.RDS")%>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age) 


# Contratrain_data
ct_metadata <- readRDS("data_new/processed_data/contratrain_metadata.RDS") %>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)


# Publicly available data
SRP102542_metadata <- readRDS("data_new/processed_data/SRP102542_metadata.RDS")%>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)


# Alpha and Omega data

A_Omega_metadata <- readRDS("data_new/processed_data/Alpha_Omega_metadata.RDS")%>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)



Relief_full_meta <- readRDS("data_new/processed_data/Relief_metadata.RDS")%>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)


# Merge the metadata into one
metadata <- rbind(copd_metadata, Vol_metadata)%>%
  rbind(ct_metadata) %>%
  rbind(A_Omega_metadata) %>%
  rbind(SRP102542_metadata) %>%
  rbind(Relief_full_meta) %>%
  mutate(across(c("age"), round, 0)) %>%
  mutate(group = case_when(age <=50 ~ "Fifty and below" ,
                           age > 50  ~ "Above fifty")) %>%
  mutate(sex = factor(sex, levels = c("female", "male")),
         group = factor(group, levels = c("Fifty and below" ,"Above fifty" )),
         time = factor(time, levels = c("PreExc", "PostExc")))
metadata$participant <- paste0(metadata$study, "_", metadata$participant)

# Standardize the age by scaling them 0 to 1
metadata$scaled_age <- round(rescale(metadata$age), digits = 2)



saveRDS(metadata, "data_new/processed_data/all_full_metadata.RDS")




# Load the full splicing data using the updated function that removes splice duplicates
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
