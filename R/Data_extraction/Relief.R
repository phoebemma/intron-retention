library(reliefdata)
library(dplyr)
library(tidyverse)
library(stringi)
source("R/Trainome_functions.R")


# Load the Relief spliceq data

Relief_df <- extract_splice_q_updated("data_new/Relief_SpliceQ_outputs/")


# Round the values to two decimal places
idx <- sapply(Relief_df, class)== "numeric"
Relief_df[, idx] <- lapply(Relief_df[, idx], round, 2)



# extraction sequence is that number following "R" on the sample name
seq_df <- as.data.frame(colnames(Relief_df[, -1])) %>%
  mutate(seq_id = as.double(stri_extract_first_regex(colnames(Relief_df[, -1]), "\\d+")))

# rename the first column to seq_sample_id to match the other datasets
colnames(seq_df)[colnames(seq_df) == "colnames(Relief_df[, -1])"] <- "seq_sample_id"

# Get sequence information from the excel file obtained from Kristian
Relief_metadata <- readxl::read_excel("Relief_sampleIDs.xlsx") %>%
  mutate(subject = as.character(subject),
         study = "ReLiEf")


# Change the "sample" column name to seq_id
colnames(Relief_metadata)[colnames(Relief_metadata) == "sample"] <- "seq_id"
# Rename subject to participant
colnames(Relief_metadata)[colnames(Relief_metadata) == "subject"] <- "participant"




Relief_metadata <- Relief_metadata %>%
  inner_join(seq_df, by = "seq_id") %>%
  inner_join(relief_participants, by =  "participant") %>%
  #dplyr::select(study, participant, time_rep, age, sex, seq_sample_id, allocation) %>%
  #Extract leg column
  mutate( 
          time = case_when(time_rep == "t1rnaL" ~ "PreExc",
                           time_rep == "t2rnaL" ~ "MidExc",
                           time_rep == "t3rnaL" ~ "PostExc",
                           time_rep == "t1rnaR"  ~ "PreExc",
                           time_rep == "t2rnaR"  ~ "MidExc",
                           time_rep == "t3rnaR"  ~ "PostExc")) %>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age, allocation) %>%
  # filter only the pre and post exercise samples
  filter(time == "PreExc" | time == "PostExc")




#  Remove the Relief samples that are controls
Relief_metadata <- Relief_metadata %>%
  filter(allocation == "int") # filters for the intervention participants







# Select the pre-and post exercise splicing data
Relief_intersect_full <- intersect(colnames(Relief_df), Relief_metadata$seq_sample_id)

Relief_df <-Relief_df %>%
  subset(select = c("transcript_ID", Relief_intersect_full))

# Save splicing data 
saveRDS(Relief_df, "data/Relief_splicing_data.RDS")


saveRDS(Relief_metadata, "data/Relief_metadata.RDS")
