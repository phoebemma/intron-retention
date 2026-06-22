library(reliefdata)
library(dplyr)
library(tidyverse)
library(stringi)
source("R/archive/Trainome_functions.R")


# Load the Relief spliceq data

Relief_df <- extract_splice_q_updated("data_new/Relief_SpliceQ_outputs/")


# Round the values to two decimal places
idx <- sapply(Relief_df, class)== "numeric"
Relief_df[, idx] <- lapply(Relief_df[, idx], round, 2)





Relief_metadata <- relief_seqsamples %>%
  select(-weight) %>%
  mutate(participant = as.character(participant),
         seq_sample_id = sub("^R", "R_", seq_sample_id),
         time = case_when(time == "t1" ~ "PreExc",
                          time == "t2" ~ "MidExc",
                          time == "t3" ~ "PostExc"),
         study = "ReLiEf") %>%
  inner_join(relief_participants, by = "participant") %>%
  
  inner_join(relief_volume %>%
               mutate(participant = as.character(participant)), by = c("participant", "leg")) %>%
  # filter only the pre and post exercise samples and intervention samples
  filter(time == "PreExc" | time == "PostExc"  & allocation == "int")  %>%

  dplyr::select(study, participant, sex, time, seq_sample_id, age) %>%
  filter(seq_sample_id %in% colnames(Relief_df))
 



# Select the pre-and post exercise splicing data
Relief_intersect_full <- intersect(colnames(Relief_df), Relief_metadata$seq_sample_id)

Relief_df <-Relief_df %>%
  subset(select = c("transcript_ID", Relief_intersect_full))

# Save splicing data 
saveRDS(Relief_df, "data/Relief_splicing_data.RDS")


saveRDS(Relief_metadata, "data/Relief_metadata.RDS")
