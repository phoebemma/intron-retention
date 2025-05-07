library(AOData)
library(dplyr)
library(tidyverse)
library(stringi)
source("R/Trainome_functions.R")


#Load the SpliceQ data, from which we would get the sequence IDs

AO_splice <- extract_splice_q("data_new/Alpha_Omega_SpliceQ_outputs/")

# Round the values to two decimal places
idx <- sapply(AO_splice, class)== "numeric"
AO_splice[, idx] <- lapply(AO_splice[, idx], round, 2)




# Extract the sequence ID to match the extraction sequence given in the metadata
# extraction sequence is that number following "s" on the sample name
seq_df <- as.data.frame(colnames(AO_splice[, -1])) %>%
  mutate(seq_id = as.double(stri_extract_first_regex(colnames(AO_splice[, -1]), "\\d+")))

# rename the first column to seq_sample_id to match the other datasets
colnames(seq_df)[colnames(seq_df) == "colnames(AO_splice[, -1])"] <- "seq_sample_id"

 A_O_seq_list <- readxl::read_excel("Alpha_Omega_sample_list_transcriptomics.xlsx") %>%
   filter(Tissue == "muscle") %>%
   mutate(time = case_when(time_rep == "T1rna1" ~ "PreExc",
          time_rep == "T4rna1" ~ "PostExc")) %>%
   dplyr::select(subject, Tissue, sample, time)


 Sequenced_samples <- seq_df %>%
   inner_join(A_O_seq_list, by = c("seq_id" = "sample"))%>%
   inner_join(idkeys, by = c("subject" = "participant")) %>%
   dplyr::select(seq_sample_id, seq_id, subject, time,  age, sex)%>%
   drop_na()




# Load the participant details

# 
# Sequenced_samples <- seq_samples %>%
#   dplyr::select(participant, time, condition, leg, extraction_seq) %>%
#   inner_join(ids, by = "participant")  %>%
#   inner_join(seq_df, by = c("extraction_seq" = "seq_id"))%>%
#    mutate(time = case_when(time == "T1" ~ "PreExc",
#                            time == "T2" ~ "PreTrain",
#                            time == "T4" ~ "PostExc"))
Sequenced_samples["participant"] <- as.character(Sequenced_samples$subject)



Sequenced_samples$study <- "Alpha/Omega"

# Change the sex identification to match the other datasets
Sequenced_samples["sex"][Sequenced_samples["sex"] == "m" ] <- "male"

Sequenced_samples["sex"][Sequenced_samples["sex"] == "f" ] <- "female"
unique(Sequenced_samples$time)
unique(Sequenced_samples$sex)

# select the pre and post exercise data
# Sequenced_samples <- Sequenced_samples %>%
#   filter(time == "PreExc" | time == "PostExc")

length(unique(Sequenced_samples$participant))

  saveRDS(Sequenced_samples, "data_new/processed_data/Alpha_Omega_metadata.RDS")

# Extract pre_exercise data  
  pre_Exc <- Sequenced_samples %>%
  filter(time == "PreExc") %>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)
 
length(unique(pre_Exc$participant))

saveRDS(pre_Exc, "data_new/Pre_Exercise/Alpha_Omega_PreExc_metadata.RDS")


# Select the pre-exercise splicing data
AO_intersect <- intersect(colnames(AO_splice), pre_Exc$seq_sample_id)

AO_pre_splicing <-AO_splice %>%
  subset(select = c("transcript_ID", AO_intersect)) %>%
  drop_na()

AO_full_intersect <- intersect(colnames(AO_splice), Sequenced_samples$seq_sample_id)

AO_splice <-AO_splice %>%
  subset(select = c("transcript_ID", AO_full_intersect)) %>%
  drop_na()

 saveRDS(AO_pre_splicing, "data_new/Pre_Exercise/Alpha_Omega_PreExc_splicing_data.RDS")
# # Save the splicing data 
 saveRDS(AO_splice, "data_new/processed_data/Alpha_Omega_splicing_data.RDS")
