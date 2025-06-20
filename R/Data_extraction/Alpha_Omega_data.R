library(AOData)
library(dplyr)
library(tidyverse)
library(stringi)
source("R/Trainome_functions.R")


#Load the SpliceQ data, from which we would get the sequence IDs

AO_splice <- extract_splice_q_updated("data_new/Alpha_Omega_SpliceQ_outputs/")

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
   mutate(participant = as.character(subject),
          study = "Alpha/Omega") %>%
   drop_na()






# Change the sex identification to match the other datasets
Sequenced_samples["sex"][Sequenced_samples["sex"] == "m" ] <- "male"

Sequenced_samples["sex"][Sequenced_samples["sex"] == "f" ] <- "female"
unique(Sequenced_samples$time)
unique(Sequenced_samples$sex)



length(unique(Sequenced_samples$participant))

  saveRDS(Sequenced_samples, "data/Alpha_Omega_metadata.RDS")



AO_full_intersect <- intersect(colnames(AO_splice), Sequenced_samples$seq_sample_id)

AO_splice <-AO_splice %>%
  subset(select = c("transcript_ID", AO_full_intersect)) 

# # Save the splicing data 
 saveRDS(AO_splice, "data/Alpha_Omega_splicing_data.RDS")
