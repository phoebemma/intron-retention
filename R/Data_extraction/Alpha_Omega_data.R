library(AOData)
library(dplyr)
library(tidyverse)
library(stringi)
source("R/Trainome_functions.R")


#Load the SpliceQ data, from which we would get the sequence IDs

AO_splice <- extract_splice_q("data/Alpha_Omega_SpliceQ_outputs/")

# Round the values to two decimal places
idx <- sapply(AO_splice, class)== "numeric"
AO_splice[, idx] <- lapply(AO_splice[, idx], round, 2)

# Save the splicing data 
saveRDS(AO_splice, "data_new/processed_data/Alpha_Omega_splicing_data.RDS")

# Extract the sequence ID to match he extraction sequence given in he metadata
# extraction sequence is that number following s on the sample name
seq_df <- as.data.frame(colnames(AO_splice[, -1])) %>%
  mutate(seq_id = as.double(stri_extract_first_regex(colnames(AO_splice[, -1]), "\\d+")))

# rename the first column to seq_sample_id to match the other datasets
colnames(seq_df)[colnames(seq_df) == "colnames(AO_splice[, -1])"] <- "seq_sample_id"



# Load the participant details
ids <- idkeys %>%
  select(participant, treat, age,  sex )

Sequenced_samples <- seq_samples %>%
  select(participant, time, condition, leg, extraction_seq) %>%
  inner_join(ids, by = "participant")  %>%
  inner_join(seq_df, by = c("extraction_seq" = "seq_id"))%>%
   mutate(time = case_when(time == "T1" ~ "PreExc",
                           time == "T2" ~ "PreTrain",
                           time == "T4" ~ "PostExc"))
Sequenced_samples["participant"] <- as.character(Sequenced_samples$participant)
# y <- c[c$extraction_seq == 920, ]
# q <- Sequenced_samples[Sequenced_samples$participant == 517, ]

Sequenced_samples$study <- "Alpha/Omega"

# Change the sex identity to match the other datasets
Sequenced_samples["sex"][Sequenced_samples["sex"] == "m" ] <- "male"

Sequenced_samples["sex"][Sequenced_samples["sex"] == "f" ] <- "female"
unique(Sequenced_samples$time)
unique(Sequenced_samples$sex)

# saveRDS(Sequenced_samples, "data_new/processed_data/Alpha_Omega_metadata.RDS")
pre_Exc <- Sequenced_samples %>%
  filter(time == "PreExc") %>%
  select(study, participant, sex, time, seq_sample_id, age)
 
# x <- seq_df %>%
#   filter(!seq_sample_id %in% Sequenced_samples$seq_sample_id)
hist(pre_Exc$age)


# saveRDS(pre_Exc, "data_new/Pre_Exercise/Alpha_Omega_PreExc_metadata.RDS")


# Select the pre-exercise metadata
AO_intersect <- intersect(colnames(AO_splice), pre_Exc$seq_sample_id)

AO_pre_splicing <-AO_splice %>%
  subset(select = c("transcript_ID", AO_intersect))


saveRDS(AO_pre_splicing, "data_new/Pre_Exercise/Alpha_Omega_PreExc_splicing_data.RDS")

 #dataframes to check out
strength

thickness

c <- tissue_samples


# # The code below runs only because I have run the preexercise_model.R in the global environment
# 
# all_pre_metadata_with_alpha <- rbind(copd_metadata, volume_metadata)%>%
#   rbind(Contratrain_metadata) %>%
#   
#   mutate(age_group = case_when(study == "copd" ~ "Old", 
#                                study == "vol" ~ "Young",
#                                study == "ct" ~ "Young")) %>%
#   rbind(SRP102542_metadata) %>%
#   rbind(SRP280348_metadata) %>% 
#   rbind(pre) %>%
#   
#   #Copd and volume and SRP102542 have age data in decimal
#   mutate(across(c("age"), round, 0)) %>%
#   # Create youps where group 1 = those below 30
#   #group 2 is those above 30 but below 51
#   # group 3 those above 50 but below 71
#   # group 4 is those above 70
#   
#   mutate(group = case_when(age <=30 ~ "<=30" ,
#                            age > 30 & age <= 40 ~ ">30 & <=40", 
#                            age > 40 & age <= 50 ~ ">40 & <=50",
#                            age > 50 & age <= 60 ~ ">50 & <=60",
#                            age > 60 & age <= 70 ~ ">60 & <=70",
#                            age > 70 ~ ">70")) %>%
#   mutate(sex = factor(sex, levels = c("female", "male")),
#          age_group = factor(age_group, levels= c("Young", "Old")),
#          group = factor(group, levels = c("<=30" , ">30 & <=40", ">40 & <=50",
#                                           ">50 & <=60",">60 & <=70", ">70"))) 
# 
# 

 