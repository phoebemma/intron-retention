library(trainomeMetaData)
library(dplyr)
library(tidyverse)


#Download metadata from TrainomeMetadata






 
 
 #Contratrain data
 ct_metadata <- ct_samples %>% 
   inner_join(ct_participants, by = c("study", "participant", "sex")) %>%
   dplyr::select(study, participant, sex, condition, age,  time, seq_sample_id)

 #saveRDS(ct_metadata, "data/processed_data/contratrain_metadata.RDS")
 