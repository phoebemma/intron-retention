library(trainomeMetaData)
library(dplyr)
library(tidyverse)


#Download metadata from TrainomeMetadata

#COPD metadaa

data("copd_samples")
data("copd_participants")

copd_metadata <- copd_samples %>%
  inner_join(copd_participants, by = c("study", "participant", "sex", "treatment"))%>%
  dplyr::select(study, participant, sex, condition, treatment, time, seq_sample_id, age, diagnosis )
colnames(copd_metadata)

#saveRDS(copd_metadata, "data/processed_data/copd_metadata.RDS")


#Volume data
 data(vol_participants)
 data(vol_samples)
 
 volume_metadata <- vol_samples %>%
   inner_join(vol_participants, by = c("study", "participant", "sex"))%>%
   select(study, participant, sex, condition, time, age, seq_sample_id)
#saveRDS(volume_metadata, "data/processed_data/volume_metadata.RDS") 
 
 
 #Contratrain data
 ct_metadata <- ct_samples %>% 
   inner_join(ct_participants, by = c("study", "participant", "sex")) %>%
   dplyr::select(study, participant, sex, condition, age,  time, seq_sample_id)

 #saveRDS(ct_metadata, "data/processed_data/contratrain_metadata.RDS")
 