library(trainomeMetaData)
library(dplyr)
library(tidyverse)
source("R/Trainome_functions.R")

# Download and extract metadata from  the TrainomeMetadata package
ct_metadata <- ct_samples %>% inner_join(ct_participants, by = c("study", "participant", "sex")) %>%
  dplyr::select(study, participant, sex, condition, time, seq_sample_id, age) %>%
  # drop the rows containing missing data
  drop_na() %>%
  # mutate to capture volume of exercise
  mutate(volume = case_when(condition == "set0"~ 0,
                            condition == "set3" ~ 3,
                            condition == "set6" ~ 6))


# Rename the time column to match those in the other dataframes
ct_metadata["time"][ct_metadata["time"] == "t1"] <- "PreExc"
#ct_strength["time"][ct_strength["time"] == "t3"] <- "t3"
ct_metadata["time"][ct_metadata["time"] == "t4"] <- "PostExc"

# Rename the condition column to match those in the other dataframes
ct_metadata["condition"][ct_metadata["condition"] == "set3"] <- "RM10"
ct_metadata["condition"][ct_metadata["condition"] == "set0"] <- "Control"
ct_metadata["condition"][ct_metadata["condition"] == "set6"] <- "RM10"




ct_metadata <- ct_metadata %>%
  # Take only the Pre and post exercise data
  subset( time == "PreExc" | time == "PostExc")  %>%
  # remove the untrained group as they are not of interest
   filter(condition != "Control" )



 saveRDS(ct_metadata, "data/contratrain_metadata.RDS")



#Load Volume splicing data
ct_data <- extract_splice_q_updated("./data_new/Contratrain_SpliceQ_output/") 
# remove everything after the underscore
colnames(ct_data)[-1] <- gsub("\\_S\\d+", "", colnames(ct_data)[-1] )
idx <- sapply(ct_data, class)== "numeric"
ct_data[, idx] <- lapply(ct_data[, idx], round, 2)


# select only the splicing details of samples captured in the metadata

ct_intersect_full <- intersect(colnames(ct_data), ct_metadata$seq_sample_id)

ct_data_full <- ct_data %>%
  subset(select = c("transcript_ID", ct_intersect_full)) 

 saveRDS(ct_data_full, "data/contratrain_splicing_data.RDS")





