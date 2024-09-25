library(trainomeMetaData)
library(dplyr)
library(tidyverse)
source("R/Trainome_functions.R")

# Download metadata from  the TrainomeMetadata package
data(ct_participants)
data(ct_samples)  


ct_metadata <- ct_samples %>% inner_join(ct_participants, by = c("study", "participant", "sex")) %>%
  dplyr::select(study, participant, sex, condition, time, seq_sample_id, age) %>%
  # drop the rows containing missing data
  drop_na() # %>%

# Rename the time column to match those in the other dataframes
ct_metadata["time"][ct_metadata["time"] == "t1"] <- "PreExc"
#ct_strength["time"][ct_strength["time"] == "t3"] <- "t3"
ct_metadata["time"][ct_metadata["time"] == "t4"] <- "PostExc"

ct_metadata <- ct_metadata %>%
  # Take only the Pre and post exercise data
  subset( time == "PreExc" | time == "PostExc") %>%
  # remove the untrained group as they are not of interest
  subset(condition == "set3" | condition == "set6")
unique(ct_metadata$time)
unique(ct_metadata$condition)
length(unique(ct_metadata$participant))

# saveRDS(ct_metadata, "data/processed_data/contratrain_metadata.RDS")



#Load Volume splicing data
ct_data <- extract_splice_q("./data/Contratrain_SpliceQ_outputs/")
idx <- sapply(ct_data, class)== "numeric"
ct_data[, idx] <- lapply(ct_data[, idx], round, 2)


# select only the splicing details of samples captured in the metadata

ct_intersect <- intersect(colnames(ct_data), ct_metadata$seq_sample_id)

ct_data <- ct_data %>%
  subset(select = c("transcript_ID", ct_intersect))

# saveRDS(ct_data, "data/processed_data/contratrain_splicing_data.RDS")




# subset the prexercise data
ct_pre_meta <- ct_metadata %>%
  subset(time == "PreExc")

ct_intersect <- intersect(colnames(ct_data), ct_pre_meta$seq_sample_id)

ct_pre_data <- ct_data %>%
  subset(select = c("transcript_ID", ct_intersect))
# saveRDS(ct_pre_meta, "data/preexercise_data/ct_PreExc_metadata.RDS")
# saveRDS(ct_pre_data, "data/preexercise_data/ct_PreExc_splicing_data.RDS")
