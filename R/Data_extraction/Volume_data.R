library(trainomeMetaData)
library(dplyr)
library(tidyverse)
source("R/Trainome_functions.R")
#Volume data
data(vol_participants)
data(vol_samples)

Vol_metadata <- vol_samples %>%
  inner_join(vol_participants, by = c("study", "participant", "sex"))%>%
  dplyr::select(study, participant, sex, condition, time, seq_sample_id, age )%>%
  mutate(volume = ifelse(condition == "single", 1, 3))%>%
  #pre exercise and post exercise are refered to as w0 and w12 respectively
  subset(time == "w0"| time == "w12") %>%
  drop_na()

unique(Vol_metadata$Volume)
#change the time variables to PreExc and PostExc
Vol_metadata["time"][Vol_metadata["time"] == "w0"] <- "PreExc"
Vol_metadata["time"][Vol_metadata["time"] == "w12"] <- "PostExc"
# Make the conditions all RM10 as that is what they are
Vol_metadata["condition"] <- "RM10"
unique(Vol_metadata$time)
unique(Vol_metadata$condition)
# saveRDS(Vol_metadata, "data/processed_data/volume_metadata.RDS") 



# Load Volume splicing data
volume_data <- extract_splice_q("./data/Volume_SpliceQ_outputs/")
idx <- sapply(volume_data, class)== "numeric"
volume_data[, idx] <- lapply(volume_data[, idx], round, 2)

#remove everything before the . in sample_id. 
colnames(volume_data) <- gsub(".*?\\.", "", colnames(volume_data) )


#select only the splicing samples captured in the metadata
volume_intersect <- intersect(colnames(volume_data), Vol_metadata$seq_sample_id)

volume_data <- volume_data %>%
  subset(select = c("transcript_ID", volume_intersect))

#saveRDS(volume_data, "./data/processed_data/volume_splicing_data.RDS")



#Subset Prexercise metadata
Vol_metadata_pre <- Vol_metadata %>%
  subset(time == "PreExc")


# saveRDS(Vol_metadata_pre, "data/preexercise_data/vol_preExc_metadata.RDS")

#select only the splicing samples captured in the metadata
volume_intersect <- intersect(colnames(volume_data), Vol_metadata_pre$seq_sample_id)

volume_pre_data <- volume_data %>%
  subset(select = c("transcript_ID", volume_intersect))

# saveRDS(volume_pre_data, "data/preexercise_data/vol_preExc_splicing_data.RDS")
