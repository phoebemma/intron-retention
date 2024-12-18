#This script shows and details the extraction of the COPD data
# Full details about the data is published in https://pubmed.ncbi.nlm.nih.gov/34229714/


library(trainomeMetaData)
library(dplyr)
library(tidyverse)
source("R/Trainome_functions.R")


#COPD metadata
data(copd_samples)
data(copd_participants)



copd_metadata <- copd_samples %>%
  inner_join(copd_participants, by = c("study", "participant", "sex", "treatment"))%>%
  dplyr::select(study, participant, sex, condition, time, seq_sample_id, age)%>%
  #select only thesubset we are interested in
  subset(time == "PreExc" | time == "PostExc")%>%
  drop_na() 


copd_metadata["volume"] <- 3
hist(copd_metadata$age)
length(unique(copd_metadata$participant))
 saveRDS(copd_metadata, "data_new/processed_data/copd_metadata.RDS")

#Load the COPD splicing data

copd_data <- extract_splice_q("./data_new/COPD_Spliceq_outputs_new/")
#Rename the ccolumn names by removing everything after the second underscore
# colnames(copd_data) <- gsub("_.*", "", colnames(copd_data) )
# #COPD sequence data contains duplicates. The code below should remove duplicate columns
# copd_data <-   copd_data[, !duplicated(colnames(copd_data))]

#the gsub functio above also removed the underscore in "transcript_ID, replace it
# colnames(copd_data)[1] <- "transcript_ID"
# #remove everything before the . in sample_id. 
# colnames(copd_data) <- gsub(".*?\\.", "", colnames(copd_data) )
# #Add an X before the sample_ID. This makes it match to the naming in metadata
# colnames(copd_data)[-1] <- paste0("X", colnames(copd_data)[-1] )

idx <- sapply(copd_data, class)== "numeric"
copd_data[, idx] <- lapply(copd_data[, idx], round, 2)



# select only the splicing samples captured in the metadata
copd_intersect <- intersect(colnames(copd_data), copd_metadata$seq_sample_id)

copd_data <- copd_data %>%
  subset(select = c("transcript_ID", copd_intersect))


 saveRDS(copd_data, "./data_new/processed_data/copd_splicing_data.RDS")



# Extract the pre-exercise data

copd_metadata_pre <- copd_metadata %>%
  subset(time == "PreExc")

 saveRDS(copd_metadata_pre, "data_new/Pre_Exercise/copd_preExc_metadata.RDS")

copd_intersect <- intersect(colnames(copd_data), copd_metadata_pre$seq_sample_id)

copd_data_pre <- copd_data %>%
  subset(select = c("transcript_ID", copd_intersect))

 saveRDS(copd_data_pre, "data_new/Pre_Exercise/copd_preExc_splicing_data.RDS")
