library(dplyr)
library(tidyverse)
library(readr)

#Load the trainome functions file
source("R/Trainome_functions.R")
#Load public data


#Bothe baseline and post-14 weeks training of old adults. To be selected only those with placebo and 
#RT
#biopsy 1 is preexercise in old, while biopsy 0 is baseline young

#see article https://pubmed.ncbi.nlm.nih.gov/33071237/
#metadata downloaded from https://www.ncbi.nlm.nih.gov/Traces/study/?page=2&acc=GSE157585&o=acc_s%3Aa


SRP280348 <- read_csv("data/SRP280348.txt")%>%
  dplyr::select("Run", "biopsy", "study_arm",  "Experiment") 

#Differentiates data into pre and post
#Source ENA
SRP280348_1 <- read_csv("public_data/sra_result-SRP280348.csv")%>%
  dplyr::select("Experiment Accession", "Experiment Title", "Sample Accession", "Study Accession") %>%
  separate("Experiment Title", c("Title", "Group", "Category"), sep = " ") #%>%
#separate("Group", c("x", "participant", "y" ), remove = FALSE) %>%
#subset(select = (-c(x, y)))

SRP280348 <- SRP280348 %>%
  inner_join(SRP280348_1, by = c( "Experiment" = "Experiment Accession")) %>%
  mutate(age_group = case_when(biopsy == 0 ~ "Young", biopsy == 1 ~ "Old",biopsy == 3 ~ "Old" ))

colnames(SRP280348)[colnames(SRP280348) == "Study Accession"] <- "study"
colnames(SRP280348)[colnames(SRP280348) == "Run"] <- "seq_sample_id"
colnames(SRP280348)[colnames(SRP280348) == "Sample Accession"] <- "participant"
colnames(SRP280348)[colnames(SRP280348) == "biopsy"] <- "time"

SRP280348$time <- as.character(SRP280348$time)

SRP280348["time"][SRP280348["time"] == 1 ] <- "PreExc"

SRP280348["time"][SRP280348["time"] == 0 ] <- "PreExc"
SRP280348["time"][SRP280348["time"] == 3 ] <- "PostExc"

unique(SRP280348_metadata$time)

saveRDS(SRP280348, "data/processed_data/SRP280348_metadata.RDS")

#saveRDS(SRP280348, "data/processed_data/SRP280348_metadata.RDS")

colnames(SRP280348)




#Load splicing data
#Contains the baseline and post-14 weeks training of old adults. To be selected only those with placebo and 
#RT
#Also contains only baseline of young individuals

#see article https://pubmed.ncbi.nlm.nih.gov/33071237/
SRP280348_data <- extract_splice_q("./data/SRP280348_GSE157585_SpliceQ_outputs/")
idx <- sapply(SRP280348_data, class)== "numeric"
SRP280348_data[, idx] <- lapply(SRP280348_data[, idx], round, 2)


#select only the splicing samples captured in the metadata
SRP280348_intersect <- intersect(colnames(SRP280348_data), SRP280348$seq_sample_id)

SRP280348_data <- SRP280348_data %>%
  subset(select = c("transcript_ID", SRP280348_intersect))


#saveRDS(SRP280348_data, "./data/processed_data/SRP280348_splicing_data.RDS")


#Subset to preexercise data

SRP280348_pre_met <- SRP280348 %>%
  subset(time == "PreExc")
#saveRDS(SRP280348_pre_met, "data/preexercise_data/SRP280348_preExc_metadata.RDS")


#select only the splicing samples captured in the metadata
SRP280348_intersect <- intersect(colnames(SRP280348_data), SRP280348_pre_met$seq_sample_id)

SRP280348_pre_data <- SRP280348_data %>%
  subset(select = c("transcript_ID", SRP280348_intersect))

#saveRDS(SRP280348_pre_data, "data/preexercise_data/SRP280348_preExc_splicing_data.RDS")
