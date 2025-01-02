library(dplyr)
library(tidyverse)
library(readr)
library(stringr)

source("R/Trainome_functions.R")
# This details the extractions and EDA fror the dataset with SRA number SRP102542
# Load and clean the publicly available data
# This contains pre and postexercise (RT) data
# see article https://pubmed.ncbi.nlm.nih.gov/28273480/
# Sequence data downloaded from https://www.ebi.ac.uk/ena/browser/view/PRJNA380649
# Metadata downloaded from https://www.ncbi.nlm.nih.gov/Traces/study/?page=2&acc=gse97084&o=acc_s%3Aa
SRP102542 <- read_csv("data/SRP102542.txt") %>%
  dplyr::select("Run", "AGE", "biopsy_timepoint", "exercise_type", "Experiment", "GEO_Accession (exp)", "SRA Study") 


# Load the metadata file downloaded from ENA
# This would be used to match the prexercise to the postexercise data
# This also contains inormation that could be used as participant_id
SRP102542_1 <- read_csv("public_data/sra_result-SRP102542.csv") %>%
  dplyr::select("Experiment Accession", "Experiment Title", "Sample Accession") %>%
  separate("Experiment Title", c("GEO_Accession (exp)", "Group", "age_group", "Exercise", "biopsy_timepoint"))

# Join the SRP102542 metadata together
SRP102542 <- SRP102542 %>%
  inner_join(SRP102542_1, by = c("Experiment" = "Experiment Accession", "GEO_Accession (exp)", "biopsy_timepoint"))
colnames(SRP102542)[colnames(SRP102542) == "SRA Study"] <- "study"
# rename "Run" to "seq_sample_id"
# This aligns it to the column name in the mtadata
colnames(SRP102542)[colnames(SRP102542) == "Run"] <- "seq_sample_id"
# colnames(SRP102542)[colnames(SRP102542) == "Experiment"] <- "participant"
colnames(SRP102542)[colnames(SRP102542) == "biopsy_timepoint"] <- "time"
SRP102542["time"][SRP102542["time"] == "PreTraining" ] <- "PreExc"
SRP102542["time"][SRP102542["time"] == "PostTraining" ] <- "PostExc"
# Extract the participant number which is unique to each participant
SRP102542["participant"] <- as.character(str_extract(SRP102542$Group, "[0-9]+"))



# Load the data obtained from the owners that include age and sex information. 

extra_data <- readxl::read_xlsx("public_data/Robinson Nair Exercise Tissue Code for Hammarstrom.xlsx")%>%
  dplyr::select("Tissue Code", "Age Group", "Sex", "Age")
extra_data$`Tissue Code` <- as.character(extra_data$`Tissue Code`)

SRP102542<- SRP102542 %>%
  inner_join(extra_data, by = c( "participant" = "Tissue Code" ))

# Rename the columns to match those in the primary data
colnames(SRP102542)[colnames(SRP102542) == "Age"] <- "age"
colnames(SRP102542)[colnames(SRP102542) == "Sex"] <- "sex"
colnames(SRP102542)[colnames(SRP102542) == "exercise_type"] <- "condition"
colnames(SRP102542)
unique(SRP102542$time)

length(unique(extra_data$`Tissue Code`))
SRP102542 <- SRP102542 %>%
  dplyr::select(study, participant, sex, condition, time, seq_sample_id, age, age_group)

SRP102542["sex"][SRP102542["sex"] == "M" ] <- "male"

SRP102542["sex"][SRP102542["sex"] == "F" ] <- "female"



length(unique(SRP102542$participant))
#Select pre-exercise data for the pre-exercise model

SRP102542_pre <- SRP102542 %>%
  subset(time == "PreExc")
 saveRDS(SRP102542_pre, "data_new/Pre_Exercise/SRP102542_preExc_metadata.RDS")
 
 # For the full data analyses, we would need only participants that performed RT

 SRP102542_full <- SRP102542 %>%
   filter(condition == "Resistance")
#  saveRDS(SRP102542_full, "./data_new/processed_data/SRP102542_metadata.RDS")


#Load the SpliceQ data
SRP102542_data <- extract_splice_q("./data_new/SRP102542_SpliceQ_outputs/")

idx <- sapply(SRP102542_data, class)== "numeric"
SRP102542_data[, idx] <- lapply(SRP102542_data[, idx], round, 2)




#select only the splicing samples captured in the metadata
SRP102542_intersect <- intersect(colnames(SRP102542_data), SRP102542$seq_sample_id)

SRP102542_data <- SRP102542_data %>%
  subset(select = c("transcript_ID", SRP102542_intersect))
# saveRDS(SRP102542_data, "./data_new/processed_data/SRP102542_splicing_data.RDS")


#Repeat for pre-exercise data
#select only the splicing samples captured in the preexercise metadata
SRP102542_intersect <- intersect(colnames(SRP102542_data), SRP102542_pre$seq_sample_id)

SRP102542_pre_data <- SRP102542_data %>%
  subset(select = c("transcript_ID", SRP102542_intersect))
# saveRDS(SRP102542_pre_data, "./data_new/Pre_Exercise/SRP102542_preExc_splicing_data.RDS")


