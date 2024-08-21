library(dplyr)
library(tidyverse)
library(readr)

#Load the trainome functions file
source("R/Trainome_functions.R")



#The following is baseline transcriptome metadata.
#See here https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58608 
#Participants were young see https://faseb.onlinelibrary.wiley.com/doi/epdf/10.1096/fj.14-255000
#metadata downloaded from https://www.ncbi.nlm.nih.gov/Traces/study/?page=2&acc=gse58608&o=acc_s%3Aa
SRP043368 <- read_csv("data/SRP043368.txt")%>%
  dplyr::select("Run", "gender", "exercise_status", "Experiment", "SRA Study") 

#Rename the Experiment to participant and SRA to study
colnames(SRP043368)[colnames(SRP043368) == "SRA Study"] <- "study"
colnames(SRP043368)[colnames(SRP043368) == "Experiment"] <- "participant"
colnames(SRP043368)[colnames(SRP043368) == "Run"] <- "seq_sample_id"

colnames(SRP043368)[colnames(SRP043368) == "exercise_status"] <- "time"
SRP043368["time"][SRP043368["time"] == "untrained" ] <- "PreExc"
#saveRDS(SRP043368, "data/preexercise_data/SRP043368_metadata.RDS")


colnames(SRP043368)

SRP043368_data <- extract_splice_q("./data/SRP043368_GSE58608_SpliceQ_outputs/")
idx <- sapply(SRP043368_data, class)== "numeric"
SRP043368_data[, idx] <- lapply(SRP043368_data[, idx], round, 2)

#saveRDS(SRP043368_data, "./data/preexercise_data/SRP043368_splicing_data.RDS")





