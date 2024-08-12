library(dplyr)
library(tidyverse)
library(readr)
library(plotly)


#Load and clean the first public data
#This contains pre and postexercise (RT) data
#see article https://pubmed.ncbi.nlm.nih.gov/28273480/
#Metadata downloaded from https://www.ncbi.nlm.nih.gov/Traces/study/?page=2&acc=gse97084&o=acc_s%3Aa
SRP102542 <- read_csv("data/SRP102542.txt") %>%
  select("Run", "AGE", "biopsy_timepoint", "exercise_type", "Experiment", "GEO_Accession (exp)", "SRA Study") 


#Load the metadata file downloaded from ENA
#This would be used to match the prexercise to the postexercise data
#This also contains inormation that could be used as participant_id
SRP102542_1 <- read_csv("public_data/sra_result-SRP102542.csv") %>%
  select("Experiment Accession", "Experiment Title", "Sample Accession") %>%
  separate("Experiment Title", c("GEO_Accession (exp)", "Group", "age_group", "Exercise", "Time"))

#JOIN the SRP102542 metadata together
SRP102542 <- SRP102542 %>%
  inner_join(SRP102542_1, by = c("Experiment" = "Experiment Accession", "GEO_Accession (exp)"))
colnames(SRP102542)[colnames(SRP102542) == "SRA Study"] <- "study"

#Extract the participant from the dataframe
SRP102542$participant <- parse_number(SRP102542$Group)

#saveRDS(SRP102542, "./data/SRP102542_metadata.RDS")

colnames(SRP102542)

#Load second public data


#The following is baseline transcriptome metadata.
#See here https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58608 
#Participants were young see https://faseb.onlinelibrary.wiley.com/doi/epdf/10.1096/fj.14-255000
#metadata downloaded from https://www.ncbi.nlm.nih.gov/Traces/study/?page=2&acc=gse58608&o=acc_s%3Aa
SRP043368 <- read_csv("data/SRP043368.txt")%>%
  select("Run", "gender", "exercise_status", "Experiment", "SRA Study") 

#Rename the Experiment to participant and SRA to study
colnames(SRP043368)[colnames(SRP043368) == "SRA Study"] <- "study"
colnames(SRP043368)[colnames(SRP043368) == "Experiment"] <- "participant"

#saveRDS(SRP043368, "data/SRP043368_metadata.RDS")


colnames(SRP043368)
#Load third public data


#Bothe baseline and post-14 weeks training of old adults. To be selected only those with placebo and 
#RT
#Contains only baseline of young individuals

#see article https://pubmed.ncbi.nlm.nih.gov/33071237/
#metadata downloaded from https://www.ncbi.nlm.nih.gov/Traces/study/?page=2&acc=GSE157585&o=acc_s%3Aa

#biopsy 1 is preexercise in old, while biopsy 0 is baseline young
SRP280348 <- read_csv("data/SRP280348.txt")%>%
  select("Run", "biopsy", "study_arm",  "Experiment") 

#Differentiates data into pre and post
#Source ENA
SRP280348_1 <- read_csv("public_data/sra_result-SRP280348.csv")%>%
  select("Experiment Accession", "Experiment Title", "Sample Accession", "Study Accession") %>%
  separate("Experiment Title", c("Title", "Group", "Category"), sep = " ") #%>%
  #separate("Group", c("x", "participant", "y" ), remove = FALSE) %>%
  #subset(select = (-c(x, y)))

SRP280348 <- SRP280348 %>%
  inner_join(SRP280348_1, by = c( "Experiment" = "Experiment Accession")) %>%
  mutate(age_group = case_when(biopsy == 0 ~ "young", biopsy == 1 ~ "old",biopsy == 3 ~ "old" ))

colnames(SRP280348)[colnames(SRP280348) == "Study Accession"] <- "study"
#saveRDS(SRP280348, "data/SRP280348_metadata.RDS")

colnames(SRP280348)

unique(SRP280348$biopsy)
