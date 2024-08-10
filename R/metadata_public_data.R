library(dplyr)
library(tidyverse)


#Load and clean the first public data
#This contains pre and postexercise (RT) data
#see article https://pubmed.ncbi.nlm.nih.gov/28273480/
#Metadata downloaded from https://www.ncbi.nlm.nih.gov/Traces/study/?page=2&acc=gse97084&o=acc_s%3Aa
SRP102542 <- read_csv("data/SRP102542.txt") %>%
  select("Run", "AGE", "biopsy_timepoint", "exercise_type") 
#saveRDS(SRP102542, "./data/SRP102542_metadata.RDS")

colnames(SRP102542)

#Load second public data


#The following is baseline transcriptome metadata.
#See here https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58608 
#Participants were young see https://faseb.onlinelibrary.wiley.com/doi/epdf/10.1096/fj.14-255000
#metadata downloaded from https://www.ncbi.nlm.nih.gov/Traces/study/?page=2&acc=gse58608&o=acc_s%3Aa
SRP043368 <- read_csv("data/SRP043368.txt")%>%
  select("Run", "gender", "exercise_status") 
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
  select("Run", "biopsy", "study_arm") 

#saveRDS(SRP280348, "data/SRP280348_metadata.RDS")

colnames(SRP280348)

unique(SRP280348$biopsy)
