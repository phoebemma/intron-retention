library(dplyr)
library(tidyverse)


#Load and clean the first public data
#This contains pre and postexercise (RT) data
#see article https://pubmed.ncbi.nlm.nih.gov/28273480/
SRP102542 <- read_csv("public_data/sra_result-SRP102542.csv") %>%
  select("Experiment Accession", "Experiment Title", "Sample Accession") %>%
  separate("Experiment Title", c("Title", "Group", "Age", "Exercise", "Time"))


unique(SRP102542$Time)
#Load second public data


#The following is baseline transcriptome metadata. But does not specifiy age
#See here https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58608 
#Participants were young see https://faseb.onlinelibrary.wiley.com/doi/epdf/10.1096/fj.14-255000
SRP043368 <- read_csv("public_data/sra_result-SRP043368.csv")%>%
  select("Experiment Accession", "Experiment Title", "Sample Accession") %>%
  separate("Experiment Title", c("Title", "Group"))


unique(SRP043368$`Experiment Title`)
#Load third public data


#Bothe baseline and post-14 weeks training of old adults. To be selected only those with placebo and 
#RT
#Contains only baseline of young individuals

#see article https://pubmed.ncbi.nlm.nih.gov/33071237/

SRP280348 <- read_csv("public_data/sra_result-SRP280348.csv")%>%
  select("Experiment Accession", "Experiment Title", "Sample Accession") %>%
  separate("Experiment Title", c("Title", "Group", "Age"), sep = " ")

SRP280348$Age
