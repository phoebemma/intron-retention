#Load the file with libraries
#source("R/libraries.R")
library(dplyr)
library(tidyverse)
library(seqwrap)

#Load the trainome functions file
source("R/Trainome_functions.R")


#Load the metadata of the the datasets

#Copd 
copd_metadata <- readRDS("data/processed_data/copd_metadata.RDS")

unique(copd_metadata$time)

#Prexercise in COPD is named "PreExc"

PreEXC_copd_meta <- copd_metadata %>%
  subset(time == "PreExc") %>%
  drop_na()

#saveRDS(PreEXC_copd_meta, "data/preexercise_data/preexc_copd_metadata.RDS")



#Volume

volume_metadata <- readRDS("data/processed_data/volume_metadata.RDS")

unique(volume_metadata$time)

#Prexercise data in volume is "w0" 

PreExc_vol_meta <- volume_metadata %>%
  subset(time == "w0") %>%
  drop_na()

#saveRDS(PreExc_vol_meta, "data/preexercise_data/preexc_volume_metadata.RDS")

#Contratrain
ct_metadata <- readRDS("data/processed_data/contratrain_metadata.RDS")

unique(ct_metadata$time)

#Subset to the Pre-exercise data

PreExc_ct_meta <- ct_metadata %>%
  subset(time == "t1") %>%
  drop_na()

#SRP280348 

SRP280348_meta <- readRDS("data/processed_data/SRP280348_metadata.RDS")

colnames(SRP280348_meta)
#Subset the preexercise  data
#This is the biopsy data 1 for old participants, and 0 for young participants

PreExc_SRP280348_meta <- SRP280348_meta %>%
  subset(biopsy == 1 | biopsy == 0)

#saveRDS(PreExc_SRP280348_meta, "data/preexercise_data/preexc_SRP280348_metadata.RDS")



#SRP043368 is baseline data alone, would be loaded from the processed data folder

SRP043368 <- readRDS("data/processed_data/SRP043368_metadata.RDS")



#SRP102542

SRP102542_meta <- readRDS("data/processed_data/SRP102542_metadata.RDS")

PreExc_SRP102542_meta <- SRP102542 %>%
  subset(biopsy_timepoint == "PreTraining")%>%
  drop_na()

#saveRDS(PreExc_SRP102542_meta, "data/preexercise_data/preexc_SRP102542_metadata.RDS")

# Post_ex_df <- full_splice_df %>%
#   subset(time == "PostExc")



#Subset only the introns wih SE 0

fully_retained_ints <- full_splice_df %>%
  subset(SE == 0)


table(fully_retained_ints$age)
table(fully_retained_ints$time)


#Plot the data
plot(fully_retained_ints$time)
