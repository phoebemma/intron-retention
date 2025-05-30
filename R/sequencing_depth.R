library(dplyr)
library(tidyverse)
library(stringr)

# calculating the average sequencing depth of each dataset
# Formula (Total reads * average read length)/ refence genome size(3.2 billion)

# Load the alpha/Omega STAR-MultiQC summary
# the depths are divided by one million

# Load metadata so as to match samplenames
all_metadata <- readRDS("data_new/processed_data/all_full_metadata.RDS")


AO_summary <- read_csv("data_new/Multiqc/AO_star_summary_table.csv") %>%
  dplyr::select( "Sample" , "Total reads", "Uniq aligned...5", "Avg. read len") %>%
  mutate(Sample = sub("(.*?)_EKRN.*", "\\1", Sample))


# extract sampes that match samplenames in metadata

AO_summary <- AO_summary %>%
  filter(Sample %in% all_metadata$seq_sample_id)

AO_summary$sequencing_depth <- (AO_summary$`Total reads` * AO_summary$`Avg. read len`) / 3200 # size GRch38 divided by million

AO_depth <- mean(AO_summary$sequencing_depth)


# Contratrain
Ct_summary <- read_csv("data_new/Multiqc/Contratrain_star_summary_table.csv")%>%
  dplyr::select( "Sample" , "Total reads", "Uniq aligned...5", "Avg. read len")

Ct_summary$sequencing_depth <- (Ct_summary$`Total reads` * Ct_summary$`Avg. read len`) / 3200

Ct_depth <- mean(Ct_summary$sequencing_depth)



# COPD
COPD_summary <- read_csv("data_new/Multiqc/COPD_star_summary_table.csv")%>%
  dplyr::select( "Sample" , "Total reads", "Uniq aligned...5", "Avg. read len")%>%
  mutate(Sample = sub("(.*?)_R1.*", "\\1", Sample))
all_metadata$seq_sample_id <- sub("^X", "", all_metadata$seq_sample_id)
COPD_summary <- COPD_summary %>%
  filter(Sample %in% all_metadata$seq_sample_id)


COPD_summary$sequencing_depth <- (COPD_summary$`Total reads` * COPD_summary$`Avg. read len`) / 3200

COPD_depth <- mean(COPD_summary$sequencing_depth)


# Relief
Relief_summary <- read_csv("data_new/Multiqc/Relief_star_summary_table.csv")%>%
  dplyr::select( "Sample" , "Total reads", "Uniq aligned...5", "Avg. read len") %>%
  mutate(Sample = sub("(.*?)_EKRN.*", "\\1", Sample))

Relief_summary <- Relief_summary %>%
  filter(Sample %in% all_metadata$seq_sample_id)

Relief_summary$sequencing_depth <- (Relief_summary$`Total reads` * Relief_summary$`Avg. read len`) / 3200

Relief_depth <- mean(Relief_summary$sequencing_depth)



# Robinson et al (2017)
SRP102542_summary <- read_csv("data_new/Multiqc/SRP102542_star_summary_table.csv")%>%
  dplyr::select( "Sample" , "Total reads", "Uniq aligned...5", "Avg. read len")%>%
  mutate(Sample = sub("(.*?)_1", "\\1", Sample))
SRP102542_summary <- SRP102542_summary %>%
  filter(Sample %in% all_metadata$seq_sample_id)

SRP102542_summary$sequencing_depth <- (SRP102542_summary$`Total reads` * SRP102542_summary$`Avg. read len`) / 3200

SRP102542_depth <- mean(SRP102542_summary$sequencing_depth)


# Volume
Volume_summary <- read_csv("data_new/Multiqc/Volume_star_summary_table.csv")%>%
  dplyr::select( "Sample" , "Total reads", "Uniq aligned...5", "Avg. read len")

Volume_summary$sequencing_depth <- (Volume_summary$`Total reads` * Volume_summary$`Avg. read len`) / 3200

Volume_depth <- mean(Volume_summary$sequencing_depth)



