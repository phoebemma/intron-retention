library(dplyr)
library(tidyverse)
library(scales)
library(biomaRt)

source("R/Trainome_functions.R")
 # load the individual metadata 
#COPD metadata
copd_metadata <- readRDS("data/copd_metadata.RDS") %>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)


# Volume_data
Vol_metadata <- readRDS("data/volume_metadata.RDS")%>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age) 


# Contratrain_data
ct_metadata <- readRDS("data/contratrain_metadata.RDS") %>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)


# Publicly available data
SRP102542_metadata <- readRDS("data/SRP102542_metadata.RDS")%>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)


# Alpha and Omega data

A_Omega_metadata <- readRDS("data/Alpha_Omega_metadata.RDS")%>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)



Relief_full_meta <- readRDS("data/Relief_metadata.RDS")%>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)


# Merge the metadata into one
metadata <- rbind(copd_metadata, Vol_metadata)%>%
  rbind(ct_metadata) %>%
  rbind(A_Omega_metadata) %>%
  rbind(SRP102542_metadata) %>%
  rbind(Relief_full_meta) %>%
  mutate(across(c("age"), round, 0)) %>%
  mutate(sex = factor(sex, levels = c("female", "male")),
         time = factor(time, levels = c("PreExc", "PostExc")),
         scaled_age = round(rescale(age), digits = 2),
         participant = paste0(study, "_", participant))




saveRDS(metadata, "data/all_full_metadata.RDS")




# Load the individual splicing data
copd_splice_df <- readRDS("data/copd_splicing_data.RDS")

vol_splice_df <- readRDS("data/volume_splicing_data.RDS")
ct_splice_df <- readRDS("data/contratrain_splicing_data.RDS")
SRP102542_splice_df <- readRDS("data/SRP102542_splicing_data.RDS")
AOD_splice_df <- readRDS("data/Alpha_Omega_splicing_data.RDS")
Relief_full_splice <- readRDS("data/Relief_splicing_data.RDS")





all_splice_df <-copd_splice_df  %>%
  inner_join(vol_splice_df, by = "transcript_ID") %>%
  inner_join(ct_splice_df, by = "transcript_ID") %>%
  inner_join(AOD_splice_df, by = "transcript_ID") %>%
  inner_join(SRP102542_splice_df, by = "transcript_ID") %>%
  inner_join(Relief_full_splice, by = "transcript_ID") 



# select only columnames in the splicing data that match sequence ids in the metadata
intersect <- intersect(colnames(all_splice_df), metadata$seq_sample_id)

all_splice_df <-all_splice_df %>%
  subset(select = c("transcript_ID", intersect)) %>%
  drop_na()


saveRDS(all_splice_df, "data/all_splice.RDS")





# Extract Ensembl gene annotation from BioMart

ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "hsapiens_gene_ensembl",
                   verbose = TRUE)

attributes <- listAttributes(ensembl)

# #extract GENE biotypes and  names
annotation<- getBM(attributes = c("transcript_biotype", "external_transcript_name",	
                                  "ensembl_transcript_id",
                                  "ensembl_gene_id", "ensembl_gene_id_version","ensembl_transcript_id_version", "external_gene_name", "transcript_length"),  mart = ensembl )


saveRDS(annotation, "data/ensembl_gene_annotation.RDS")
