source("R/Trainome_functions.R")
library(dplyr)
library(tidyverse)
library(ggplot2)
library(trainomeMetaData)
library(edgeR)




alpha_omega_counts <- extract_rsem_isoform_counts("data_new/gene_counts/Alpha_Omega_RSEM_outputs/")

# load the alpha omega metadata 
alpha_omega_metadata <- readRDS("data_new/processed_data/Alpha_Omega_metadata.RDS")

alpha_omega <- alpha_omega_counts %>% 
   separate(transcript_id, c("transcript_id", "transcript_name"), sep = "_", extra = "merge") %>%
  
  # drop transcript_id, select transcripte_name and any of the sample names that match sample name in metadata
  dplyr::select(transcript_id, all_of(alpha_omega_metadata$seq_sample_id))


# Create dge lists and calculate norm factors
dge_ao   <- DGEList(alpha_omega[,-1])


dge_ao  <- calcNormFactors(dge_ao, method = "TMM")


# get the normalised counts

normalised_ao_counts <- as.data.frame(cpm(dge_ao, normalized.lib.sizes = T))

normalised_ao_counts$transcript_id <- alpha_omega$transcript_id

# Make the transcript_ID the first column

normalised_ao_counts <- normalised_ao_counts %>%
  dplyr::select(transcript_id, everything())


saveRDS(normalised_ao_counts, "data_new/gene_counts/RSEM_cpm_normalised_AO_counts.RDS")




# Repeat normalisation step for Volume data


Vol_counts <- extract_rsem_isoform_counts("data_new/gene_counts/Volume_RSEM_outputs_new/")
#remove everything before the . in sample_id. 
colnames(Vol_counts) <- gsub(".*?\\.", "", colnames(Vol_counts) )


# load the volume metadata 
volume_metadata <- readRDS("data_new/processed_data/volume_metadata.RDS")



volume <- Vol_counts %>% 
   separate(transcript_id, c("transcript_id", "transcript_name"), sep = "_", extra = "merge") %>%
  
  # drop transcript_id, select transcripte_name and any of the sample names that match sample name in metadata
  dplyr::select(transcript_id, all_of(volume_metadata$seq_sample_id))







# Create dge lists and calculate norm factors
dge_vol   <- DGEList(volume[,-1])


dge_vol  <- calcNormFactors(dge_vol, method = "TMM")


# get the normalised counts

normalised_vol_counts <- as.data.frame(cpm(dge_vol, normalized.lib.sizes = T))

normalised_vol_counts$transcript_id<- volume$transcript_id

# Make the transcript_ID the first column

normalised_vol_counts <- normalised_vol_counts %>%
  dplyr::select(transcript_id, everything())


saveRDS(normalised_vol_counts, "data_new/gene_counts/RSEM_cpm_normalised_vol_counts.RDS")






# Contratrain data

Ct_counts <- extract_rsem_isoform_counts("data_new/gene_counts/Contratrain_RSEM_outputs_new/")


# load the volume metadata 
ct_metadata <- readRDS("data_new/processed_data/contratrain_metadata.RDS")



Contratrain <- Ct_counts %>% 
   separate(transcript_id, c("transcript_id", "transcript_name"), sep = "_", extra = "merge") %>%
  
  # drop transcript_id, select transcripte_name and any of the sample names that match sample name in metadata
  dplyr::select(transcript_id, all_of(ct_metadata$seq_sample_id))




# Create dge lists and calculate norm factors
dge_ct   <- DGEList(Contratrain[,-1])


dge_ct  <- calcNormFactors(dge_ct, method = "TMM")


# get the normalised counts

normalised_ct_counts <- as.data.frame(cpm(dge_ct, normalized.lib.sizes = T))

normalised_ct_counts$transcript_id <- Contratrain$transcript_id

# Make the transcript_ID the first column

normalised_ct_counts <- normalised_ct_counts %>%
  dplyr::select(transcript_id, everything())


saveRDS(normalised_ct_counts, "data_new/gene_counts/RSEM_cpm_normalised_contratrain_counts.RDS")




# COPD


Copd_counts <-  extract_rsem_isoform_counts("data_new/gene_counts/COPD_RSEM_outputs_new/")

# load COPD metadata
copd_metadata <- readRDS("data_new/processed_data/copd_metadata.RDS")


Copd <- Copd_counts %>% 
   separate(transcript_id, c("transcript_id", "transcript_name"), sep = "_", extra = "merge") %>%
  
  # drop transcript_id, select transcripte_name and any of the sample names that match sample name in metadata
  dplyr::select(transcript_id, all_of(copd_metadata$seq_sample_id))




# Create dge lists and calculate norm factors
dge_copd   <- DGEList(Copd[,-1])


dge_copd  <- calcNormFactors(dge_copd, method = "TMM")


# get the normalised counts

normalised_copd_counts <- as.data.frame(cpm(dge_copd, normalized.lib.sizes = T))

normalised_copd_counts$transcript_id <- Copd$transcript_id

# Make the transcript_ID the first column

normalised_copd_counts <- normalised_copd_counts %>%
  dplyr::select(transcript_id, everything())


saveRDS(normalised_copd_counts, "data_new/gene_counts/RSEM_cpm_normalised_COPD_counts.RDS")









Relief_counts <- extract_rsem_isoform_counts("data_new/gene_counts/Relief_RSEM_outputs/")


# Load the Relief metadata
Relief_metadata <- readRDS("data_new/processed_data/Relief_metadata.RDS")


Relief <- Relief_counts %>% 
   separate(transcript_id, c("transcript_id", "transcript_name"), sep = "_", extra = "merge") %>%
  
  # drop transcript_id, select transcripte_name and any of the sample names that match sample name in metadata
  dplyr::select(transcript_id, all_of(Relief_metadata$seq_sample_id))




# Create dge lists and calculate norm factors
dge_Relief   <- DGEList(Relief[,-1])


dge_Relief  <- calcNormFactors(dge_Relief, method = "TMM")


# get the normalised counts

normalised_Relief_counts <- as.data.frame(cpm(dge_Relief, normalized.lib.sizes = T))

normalised_Relief_counts$transcript_id <- Relief$transcript_id

# Make the transcript_ID the first column

normalised_Relief_counts <- normalised_Relief_counts %>%
  dplyr::select(transcript_id, everything())


saveRDS(normalised_Relief_counts, "data_new/gene_counts/RSEM_cpm_normalised_Relief_counts.RDS")







# SRP102542 data

SRP102542_counts <- extract_rsem_isoform_counts("data_new/gene_counts/SRP102542_RSEM_outputs/")




# Load the Relief metadata
SRP102542_metadata <- readRDS("data_new/processed_data/SRP102542_metadata.RDS")


SRP102542 <- SRP102542_counts %>% 
  separate(transcript_id, c("transcript_id", "transcript_name"), sep = "_", extra = "merge") %>%
  
  # drop transcript_id, select transcripte_name and any of the sample names that match sample name in metadata
  dplyr::select(transcript_id, all_of(SRP102542_metadata$seq_sample_id))



# Create dge lists and calculate norm factors
dge_SRP102542   <- DGEList(SRP102542[,-1])


dge_SRP102542  <- calcNormFactors(dge_SRP102542, method = "TMM")


# get the normalised counts

normalised_SRP102542_counts <- as.data.frame(cpm(dge_SRP102542, normalized.lib.sizes = T))

normalised_SRP102542_counts$transcript_id <- SRP102542$transcript_id

# Make the transcript_ID the first column

normalised_SRP102542_counts <- normalised_SRP102542_counts %>%
  dplyr::select(transcript_id, everything())


saveRDS(normalised_SRP102542_counts, "data_new/gene_counts/RSEM_cpm_normalised_SRP102542_counts.RDS")


###### Batch effects correction
# Merge the normalised gene counts into one

all_normalised_transcripts <- normalised_copd_counts %>%
  full_join(normalised_vol_counts, by = "transcript_id") %>%
  full_join(normalised_ct_counts, by = "transcript_id") %>%
  full_join(normalised_ao_counts, by = "transcript_id") %>%
  full_join(normalised_SRP102542_counts, by = "transcript_id") %>%
  full_join(normalised_Relief_counts, by = "transcript_id") %>%
  drop_na() 




all_metadata <- readRDS("data_new/processed_data/all_full_metadata.RDS")

# Check if everything matches except the transcript_id
match(colnames(all_normalised_transcripts), all_metadata$seq_sample_id)

all_transcripts_reordered <- all_normalised_transcripts[ , c("transcript_id",all_metadata$seq_sample_id)] 

all_transcripts_reordered[,-1] <- round(all_transcripts_reordered[,-1], 2)

# recheck# recheckround()
match(colnames(all_transcripts_reordered), all_metadata$seq_sample_id)

#saveRDS(all_transcripts_reordered, "data_new/gene_counts/all_normalised_counts.RDS")
# all_transcripts_reordered <- readRDS("data_new/gene_counts/all_normalised_counts.RDS")
saveRDS(all_transcripts_reordered, "data_new/gene_counts/all_normalised_RSEM_counts.RDS")

# Define the batch, which in this case is the styd
batch <- factor(all_metadata$study)

# correct batch effect
corrected_counts <- removeBatchEffect(all_transcripts_reordered[,-1], batch = batch, group = all_metadata$time)
corrected_counts$transcript_ID <- all_transcripts_reordered$transcript_id

#boxplot(corrected_counts, main = "batch corrected data")

saveRDS(corrected_counts, "data_new/gene_counts/batch_corrected_RSEM_counts")
