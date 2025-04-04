# This script contains the extraction and processing of the gene expression data from the 6 studies
source("R/Trainome_functions.R")
library(dplyr)
library(tidyverse)
library(ggplot2)
library(trainomeMetaData)
library(edgeR)

# Download the gene level data from trainomeMetaData
# Copd
# download_ome(download = "copd_gene_rsem")
# 
# # Contratrain
# download_ome(download = "ct_gene_rsem")
# 
# # Volume
# download_ome(download = "vol_gene_rsem")
# 
# 
# # Alpha/Omega
# download_ome(download = "alphaomega_gene_rsem")
 


# load the alpha omega metadata 
alpha_omega_metadata <- readRDS("data_new/processed_data/Alpha_Omega_metadata.RDS")


# Read the AO counts
# This uses a function created for extracting RSEM output as expected counts
alpha_omega <-  extract_rsem_gene_counts("data_new/gene_counts/Alpha_Omega_RSEM_outputs/") %>%
  separate(gene_id, c("gene_id", "gene_name"), sep = "_", extra = "merge") %>%
  
  # select gene_name and any of the sample names that match sample name in metadata
  dplyr::select(gene_id, all_of(alpha_omega_metadata$seq_sample_id))



# Create dge lists and calculate norm factors
dge_ao   <- DGEList(alpha_omega[,-1])


dge_ao  <- calcNormFactors(dge_ao, method = "TMM")


# get the normalised counts

normalised_ao_counts <- as.data.frame(cpm(dge_ao, normalized.lib.sizes = T))

normalised_ao_counts$gene_id <- alpha_omega$gene_id

# Make the gene_Id the first column

normalised_ao_counts <- normalised_ao_counts %>%
  dplyr::select(gene_id, everything())




# Repeat normalisation step for Volume data

volume <- extract_rsem_gene_counts("data_new/gene_counts/Volume_RSEM_outputs_new/") 
#remove everything before the . in sample_id. 
colnames(volume) <- gsub(".*?\\.", "", colnames(volume) )

# load the volume metadata 
volume_metadata <- readRDS("data_new/processed_data/volume_metadata.RDS") 


volume <- volume %>%
  separate(gene_id, c("gene_id", "gene_name"), sep = "_", extra = "merge") %>%
  dplyr::select(gene_id, all_of(volume_metadata$seq_sample_id))

dge_vol   <- DGEList(volume[,-1])


dge_vol  <- calcNormFactors(dge_vol, method = "TMM")


# get the normalised counts

normalised_vol_counts <- as.data.frame(cpm(dge_vol, normalized.lib.sizes = T))

normalised_vol_counts$gene_id <- volume$gene_id

normalised_vol_counts <- normalised_vol_counts %>%
  dplyr::select(gene_id, everything())



#  Repeat for Contratrain data

ct_metadata <- readRDS("data_new/processed_data/contratrain_metadata.RDS")

Contratrain <- extract_rsem_gene_counts("data_new/gene_counts/Contratrain_RSEM_outputs_new/") %>%
  separate(gene_id, c("gene_id", "gene_name"), sep = "_", extra = "merge") %>%
  dplyr::select(gene_id, all_of(ct_metadata$seq_sample_id)) 


dge_ct   <- DGEList(Contratrain[,-1])


dge_ct  <- calcNormFactors(dge_ct, method = "TMM")

normalised_ct_counts <- as.data.frame(cpm(dge_ct, normalized.lib.sizes = T))

normalised_ct_counts$gene_id <- Contratrain$gene_id

normalised_ct_counts <- normalised_ct_counts %>%
  dplyr::select(gene_id, everything())




# COPD data normalisation

# load COPD metadata
copd_metadata <- readRDS("data_new/processed_data/copd_metadata.RDS")


Copd <-  extract_rsem_gene_counts("data_new/gene_counts/COPD_RSEM_outputs_new/") %>%

  separate(gene_id, c("gene_id", "gene_name"), sep = "_", extra = "merge") %>%
  
  dplyr::select(gene_id, all_of(copd_metadata$seq_sample_id))


dge_copd   <- DGEList(Copd[,-1])


dge_copd  <- calcNormFactors(dge_copd, method = "TMM")

normalised_copd_counts <- as.data.frame(cpm(dge_copd, normalized.lib.sizes = T))

normalised_copd_counts$gene_id <- Copd$gene_id


normalised_copd_counts <- normalised_copd_counts %>%
  dplyr::select(gene_id, everything())


# Relief data
# Load the Relief metadata
Relief_metadata <- readRDS("data_new/processed_data/Relief_metadata.RDS")

Relief <- extract_rsem_gene_counts("data_new/gene_counts/Relief_RSEM_outputs/") %>%
  separate(gene_id, c("gene_id", "gene_name"), sep = "_", extra = "merge") %>%
  dplyr::select(gene_id, all_of(Relief_metadata$seq_sample_id))

dge_Relief   <- DGEList(Relief[,-1])


dge_Relief  <- calcNormFactors(dge_Relief, method = "TMM")

normalised_Relief_counts <- as.data.frame(cpm(dge_Relief, normalized.lib.sizes = T))

normalised_Relief_counts$gene_id <- Relief$gene_id

normalised_Relief_counts <- normalised_Relief_counts %>%
  dplyr::select(gene_id, everything())



# SRP102542 data
# Load the  metadata
SRP102542_metadata <- readRDS("data_new/processed_data/SRP102542_metadata.RDS")


SRP102542 <- extract_rsem_gene_counts("data_new/gene_counts/SRP102542_RSEM_outputs/") %>%
  separate(gene_id, c("gene_id", "gene_name"), sep = "_", extra = "merge") %>%
 
  dplyr::select(gene_id, all_of(SRP102542_metadata$seq_sample_id))

dge_SRP102542   <- DGEList(SRP102542[,-1])


dge_SRP102542  <- calcNormFactors(dge_SRP102542, method = "TMM")

normalised_SRP102542_counts <- as.data.frame(cpm(dge_SRP102542, normalized.lib.sizes = T))

normalised_SRP102542_counts$gene_id <- SRP102542$gene_id


normalised_SRP102542_counts <- normalised_SRP102542_counts %>%
  dplyr::select(gene_id, everything())



###### Batch effects correction
 # Merge the normalised gene counts into one

all_normalised_transcripts <- normalised_copd_counts %>%
  full_join(normalised_vol_counts, by = "gene_id") %>%
  full_join(normalised_ct_counts, by = "gene_id") %>%
  full_join(normalised_ao_counts, by = "gene_id") %>%
  full_join(normalised_SRP102542_counts, by = "gene_id") %>%
  full_join(normalised_Relief_counts, by = "gene_id") 

# round to two decimal places
all_normalised_transcripts[,-1] <- round(all_normalised_transcripts[,-1], 2)


all_metadata <- readRDS("data_new/processed_data/all_full_metadata.RDS")

# Check if everything matches except the transcript_id
match(colnames(all_normalised_transcripts), all_metadata$seq_sample_id)

# If it doesnt, uncomment and run the next code
# all_transcripts_reordered <- all_normalised_transcripts[ , c("transcript_ID",all_metadata$seq_sample_id)] 

# recheck# recheckround()
# match(colnames(all_transcripts_reordered), all_metadata$seq_sample_id)

# save normalised gene counts
 saveRDS(all_normalised_transcripts, "data_new/gene_counts/all_normalised_gene_counts.RDS")
 
 
 # Remove rows with all zeros
 nonzero <- all_normalised_transcripts %>%
      dplyr::filter(rowSums(all_normalised_transcripts[,-1]) != 0)
 
 # save 
 saveRDS(nonzero, "data_new/gene_counts/nonzero_rows_normalised_genecounts.RDS")

 
 # Define the batch, which in this case is the study
batch <- factor(all_metadata$study)

# correct batch effect
corrected_counts <- as.data.frame(removeBatchEffect(nonzero[,-1], batch = batch, group = all_metadata$time))

corrected_counts$gene_id <-  nonzero$gene_id

#bring the gene counts to the fore

corrected_counts <- corrected_counts %>%
  dplyr::select(gene_id, everything())

saveRDS(corrected_counts, "data_new/gene_counts/batch_corrected_genecounts.RDS")
