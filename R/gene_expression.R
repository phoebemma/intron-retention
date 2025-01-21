source("R/Trainome_functions.R")
library(dplyr)
library(tidyverse)
library(ggplot2)
library(trainomeMetaData)
library(edgeR)

# Download the datasets available on the trainomemetadata
# would be using the transcript level data since the splicing information is contained at the transcript level



# download the counts available on trainomeMetaData 
# download_ome(download = "vol_transcript_rsem")
# download_ome(download = "ct_transcript_rsem")
# download_ome(download = "copd_transcript_rsem")
# download_ome(download = "alphaomega_transcript_rsem")

# to check available files
download_ome(download = "none")



# Vol_counts <- read_csv("ome-data/vol_transcript_rsem.csv")
# 
# Ct_counts <- read_csv("ome-data/ct_transcript_rsem.csv")
# 
# Copd_counts <- read_csv("ome-data/copd_transcript_rsem.csv")

# alpha_omega_counts <- read_csv("ome-data/alphaomega_transcript_rsem.csv")


alpha_omega_counts <- extract_rsem_isoform_counts("data_new/gene_counts/Alpha_Omega_gene_counts/")

# load the alpha omega metadata 
alpha_omega_metadata <- readRDS("data_new/processed_data/Alpha_Omega_metadata.RDS")

alpha_omega <- alpha_omega_counts %>% 
  separate(transcript_id, c("transcript_id", "transcript_name"), sep = "_", extra = "merge") %>%
  
  # drop transcript_id, select transcripte_name and any of the sample names that match sample name in metadata
  dplyr::select(transcript_id, all_of(alpha_omega_metadata$seq_sample_id))


## Keep nonzero rows
nonzero <- alpha_omega %>%
  dplyr::filter(rowSums(alpha_omega[, -1]) != 0)

# saveRDS(nonzero, file = "data_new/gene_counts/nonzero_alphaomega_counts.RDS")
# Keep filtered genes (based on time)
# time filtering assures that genes are kept that are expressed
# in all times as that is the main comparing factor in subsequent 
# models.

filtered <- nonzero %>%
  dplyr::filter(filterByExpr(nonzero[,-1], # exclude the first column
                             min.count = 1,
                          #   min.total.count = 15,
                             large.n = 10, 
                             min.prop = 0.7,
                             group = paste(alpha_omega_metadata$time)))  



# saveRDS(filtered, "data_new/gene_counts/filtered_by_time_alphaomega.RDS")



# Create dge lists and calculate norm factors
dge_ao   <- DGEList(filtered[,-1])


dge_ao  <- calcNormFactors(dge_ao, method = "TMM")


# get the normalised counts

normalised_ao_counts <- as.data.frame(cpm(dge_ao, normalized.lib.sizes = T))

normalised_ao_counts$transcript_ID <- filtered$transcript_id

# Make the transcript_ID the first column

normalised_ao_counts <- normalised_ao_counts %>%
  dplyr::select(transcript_ID, everything())


saveRDS(normalised_ao_counts, "data_new/gene_counts/cpm_normalised_AO_counts.RDS")


## Add effective library size to 
# alpha_omega_metadata <- dge_ct$samples %>%
#   rownames_to_column(var = "seq_sample_id") %>%
#    dplyr::select(-group) %>%
#   inner_join(alpha_omega_metadata, by = "seq_sample_id") %>%
#   mutate(efflibsize = (lib.size * norm.factors) / median(lib.size * norm.factors))
# 




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


## Keep nonzero rows
nonzero_vol <- volume %>%
  dplyr::filter(rowSums(volume[,-1]) != 0)


# Keep filtered genes (based on time)
# time filtering assures that genes are kept that are expressed
# in all times as that is the main comparing factor in subsequent 
# models.

filtered_vol <- nonzero_vol %>%
  dplyr::filter(filterByExpr(nonzero_vol[,-1], 
                             min.count = 1,
                             #   min.total.count = 15,
                             large.n = 10, 
                             min.prop = 0.7,
                             group = paste(volume_metadata$time)))  






# Create dge lists and calculate norm factors
dge_vol   <- DGEList(filtered_vol[,-1])


dge_vol  <- calcNormFactors(dge_vol, method = "TMM")


# get the normalised counts

normalised_vol_counts <- as.data.frame(cpm(dge_vol, normalized.lib.sizes = T))

normalised_vol_counts$transcript_ID <- filtered_vol$transcript_id

# Make the transcript_ID the first column

normalised_vol_counts <- normalised_vol_counts %>%
  dplyr::select(transcript_ID, everything())


saveRDS(normalised_vol_counts, "data_new/gene_counts/cpm_normalised_vol_counts.RDS")






# Contratrain data

Ct_counts <- extract_rsem_isoform_counts("data_new/gene_counts/Contratrain_RSEM_outputs_new/")


# load the volume metadata 
ct_metadata <- readRDS("data_new/processed_data/contratrain_metadata.RDS")



Contratrain <- Ct_counts %>% 
  separate(transcript_id, c("transcript_id", "transcript_name"), sep = "_", extra = "merge") %>%
  
  # drop transcript_id, select transcripte_name and any of the sample names that match sample name in metadata
  dplyr::select(transcript_id, all_of(ct_metadata$seq_sample_id))


## Keep nonzero rows
nonzero_ct <- Contratrain %>%
  dplyr::filter(rowSums(Contratrain[,-1]) != 0)


# Keep filtered genes (based on time)
# time filtering assures that genes are kept that are expressed
# in all times as that is the main comparing factor in subsequent 
# models.

filtered_ct <- nonzero_ct %>%
  dplyr::filter(filterByExpr(nonzero_ct[,-1], 
                             min.count = 1,
                             #   min.total.count = 15,
                             large.n = 10, 
                             min.prop = 0.7,
                             group = paste(ct_metadata$time)))  






# Create dge lists and calculate norm factors
dge_ct   <- DGEList(filtered_ct[,-1])


dge_ct  <- calcNormFactors(dge_ct, method = "TMM")


# get the normalised counts

normalised_ct_counts <- as.data.frame(cpm(dge_ct, normalized.lib.sizes = T))

normalised_ct_counts$transcript_ID <- filtered_ct$transcript_id

# Make the transcript_ID the first column

normalised_ct_counts <- normalised_ct_counts %>%
  dplyr::select(transcript_ID, everything())


saveRDS(normalised_ct_counts, "data_new/gene_counts/cpm_normalised_contratrain_counts.RDS")







Copd_counts <- extract_rsem_isoform_counts("data_new/gene_counts/COPD_RSEM_outputs_new/")

# load COPD metadata
copd_metadata <- readRDS("data_new/processed_data/copd_metadata.RDS")








Relief_counts <- extract_rsem_isoform_counts("data_new/gene_counts/Relief_RSEM_outputs/")

SRP102542_counts <- extract_rsem_isoform_counts("data_new/gene_counts/SRP102542_RSEM_outputs/")








