source("R/Trainome_functions.R")
library(dplyr)
library(tidyverse)
library(ggplot2)
library(trainomeMetaData)

Vol_counts <- extract_rsem_gene_counts("data_new/gene_counts/Volume_RSEM_outputs_new/") 

Ct-counts <- extract_rsem_gene_counts("data_new/gene_counts/Contratrain_RSEM_outputs_new/")

Copd_counts <- extract_rsem_gene_counts("data_new/gene_counts/COPD_RSEM_outputs_new/")

# download_ome(download= "alphaomega_gene_rsem")

# Alpha_omega <- read_csv("ome-data/alphaomega_gene_rsem.csv")
