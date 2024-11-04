# Explore how RT generally and its conditions affect the DS introns
library(dplyr)
library(tidyverse)
library(seqwrap)
library(gridExtra)
library(ggpubr)
library(cowplot)


# Load the Trainome functions
# Contains functions needed for model building
source("R/Trainome_functions.R")


#COPD metadata
copd_metadata <- readRDS("data/processed_data/copd_metadata.RDS")


# hist(copd_metadata$age)
# range(copd_metadata$age)
 unique(copd_metadata$time)
 
length(unique(copd_metadata$participant)) 
length(unique(copd_metadata$seq_sample_id))
# Volume_data
Vol_metadata <- readRDS("data/processed_data/volume_metadata.RDS")

# unique(Vol_metadata$time)
length(unique(Vol_metadata$participant))
length(Vol_metadata$seq_sample_id)
range(Vol_metadata$age)


ct_metadata <- readRDS("data/processed_data/contratrain_metadata.RDS")
