# This is the updated post-excercise model
# Here, the models will be built after removing the data that score perfectly across all samples

library(dplyr)
library(tidyverse)
library(seqwrap)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(effects)


# Load the Trainome functions
# Contains functions needed for model building
source("R/Trainome_functions.R")



all_splice_df <- readRDS("data_new/processed_data/all_splice_data.RDS")


# Extract the introns that are perfectly spliced across all samples

# Get the perfectly spliced introns across all datasets
High_SE_df <- all_splice_df %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID)) %>%
  summarise(.by = transcript_ID, 
            mode = getmode(SE),
            min = min(SE), 
            max = max(SE)) %>%
  filter(min == 1) 


non_perfect_SE <- all_splice_df %>%
  filter(!(transcript_ID %in% High_SE_df$transcript_ID))


all_full_metadata <- readRDS("data_new/processed_data/all_full_metadata.RDS")

all_full_metadata$scaled_age <- round(rescale(all_full_metadata$age), digits = 2)




# Select only the splicing data whose metadata is available
all_full_metadata <- all_full_metadata %>%
  filter((seq_sample_id %in% colnames(non_perfect_SE [,-1]))) 


# REORDER THE SEQUENCE ID TO MATCH BOTH DATAFRAMMES
all_splice_reordered <- non_perfect_SE[, c("transcript_ID",all_full_metadata$seq_sample_id)]


# Check if everything matches except the transcript_id
match(colnames(all_splice_reordered), all_full_metadata$seq_sample_id)

# convert the 1.0 to 0.999. This is becasue beta-model accepts only values between 0 and one
all_splice_reordered[all_splice_reordered == 1 ] <- 0.999


# Since only few of the ds by age introns were contained in the full data, might be nice to look at the 
# Impact of age and exercise in one go


args_full <-list(formula = y ~  scaled_age*time + sex*time + (1|study) + (scaled_age+0|study) +(1|participant), 
                 family = glmmTMB::beta_family())



full_RT_model <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                         arguments = args_full,
                         data = all_splice_reordered,
                         metadata = all_full_metadata,
                         samplename = "seq_sample_id",
                         summary_fun = sum_fun2,
                         eval_fun = eval_mod,
                         exported = list(),
                         save_models = FALSE,
                         return_models = FALSE,
                         # subset = 1:100,
                         cores = ncores-2)
