library(dplyr)
library(tidyverse)
library(seqwrap)
library(effects)
library(scales)
library(glmmTMB)


# Load the Trainome functions
# Contains functions needed for model building
source("R/Trainome_functions.R")



all_splice_df <- readRDS("data/all_splice.RDS") %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  drop_na() 



all_full_metadata <- readRDS("data_new/processed_data/all_full_metadata.RDS")
colnames(all_full_metadata)


# REORDER THE SEQUENCE ID TO MATCH BOTH DATAFRAMMES
all_splice_reordered <- all_splice_df[, c("transcript_ID",all_full_metadata$seq_sample_id)]

# Check if everything matches except the transcript_id
match(colnames(all_splice_reordered), all_full_metadata$seq_sample_id)


# convert the 1.0 to 0.999. This is becasue beta-model accepts only values between 0 and one
all_splice_reordered[all_splice_reordered == 1 ] <- 0.999



# initialise the argument
args_full <-list(formula = y ~  scaled_age * time + sex + (1|study) +(1|participant), 
                 family = glmmTMB::beta_family())


# check the functions and datasets
container <- seqwrap_compose(data = all_splice_reordered,
                             metadata = all_full_metadata,
                             samplename = "seq_sample_id",
                             modelfun = glmmTMB::glmmTMB,
                             arguments = args_full,
                             summary_fun = sum_with_pred,
                             eval_fun = eval_mod )


# build model
model <- seqwrap(container,
                 summary_fun = sum_with_pred,
                 eval_fun = eval_mod,
                 return_models = F,
                 # subset = 1:15,
                 cores = ncores-2)


model_summary <- seqwrap_summarise(model)
head(model_summary)
str(model)
saveRDS(model_summary, "data/RT_model_summary.RDS")
