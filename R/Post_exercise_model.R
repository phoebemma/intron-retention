# This is the updated post-excercise model
# Here, the models will be built after removing the data that score perfectly across all samples

library(dplyr)
library(tidyverse)
library(seqwrap)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(effects)
library(scales)
library(glmmTMB)


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


saveRDS(High_SE_df, "data_new/processed_data/perfectly_spliced_across_all.RDS")

non_perfect_SE <- all_splice_df %>%
  filter(!(transcript_ID %in% High_SE_df$transcript_ID))


all_full_metadata <- readRDS("data_new/processed_data/all_full_metadata.RDS")

all_full_metadata$scaled_age <- round(rescale(all_full_metadata$age), digits = 2)




# Select only the splicing data whose metadata is available
all_full_metadata <- all_full_metadata %>%
  filter((seq_sample_id %in% colnames(non_perfect_SE [,-1]))) 

saveRDS(all_full_metadata, "data_new/processed_data/all_full_metadata.RDS")
# REORDER THE SEQUENCE ID TO MATCH BOTH DATAFRAMMES
all_splice_reordered <- non_perfect_SE[, c("transcript_ID",all_full_metadata$seq_sample_id)]


# Check if everything matches except the transcript_id
match(colnames(all_splice_reordered), all_full_metadata$seq_sample_id)

# convert the 1.0 to 0.999. This is becasue beta-model accepts only values between 0 and one
all_splice_reordered[all_splice_reordered == 1 ] <- 0.999




args_full <-list(formula = y ~  scaled_age * time + sex + (1|study) +(1|participant), 
                 family = glmmTMB::beta_family())



full_RT_model <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                         arguments = args_full,
                         data = all_splice_reordered,
                         metadata = all_full_metadata,
                         samplename = "seq_sample_id",
                         summary_fun = sum_with_pred,
                         eval_fun = eval_mod,
                         exported = list(),
                         save_models = F,
                         return_models = F,
                          #subset = 1:10,
                         cores = ncores-2)

ex <- full_RT_model$summaries$ENST00000001008.6_5_12 %>%
  select(scaled_age, time, fit)
full_RT_model$summaries



saveRDS(full_RT_model, "data_new/simpler_full_RT_model.RDS")
#saveRDS(full_RT_model, "data_new/full_RT_model.RDS")



missing_full <- names(which(full_RT_model$summaries == "NULL"))
avail_full <- names(which(full_RT_model$summaries != "NULL"))



mod_sum <- bind_rows(within(full_RT_model$summaries, rm(missing_full))) %>%
  mutate(target = rep(avail_full, each = 27)) %>%
 # subset(coef != "(Intercept)")  %>%
  mutate(adj.p = p.adjust( Pr...z.., method = "fdr"),
         log2fc = Estimate/log(2),
         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns"),
         odds_ratio = exp(Estimate),
         type = "model coefficient")


saveRDS(mod_sum, "data_new/simpler_model_summary.RDS")
#saveRDS(mod_sum, "data_new/model_summary_RT_model.RDS")
model_df <- baseline_predictions %>%
  inner_join(mod_sum, by = "target")
