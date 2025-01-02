# Explore the ds introns on Alpha-Omega data
library(dplyr)
library(tidyverse)
library(seqwrap)

source("R/Trainome_functions.R")

A_Omega_metadata <- readRDS("data_new/processed_data/Alpha_Omega_metadata.RDS") %>%
  mutate(sex = factor(sex, levels = c("female", "male")),
         time = factor(time, levels = c("PreExc", "PostExc")))
unique(A_Omega_metadata$time)
AOD_splice_df <- readRDS("data_new/processed_data/Alpha_Omega_splicing_data.RDS")

AOD_intersect <- intersect(colnames(AOD_splice_df), A_Omega_metadata$seq_sample_id)

AOD_splice_df <- AOD_splice_df%>%
  subset(select = c("transcript_ID", AOD_intersect))%>%
  drop_na()

# Get the ds introns by group and sex without interaction
# That proved to be the best performing model
introns_of_int <- readRDS("data_new/models/filt_preExc_group_sex_no_int.RDS")%>%
  pull(target)

full_splice_df <- AOD_splice_df [AOD_splice_df $transcript_ID %in% introns_of_int,]

# Check if everything matches except the transcript_id
match(colnames(full_splice_df), A_Omega_metadata$seq_sample_id)

full_splice_reordered <- full_splice_df[ , c("transcript_ID",A_Omega_metadata$seq_sample_id)]

# Check if everything matches except the transcript_id
match(colnames(full_splice_reordered), A_Omega_metadata$seq_sample_id)


full_splice_reordered[full_splice_reordered == 1 ] <- 0.999

# model  age_group, time, volume and its interaction with condition

args<- list(formula = y ~  time  +(1|participant), 
            family = glmmTMB::beta_family())



RT_AO_model <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                           arguments = args,
                           data = full_splice_reordered,
                           metadata = A_Omega_metadata,
                           samplename = "seq_sample_id",
                           summary_fun = sum_fun,
                           eval_fun = eval_mod,
                           exported = list(),
                           save_models = FALSE,
                           return_models = FALSE,
                           # subset = 1:10,
                           cores = ncores-2)
RT_AO_model$summaries


mod_sum <- model_sum(RT_AO_model, 2)



mod_eval <- model_eval(RT_AO_model)


model_cont <- mod_sum %>%
  inner_join(mod_eval, by = "target")  %>%
  filter(adj.p <= 0.05)
