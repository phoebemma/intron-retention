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

colnames(copd_metadata)
# hist(copd_metadata$age)
# range(copd_metadata$age)
 unique(copd_metadata$time)
 
length(unique(copd_metadata$participant)) 
length(unique(copd_metadata$seq_sample_id))
# Volume_data
Vol_metadata <- readRDS("data/processed_data/volume_metadata.RDS")
colnames(Vol_metadata)

# unique(Vol_metadata# unique(Vol_metadata$time)
length(unique(Vol_metadata$participant))
length(Vol_metadata$seq_sample_id)
range(Vol_metadata$age)


ct_metadata <- readRDS("data/processed_data/contratrain_metadata.RDS")
colnames(ct_metadata)

all_full_metadata <- rbind(copd_metadata, Vol_metadata)%>%
  rbind(ct_metadata) %>%
  
  mutate(age_group = case_when(study == "copd" ~ "Old", 
                               study == "vol" ~ "Young",
                               study == "ct" ~ "Young")) %>%
  mutate(across(c("age"), round, 0)) %>%
  mutate(group = case_when(age <=25 ~ "<=25" ,
                           age > 25 & age <= 50 ~ ">25 & <=50", 
                           age > 50 & age <= 70 ~ ">50 & <=70",
                           age > 70 ~ ">70")) %>%
  mutate(sex = factor(sex, levels = c("female", "male")),
         time = factor(time, levels = c("PreExc", "PostExc")),
         condition = factor(condition, levels = c("RM10", "RM30")),
         age_group = factor(age_group, levels= c("Young", "Old")),
         group = factor(group, levels = c("<=25", ">25 & <=50", ">50 & <=70", ">70"))) 

unique(all_full_metadata$sex)
unique(all_full_metadata$study)
unique(all_full_metadata$group)
unique(all_full_metadata$age_group)
unique(all_full_metadata$time)
unique(all_full_metadata$volume)
unique(all_full_metadata$condition)

copd_splice_df <- readRDS("data/processed_data/copd_splicing_data.RDS")

vol_splice_df <- readRDS("data/processed_data/volume_splicing_data.RDS")
ct_splice_df <- readRDS("data/processed_data/contratrain_splicing_data.RDS")


all_splice_df <- copd_splice_df%>%
  inner_join(vol_splice_df, by = "transcript_ID") %>%
  inner_join(ct_splice_df, by = "transcript_ID") %>%
  drop_na()

# Get the ds introns by age
introns_of_int <- readRDS("data/Trainome_data_models/preExc_model_age_group.RDS")%>%
  pull(target)


full_splice_df <- all_splice_df[all_splice_df$transcript_ID %in% introns_of_int,]

# Check if everything matches except the transcript_id
match(colnames(full_splice_df), all_full_metadata$seq_sample_id)

full_splice_reordered <- full_splice_df[ , c("transcript_ID",all_full_metadata$seq_sample_id)]

# Check if everything matches except the transcript_id
match(colnames(full_splice_reordered), all_full_metadata$seq_sample_id)



# convert the 1.0 to 0.999. This is becasue beta-model accepts only values between 0 and one
full_splice_reordered[full_splice_reordered == 1 ] <- 0.999

# model  age_group, time, volume and its interaction with condition

args<- list(formula = y ~  age_group*time + condition + volume + (1|study) +(1|participant), 
            family = glmmTMB::beta_family())



RT_impact_model <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                        arguments = args,
                        data = full_splice_reordered,
                        metadata = all_full_metadata,
                        samplename = "seq_sample_id",
                        summary_fun = sum_fun,
                        eval_fun = eval_mod,
                        exported = list(),
                        save_models = FALSE,
                        return_models = FALSE,
                        # subset = 1:10,
                        cores = ncores-2)
RT_impact_model$summaries


args_2 <- list(formula = y ~  time + condition + volume + (1|study) +(1|participant), 
               family = glmmTMB::beta_family())


RT_impact_model_2 <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                           arguments = args_2,
                           data = full_splice_reordered,
                           metadata = all_full_metadata,
                           samplename = "seq_sample_id",
                           summary_fun = sum_fun,
                           eval_fun = eval_mod,
                           exported = list(),
                           save_models = FALSE,
                           return_models = FALSE,
                           # subset = 1:10,
                           cores = ncores-2)

