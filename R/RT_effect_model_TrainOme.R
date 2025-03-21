# Explore how RT generally and its conditions affect the DS introns
library(dplyr)
library(tidyverse)
library(seqwrap)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(scales)


# Load the Trainome functions
# Contains functions needed for model building
source("R/Trainome_functions.R")


# Load the full metadata 
all_full_metadata <- readRDS("data_new/processed_data/all_full_metadata.RDS")%>%
  mutate(group = case_when(age <=20 ~ "20 and below" ,
                           age > 20 & age < 30 ~ "21 to 29",
                           age >= 30 & age < 40 ~ "30 to 39", 
                           age >= 40 & age < 50 ~ "40 to 49",
                           age >= 50 & age < 60 ~ "50 to 59",
                           age >= 60 & age < 70 ~ "60 to 69",
                           age >= 70 & age < 80 ~ "70 to 79",
                           age >= 80 ~ "80 and above")) %>%
  mutate(group = factor(group, levels = c("20 and below", "21 to 29", "30 to 39",
                                          "40 to 49", "50 to 59",  "60 to 69",
                                          "70 to 79", "80 and above" )))

length(unique(all_pre_metadata$participant))
# Standardize the age by scaling them 0 to 1

all_full_metadata$scaled_age <- round(rescale(all_full_metadata$age), digits = 2)
# Load full splice data
all_splice_df <- readRDS("data_new/processed_data/all_splice_data.RDS")



# Evaluate the impact of RT on the full dataset


# Select only the splicing data whose metadata is available
all_full_metadata <- all_full_metadata %>%
  filter((seq_sample_id %in% colnames(all_splice_df[,-1]))) 


# REORDER THE SEQUENCE ID TO MATCH BOTH DATAFRAMMES
all_splice_reordered <- all_splice_df[, c("transcript_ID",all_full_metadata$seq_sample_id)]


# Check if everything matches except the transcript_id
match(colnames(all_splice_reordered), all_full_metadata$seq_sample_id)

# convert the 1.0 to 0.999. This is becasue beta-model accepts only values between 0 and one
all_splice_reordered[all_splice_reordered == 1 ] <- 0.999


# Since only few of the ds by age introns were contained in the full data, might be nice to look at the 
# Impact of age and exercise in one go


args_full <-list(formula = y ~  scaled_age*time + (1|study) + (scaled_age+0|study) +(1|participant), 
                        family = glmmTMB::beta_family())



full_RT_model <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                    arguments = args_full,
                    data = all_splice_reordered,
                    metadata = all_full_metadata,
                    samplename = "seq_sample_id",
                    summary_fun = sum_fun,
                    eval_fun = eval_mod,
                    exported = list(),
                    save_models = FALSE,
                    return_models = FALSE,
                    # subset = 1:10,
                    cores = ncores-2)


full_RT_model$summaries


missing_full <- names(which(full_RT_model$summaries == "NULL"))
avail_full <- names(which(full_RT_model$summaries != "NULL"))



mod_sum <- bind_rows(within(full_RT_model$summaries, rm(missing_full))) %>%
  mutate(target = rep(avail_full, each = 4)) %>%
  subset(coef != "(Intercept)")  %>%
   mutate(.by = coef,
          adj.p = p.adjust(Pr...z.., method = "fdr"),
          log2fc = Estimate/log(2),
          fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) 



mod_eval <- bind_rows(within(full_RT_model$evaluations, rm(missing_full)))%>%
  mutate(target = avail_full)



model_cont <- mod_sum %>%
  inner_join(mod_eval, by = "target") # %>%
# filter(Pr...z..<= 0.05)

# visualise
model_cont %>%
  ggplot(aes(x = coef)) +
  geom_bar()+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))


saveRDS(model_cont, "data_new/models/full_data_RT_scaled_age_int_model.RDS")




# Check how age group affects this

args_group <-list(formula = y ~  group*time + (1|study) + (1|participant), 
                 family = glmmTMB::beta_family())



group_RT_model <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                         arguments = args_group,
                         data = all_splice_reordered,
                         metadata = all_full_metadata,
                         samplename = "seq_sample_id",
                         summary_fun = sum_fun,
                         eval_fun = eval_mod,
                         exported = list(),
                         save_models = FALSE,
                         return_models = FALSE,
                         # subset = 1:10,
                         cores = ncores-2)


group_RT_model$summaries


missing_group <- names(which(group_RT_model$summaries == "NULL"))
avail_group <- names(which(group_RT_model$summaries != "NULL"))



mod_sum_group <- bind_rows(within(group_RT_model$summaries, rm(missing_group))) %>%
  mutate(target = rep(avail_group, each = 16)) %>%
  subset(coef != "(Intercept)")  %>%
  mutate(.by = coef,
         adj.p = p.adjust(Pr...z.., method = "fdr"),
         log2fc = Estimate/log(2),
         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) 



mod_eval_group <- bind_rows(within(group_RT_model$evaluations, rm(missing_group)))%>%
  mutate(target = avail_group)



model_cont_group <- mod_sum_group %>%
  inner_join(mod_eval_group, by = "target") 

saveRDS(model_cont_group, "data_new/models/agegroup_RT_model.RDS")