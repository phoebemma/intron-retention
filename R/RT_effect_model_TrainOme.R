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
copd_metadata <- readRDS("data_new/processed_data/copd_metadata.RDS") %>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)

colnames(copd_metadata)
# hist(copd_metadata$age)
# range(copd_metadata$age)
 unique(copd_metadata$time)
 unique(copd_metadata$volume)
 
length(unique(copd_metadata$participant)) 
length(unique(copd_metadata$seq_sample_id))
# Volume_data
Vol_metadata <- readRDS("data_new/processed_data/volume_metadata.RDS")%>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)
colnames(Vol_metadata)

unique(Vol_metadata$time)
# unique(Vol_metadata# unique(Vol_metadata$time)
length(unique(Vol_metadata$participant))
length(Vol_metadata$seq_sample_id)
range(Vol_metadata$age)


ct_metadata <- readRDS("data_new/processed_data/contratrain_metadata.RDS") %>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)
colnames(ct_metadata)

unique(ct_metadata$time)

SRP102542_metadata <- readRDS("data_new/processed_data/SRP102542_metadata.RDS")%>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)
unique(SRP102542_metadata$time)



# Alpha and Omega data

A_Omega_metadata <- readRDS("data_new/processed_data/Alpha_Omega_metadata.RDS")%>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)







all_full_metadata <- rbind(copd_metadata, Vol_metadata)%>%
  rbind(ct_metadata) %>%
  rbind(A_Omega_metadata) %>%
  rbind(SRP102542_metadata) %>%
  
  mutate(across(c("age"), round, 0)) %>%
  
  mutate(group = case_when(age <=20 ~ "<=20" ,
                           age > 20 & age <= 30 ~ ">20 & <=30",
                           age > 30 & age <= 40 ~ ">30 & <=40", 
                           age > 40 & age <= 50 ~ ">40 & <=50",
                           age > 50 & age <= 60 ~ ">50 & <=60",
                           age > 60 & age <= 70 ~ ">60 & <=70",
                           age > 70 ~ ">70")) %>%
  mutate(sex = factor(sex, levels = c("female", "male")),
         time = factor(time, levels = c("PreExc", "PostExc")),
         group = factor(group, levels = c("<=20" ,">20 & <=30", ">30 & <=40", ">40 & <=50",
                                          ">50 & <=60",">60 & <=70", ">70"))) 

unique(all_full_metadata$sex)
unique(all_full_metadata$study)
unique(all_full_metadata$time)
unique(all_full_metadata$group)


copd_splice_df <- readRDS("data_new/processed_data/copd_splicing_data.RDS")

vol_splice_df <- readRDS("data_new/processed_data/volume_splicing_data.RDS")
ct_splice_df <- readRDS("data_new/processed_data/contratrain_splicing_data.RDS")
SRP102542_splice_df <- readRDS("data_new/processed_data/SRP102542_splicing_data.RDS")
AOD_splice_df <- readRDS("data_new/processed_data/Alpha_Omega_splicing_data.RDS")

all_splice_df <-copd_splice_df  %>%
  inner_join(vol_splice_df, by = "transcript_ID") %>%
  inner_join(ct_splice_df, by = "transcript_ID") %>%
  inner_join(AOD_splice_df, by = "transcript_ID") %>%
  inner_join(SRP102542_splice_df, by = "transcript_ID") 


# select only the splicing samples captured in the metadata
all_intersect <- intersect(colnames(all_splice_df), all_full_metadata$seq_sample_id)

all_splice_df <- all_splice_df %>%
  subset(select = c("transcript_ID", all_intersect))%>%
  drop_na()

# Get the ds introns by group and sex without interaction
# That proved to be the best performing model
introns_of_int <- readRDS("data_new/models/filt_preExc_group_sex_no_int.RDS")%>%
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

args<- list(formula = y ~  time + (1|study) +(1|participant), 
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


# missing <- names(which(RT_impact_model$summaries == "NULL"))
# avail <- names(which(RT_impact_model$summaries != "NULL"))

#Remove all that have output NULL



mod_sum <- model_sum(RT_impact_model, 2)



mod_eval <- model_eval(RT_impact_model)


model_cont <- mod_sum %>%
  inner_join(mod_eval, by = "target")  %>%
  filter(adj.p <= 0.05)




# Evaluate the impact of RT on the full dataset


all_splice_reordered <- all_splice_df[ , c("transcript_ID",all_full_metadata$seq_sample_id)]

colnames(all_pre_splice_reordered)
rownames(all_pre_metadata)

# Check if everything matches except the transcript_id
match(colnames(all_splice_reordered), all_full_metadata$seq_sample_id)

# convert the 1.0 to 0.999. This is becasue beta-model accepts only values between 0 and one
all_splice_reordered[all_splice_reordered == 1 ] <- 0.999

args<- list(formula = y ~  time + (1|study) +(1|participant), 
            family = glmmTMB::beta_family())



RT_model <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                           arguments = args,
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

RT_model$summaries


missing <- names(which(RT_model$summaries == "NULL"))
avail <- names(which(RT_model$summaries != "NULL"))



mod_sum <- bind_rows(within(RT_model$summaries, rm(missing))) %>%
  mutate(target = rep(avail, each = 2)) %>%
  subset(coef != "(Intercept)")  %>%
   mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
 log2fc = Estimate/log(2),
 fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) 



mod_eval <- bind_rows(within(RT_model$evaluations, rm(missing)))%>%
  mutate(target = avail)


model_cont <- mod_sum %>%
  inner_join(mod_eval, by = "target")# %>%
#  filter(adj.p <= 0.05)

hist(model_cont$Estimate)

saveRDS(model_cont, "data_new/models/full_data_RT_model.RDS")






# Since only few of the ds by age introns were contained in the full data, might be nice to look at the 
# Impact of age and exercise in one go


args_full <- args<- list(formula = y ~  group*time + sex + (1|study) +(1|participant), 
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
  mutate(target = rep(avail_full, each = 15)) %>%
  subset(coef != "(Intercept)")  %>%
   mutate(.by = coef,
          adj.p = p.adjust(Pr...z.., method = "fdr"),
          log2fc = Estimate/log(2),
          fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) 



mod_eval <- bind_rows(within(full_RT_model$evaluations, rm(missing_full)))%>%
  mutate(target = avail_full)


model_cont <- mod_sum %>%
  inner_join(mod_eval, by = "target")# %>%
  #filter(adj.p <= 0.05)

length(unique(model_cont$target))
hist(model_cont$Estimate)

saveRDS(model_cont, "data_new/models/full_data_RT_and_age_model.RDS")
