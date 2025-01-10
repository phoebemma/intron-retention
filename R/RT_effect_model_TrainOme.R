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


# Load the full metadata 
all_full_metadata <- readRDS("data_new/processed_data/all_full_metadata.RDS")

# Load full splice data
all_splice_df <- readRDS("data_new/processed_data/all_splice_data.RDS")

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

colnames(all_splice_reordered)
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
  inner_join(mod_eval, by = "target") # %>%
  #filter(adj.p <= 0.05)

hist(model_cont$Estimate)
length(unique(model_cont$target))
saveRDS(model_cont, "data_new/models/full_data_RT_model.RDS")






# Since only few of the ds by age introns were contained in the full data, might be nice to look at the 
# Impact of age and exercise in one go


args_full <-list(formula = y ~  group*time + sex + (1|study) +(1|participant), 
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
  inner_join(mod_eval, by = "target") # %>%
#  filter(adj.p <= 0.05)

length(unique(model_cont$target))
hist(model_cont$Estimate)

model_cont %>%
  ggplot(aes(x = coef)) +
  geom_bar()+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

unique(model_cont$coef)
saveRDS(model_cont, "data_new/models/full_data_RT_and_age_model.RDS")




# Check the interaction of sex and time in addition to group and time
args_full_2 <-list(formula = y ~  group*time + sex*time + (1|study) +(1|participant), 
                 family = glmmTMB::beta_family())



full_RT_model2 <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                         arguments = args_full_2,
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

full_RT_model2$summaries


missing_full2 <- names(which(full_RT_model2$summaries == "NULL"))
avail_full2 <- names(which(full_RT_model2$summaries != "NULL"))



mod_sum2 <- bind_rows(within(full_RT_model2$summaries, rm(missing_full2))) %>%
  mutate(target = rep(avail_full2, each = 16)) %>%
  subset(coef != "(Intercept)")  %>%
  mutate(.by = coef,
         adj.p = p.adjust(Pr...z.., method = "fdr"),
         log2fc = Estimate/log(2),
         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) 



mod_eval2 <- bind_rows(within(full_RT_model2$evaluations, rm(missing_full2)))%>%
  mutate(target = avail_full2)


model_cont2 <- mod_sum2 %>%
  inner_join(mod_eval2, by = "target") # %>%
#  filter(adj.p <= 0.05)

length(unique(model_cont$target))
hist(model_cont$Estimate)

saveRDS(model_cont2, "data_new/models/full_data_RT_interction_group_sex_model.RDS")
