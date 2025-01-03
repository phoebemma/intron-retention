# Load the file with libraries
# source("R/libraries.R")
library(dplyr)
library(tidyverse)
library(seqwrap)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(trainomeMetaData)


# Load the Trainome functions
# Contains functions needed for model building
source("R/Trainome_functions.R")

# Load metadata
all_pre_metadata <- readRDS("data_new/Pre_Exercise/all_prexercise_metadata.RDS")




# splicing data
all_pre_splice <- readRDS("data_new/Pre_Exercise/all_pre_Exc_splicing_data.RDS")

# Load Splicing data

# reorder the column name to match how they occur in the metadata

# Filter to remove missing values
all_pre_metadata <- all_pre_metadata %>%
  filter((seq_sample_id %in% colnames(all_pre_splice[,-1]))) %>%
  print()


all_pre_splice_reordered <- all_pre_splice[ , c("transcript_ID",all_pre_metadata$seq_sample_id)]

colnames(all_pre_splice_reordered)
rownames(all_pre_metadata)

# Check if everything matches except the transcript_id
 match(colnames(all_pre_splice_reordered), all_pre_metadata$seq_sample_id)


# invert the data to create a dataframe suitable for zero inflated analyses
# splice_df_inverted <- all_pre_splice_reordered %>%
#   mutate(across(X102PreExcVLR12:X134.subj8sample4, function(x)1-x))



# convert the 1.0 to 0.999. This is becasue beta-model accepts only values between 0 and one
all_pre_splice_reordered[all_pre_splice_reordered == 1 ] <- 0.999

# This argument models for age as a continous variable

args<- list(formula = y ~  age*sex + (1|study) +(1|participant), 
            family = glmmTMB::beta_family())

 

SE_model <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                    arguments = args,
                    data = all_pre_splice_cont,
                    metadata = all_pre_metadata,
                    samplename = "seq_sample_id",
                    summary_fun = sum_fun,
                    eval_fun = eval_mod,
                    exported = list(),
                    save_models = FALSE,
                    return_models = FALSE,
                   # subset = 1:10,
                    cores = ncores-2)

SE_model$summaries

SE_model$summaries$ENST00000296098.4_6_2


excl_1 <- names(which(SE_model$summaries == "NULL"))
geneids_1 <- names(which(SE_model$summaries != "NULL"))

#Remove all that have output NULL



mod_sum <- bind_rows(within(SE_model$summaries, rm(excl_1))) %>%
  mutate(target = rep(geneids_1, each = 4)) %>%
  subset(coef != "(Intercept)") # %>%
 #   mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
 # log2fc = Estimate/log(2),
 # fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) 



mod_eval <- bind_rows(within(SE_model$evaluations, rm(excl_1)))%>%
  mutate(target = geneids_1)


model_cont <- mod_sum %>%
  inner_join(mod_eval, by = "target")

hist(model_cont$Pr...z..)

# saveRDS(model_cont, "data/re_models/primary_preExc_count_interaction_model.RDS")




# Model using age_group


args<- list(formula = y ~  group*sex + (1|study) +(1|participant), 
            family = glmmTMB::beta_family())

SE_count_model <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                    arguments = args,
                    data = all_pre_splice_reordered,
                    metadata = all_pre_metadata,
                    samplename = "seq_sample_id",
                    summary_fun = sum_fun,
                    eval_fun = eval_mod,
                    exported = list(),
                    save_models = FALSE,
                    return_models = FALSE,
                    cores = ncores-2)

SE_count_model$summaries$ENST00000007516.8_2_16

excl <- names(which(SE_count_model$summaries == "NULL"))
geneids <- names(which(SE_count_model$summaries != "NULL"))


mod_sum_group <- bind_rows(within(SE_count_model$summaries, rm(excl))) %>%
  mutate(target = rep(geneids, each = 14))  %>%
  subset(coef != "(Intercept)") # %>%
  #  mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
  #         log2fc = Estimate/log(2),
  #         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) %>%
  # filter(adj.p <= 0.05)



mod_eval_group <- bind_rows(within(SE_count_model$evaluations, rm(excl)))%>%
  mutate(target = geneids)


model_cont_group <- mod_sum_group %>%
  inner_join(mod_eval_group, by = "target")

# This contains the interaction between age-group and sex
saveRDS(model_cont_group, "data_new/models/preExc_group_and_sex_model.RDS")




# Model for age group alone. Without sex interaction

args<- list(formula = y ~  group + (1|study) +(1|participant), 
            family = glmmTMB::beta_family())

SE_group_model <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                          arguments = args,
                          data = all_pre_splice_reordered,
                          metadata = all_pre_metadata,
                          samplename = "seq_sample_id",
                          summary_fun = sum_fun,
                          eval_fun = eval_mod,
                          exported = list(),
                          save_models = FALSE,
                          return_models = FALSE,
                          cores = ncores-2)

SE_group_model$summaries$ENST00000007516.8_2_16

excl_2 <- names(which(SE_group_model$summaries == "NULL"))
geneids_2 <- names(which(SE_group_model$summaries != "NULL"))


mod_sum_group <- bind_rows(within(SE_group_model$summaries, rm(excl_2))) %>%
  mutate(target = rep(geneids_2, each = 7))  %>%
  subset(coef != "(Intercept)") # %>%
#  mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
#         log2fc = Estimate/log(2),
#         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) %>%
# filter(adj.p <= 0.05)



mod_eval_group <- bind_rows(within(SE_group_model$evaluations, rm(excl_2)))%>%
  mutate(target = geneids_2)


model_cont_group <- mod_sum_group %>%
  inner_join(mod_eval_group, by = "target")

saveRDS(model_cont_group, "data_new/models/preExc_group_only_model.RDS")




# Model age and sex without looking at interaction

args<- list(formula = y ~  age+sex + (1|study) +(1|participant), 
            family = glmmTMB::beta_family())



SE_model_2 <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                    arguments = args,
                    data = all_pre_splice_cont,
                    metadata = all_pre_metadata,
                    samplename = "seq_sample_id",
                    summary_fun = sum_fun,
                    eval_fun = eval_mod,
                    exported = list(),
                    save_models = FALSE,
                    return_models = FALSE,
                    # subset = 1:10,
                    cores = ncores-2)

SE_model_2$summaries

SE_model_2$summaries$ENST00000296098.4_6_2


excl__3 <- names(which(SE_model_2$summaries == "NULL"))
geneids_3 <- names(which(SE_model_2$summaries != "NULL"))


mod_sum <- bind_rows(within(SE_model_2$summaries, rm(excl_3))) %>%
  mutate(target = rep(geneids_3, each = 3)) %>%
  subset(coef != "(Intercept)") # %>%
#   mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
# log2fc = Estimate/log(2),
# fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) 



mod_eval <- bind_rows(within(SE_model_2$evaluations, rm(excl_3)))%>%
  mutate(target = geneids_3)


model_cont <- mod_sum %>%
  inner_join(mod_eval, by = "target")





args<- list(formula = y ~  group + sex  + (1|study) +(1|participant), 
            family = glmmTMB::beta_family())

SE_group_model_2 <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                          arguments = args,
                          data = all_pre_splice_reordered,
                          metadata = all_pre_metadata,
                          samplename = "seq_sample_id",
                          summary_fun = sum_fun,
                          eval_fun = eval_mod,
                          exported = list(),
                          save_models = FALSE,
                          return_models = FALSE,
                          cores = ncores-2)

SE_group_model_2$summaries$ENST00000007516.8_2_16

excl_4 <- names(which(SE_group_model_2$summaries == "NULL"))
geneids_4 <- names(which(SE_group_model_2$summaries != "NULL"))


mod_sum_group <- bind_rows(within(SE_group_model_2$summaries, rm(excl_4))) %>%
  mutate(target = rep(geneids_4, each = 8))  %>%
  subset(coef != "(Intercept)") # %>%
#  mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
#         log2fc = Estimate/log(2),
#         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) %>%
# filter(adj.p <= 0.05)



mod_eval_group <- bind_rows(within(SE_group_model_2$evaluations, rm(excl_4)))%>%
  mutate(target = geneids_4)


model_cont_group <- mod_sum_group %>%
  inner_join(mod_eval_group, by = "target")

saveRDS(model_cont_group, "data_new/models/preExc_group_and_sex_without_interaction_model.RDS")