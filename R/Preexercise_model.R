# Load the file with libraries
# source("R/libraries.R")
library(dplyr)
library(tidyverse)
library(seqwrap)
library(gridExtra)
library(ggpubr)
library(cowplot)



# Load the Trainome functions
# Contains functions needed for model building
source("R/Trainome_functions.R")

# Load metadata
all_pre_metadata <- readRDS("data_new/Pre_Exercise/all_prexercise_metadata.RDS")




# splicing data
all_pre_splice <- readRDS("data_new/Pre_Exercise/all_pre_Exc_splicing_data.RDS")

colnames(all_pre_metadata)
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

arg_1<- list(formula = y ~  age*sex + (1|study) +(1|participant), 
            family = glmmTMB::beta_family())

# Argument for age group and interaction with sex
arg_2<- list(formula = y ~  group*sex + (1|study) +(1|participant), 
            family = glmmTMB::beta_family())
 

# Argument for age group alone
arg_3<- list(formula = y ~  group + (1|study) +(1|participant), 
            family = glmmTMB::beta_family())


# Model age and sex without looking at interaction

arg_4<- list(formula = y ~  age+sex + (1|study) +(1|participant), 
            family = glmmTMB::beta_family())


arg_5<- list(formula = y ~  group + sex  + (1|study) +(1|participant), 
            family = glmmTMB::beta_family())



args_6<- list(formula = y ~  group + sex  + (1 + group|study) +(1|participant), 
             family = glmmTMB::beta_family())


args_7<- list(formula = y ~  age_class + sex  + (1 + age_class|study) +(1|participant), 
             family = glmmTMB::beta_family())


args_8 <- list(formula = y ~  age_class + sex  + (1|study) +(1|participant), 
              family = glmmTMB::beta_family())



model_1 <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                    arguments = arg_1,
                    data = all_pre_splice_reordered,
                    metadata = all_pre_metadata,
                    samplename = "seq_sample_id",
                    summary_fun = sum_fun,
                    eval_fun = eval_mod,
                    exported = list(),
                    save_models = FALSE,
                    return_models = FALSE,
                   # subset = 1:10,
                    cores = ncores-2)

model_1$summaries

model_1$summaries$ENST00000296098.4_6_2


excl_1 <- names(which(model_1$summaries == "NULL"))
geneids_1 <- names(which(model_1$summaries != "NULL"))

#Remove all that have output NULL



mod_sum_1 <- bind_rows(within(model_1$summaries, rm(excl_1))) %>%
  mutate(target = rep(geneids_1, each = 4))# %>%
#  subset(coef != "(Intercept)") # %>%
 #   mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
 # log2fc = Estimate/log(2),
 # fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) 



mod_eval_1 <- bind_rows(within(model_1$evaluations, rm(excl_1)))%>%
  mutate(target = geneids_1)


model_cont_1 <- mod_sum_1 %>%
  inner_join(mod_eval_1, by = "target")

hist(model_cont_1$Pr...z..)

 saveRDS(model_cont_1, "data_new/models/preExc_count_interaction_model.RDS")



 
 
 
 

# Model using age_group and interaction with sex

 model_2 <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                    arguments = arg_2,
                    data = all_pre_splice_reordered,
                    metadata = all_pre_metadata,
                    samplename = "seq_sample_id",
                    summary_fun = sum_fun,
                    eval_fun = eval_mod,
                    exported = list(),
                    save_models = FALSE,
                    return_models = FALSE,
                    # subset = 1:10,
                    cores = ncores-2)
 
 model_2$summaries
 
 model_2$summaries$ENST00000296098.4_6_2
 
 
 excl_2 <- names(which(model_2$summaries == "NULL"))
 geneids_2 <- names(which(model_2$summaries != "NULL"))
 
 #Remove all that have output NULL
 
 
 
 mod_sum_2 <- bind_rows(within(model_2$summaries, rm(excl_2))) %>%
   mutate(target = rep(geneids_2, each = 14)) %>%
   subset(coef != "(Intercept)") # %>%
 #   mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
 # log2fc = Estimate/log(2),
 # fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) 
 
 
 
 mod_eval_2 <- bind_rows(within(model_2$evaluations, rm(excl_2)))%>%
   mutate(target = geneids_2)
 
 
 model_cont_2 <- mod_sum_2 %>%
   inner_join(mod_eval_2, by = "target")


# This contains the interaction between age-group and sex
saveRDS(model_cont_2, "data_new/models/preExc_group_interaction_model.RDS")











# Model for age group alone. Without sex interaction



model_3 <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                          arguments = arg_3,
                          data = all_pre_splice_reordered,
                          metadata = all_pre_metadata,
                          samplename = "seq_sample_id",
                          summary_fun = sum_fun,
                          eval_fun = eval_mod,
                          exported = list(),
                          save_models = FALSE,
                          return_models = FALSE,
                          cores = ncores-2)

model_3$summaries$ENST00000007516.8_2_16

excl_3 <- names(which(model_3$summaries == "NULL"))
geneids_3 <- names(which(model_3$summaries != "NULL"))


mod_sum_3 <- bind_rows(within(model_3$summaries, rm(excl_3))) %>%
  mutate(target = rep(geneids_3, each = 7)) %>%
  subset(coef != "(Intercept)") # %>%
#   mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
# log2fc = Estimate/log(2),
# fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) 



mod_eval_3 <- bind_rows(within(model_3$evaluations, rm(excl_3)))%>%
  mutate(target = geneids_3)


model_cont_3 <- mod_sum_3 %>%
  inner_join(mod_eval_3, by = "target")


saveRDS(model_cont_3, "data_new/models/preExc_group_only_model.RDS")



# Model age and sex without looking at interaction


model_4 <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                    arguments = arg_4,
                    data = all_pre_splice_reordered,
                    metadata = all_pre_metadata,
                    samplename = "seq_sample_id",
                    summary_fun = sum_fun,
                    eval_fun = eval_mod,
                    exported = list(),
                    save_models = FALSE,
                    return_models = FALSE,
                    # subset = 1:10,
                    cores = ncores-2)

model_4$summaries$ENST00000007516.8_2_16

excl_4 <- names(which(model_4$summaries == "NULL"))
geneids_4 <- names(which(model_4$summaries != "NULL"))


mod_sum_4 <- bind_rows(within(model_4$summaries, rm(excl_4))) %>%
  mutate(target = rep(geneids_4, each = 3)) %>%
  subset(coef != "(Intercept)") # %>%
#   mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
# log2fc = Estimate/log(2),
# fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) 



mod_eval_4 <- bind_rows(within(model_4$evaluations, rm(excl_4)))%>%
  mutate(target = geneids_4)


model_cont_4 <- mod_sum_4 %>%
  inner_join(mod_eval_4, by = "target")


saveRDS(model_cont_4, "data_new/models/PreExc_count_no_interaction_model.RDS")




# Model age group without interaction


model_5 <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                          arguments = arg_5,
                          data = all_pre_splice_reordered,
                          metadata = all_pre_metadata,
                          samplename = "seq_sample_id",
                          summary_fun = sum_fun,
                          eval_fun = eval_mod,
                          exported = list(),
                          save_models = FALSE,
                          return_models = FALSE,
                          cores = ncores-2)

model_5$summaries$ENST00000007516.8_2_16

excl_5 <- names(which(model_5$summaries == "NULL"))
geneids_5 <- names(which(model_5$summaries != "NULL"))


mod_sum_5 <- bind_rows(within(model_5$summaries, rm(excl_5))) %>%
  mutate(target = rep(geneids_5, each = 8)) %>%
  subset(coef != "(Intercept)") # %>%
#   mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
# log2fc = Estimate/log(2),
# fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) 



mod_eval_5 <- bind_rows(within(model_5$evaluations, rm(excl_5)))%>%
  mutate(target = geneids_5)


model_cont_5 <- mod_sum_5 %>%
  inner_join(mod_eval_5, by = "target")

saveRDS(model_cont_5, "data_new/models/preExc_group_and_sex_without_interaction_model.RDS")









# Model group and sex (no int) with each study having a unique group intercept


model_6 <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                            arguments = args_6,
                            data = all_pre_splice_reordered,
                            metadata = all_pre_metadata,
                            samplename = "seq_sample_id",
                            summary_fun = sum_fun,
                            eval_fun = eval_mod,
                            exported = list(),
                            save_models = FALSE,
                            return_models = FALSE,
                            cores = ncores-2)

model_6$summaries$ENST00000007516.8_2_16

excl_6 <- names(which(model_6$summaries == "NULL"))
geneids_6 <- names(which(model_6$summaries != "NULL"))

#  Model is empty
# mod_sum_6 <- bind_rows(within(model_6$summaries, rm(excl_6))) %>%
#   mutate(target = rep(geneids_6, each = 8)) %>%
#   subset(coef != "(Intercept)") # %>%
# #   mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
# # log2fc = Estimate/log(2),
# # fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) 
# 
# 
# 
# mod_eval_6 <- bind_rows(within(model_6$evaluations, rm(excl_6)))%>%
#   mutate(target = geneids_6)
# 
# 
# model_cont_6 <- mod_sum_6 %>%
#   inner_join(mod_eval_6, by = "target")

# Having an intercept unique to the study gave all NANs in summary
# saveRDS(model_cont_6, "data_new/models/PreExc_group_intercept.RDS")



# Model age class and sex without interaction but with each study having a unique intercept

model_7 <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                   arguments = args_7,
                   data = all_pre_splice_reordered,
                   metadata = all_pre_metadata,
                   samplename = "seq_sample_id",
                   summary_fun = sum_fun,
                   eval_fun = eval_mod,
                   exported = list(),
                   save_models = FALSE,
                   return_models = FALSE,
                   cores = ncores-2)

model_7$summaries$ENST00000216330.7_6_14

excl_7 <- names(which(model_7$summaries == "NULL"))
geneids_7 <- names(which(model_7$summaries != "NULL"))


mod_sum_7 <- bind_rows(within(model_7$summaries, rm(excl_7))) %>%
  mutate(target = rep(geneids_7, each = 4)) %>%
  subset(coef != "(Intercept)") # %>%
#   mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
# log2fc = Estimate/log(2),
# fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) 



mod_eval_7 <- bind_rows(within(model_7$evaluations, rm(excl_7)))%>%
  mutate(target = geneids_7)


model_cont_7 <- mod_sum_7 %>%
  inner_join(mod_eval_7, by = "target")


saveRDS(model_cont_7, "data_new/models/Preexc_age_class_no_interaction_unique_intercept.RDS")

args_8 <- list(formula = y ~  age_class + sex  + (1|study) +(1|participant), 
               family = glmmTMB::beta_family())

# Model age class and sex without interaction

model_8 <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                   arguments = args_8,
                   data = all_pre_splice_reordered,
                   metadata = all_pre_metadata,
                   samplename = "seq_sample_id",
                   summary_fun = sum_fun,
                   eval_fun = eval_mod,
                   exported = list(),
                   save_models = FALSE,
                   return_models = FALSE,
                   cores = ncores-2)

model_8$summaries$ENST00000007516.8_2_16

excl_8 <- names(which(model_8$summaries == "NULL"))
geneids_8 <- names(which(model_8$summaries != "NULL"))


mod_sum_8 <- bind_rows(within(model_8$summaries, rm(excl_8))) %>%
  mutate(target = rep(geneids_8, each = 4)) %>%
  subset(coef != "(Intercept)") # %>%
#   mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
# log2fc = Estimate/log(2),
# fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) 



mod_eval_8 <- bind_rows(within(model_8$evaluations, rm(excl_8)))%>%
  mutate(target = geneids_8)


model_cont_8 <- mod_sum_8 %>%
  inner_join(mod_eval_8, by = "target")

saveRDS(model_cont_8, "data_new/models/PreExc_age_class_sex_no_interaction.RDS")