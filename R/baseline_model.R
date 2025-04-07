# This is the updated baseline model
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

# Load metadata
all_pre_metadata <- readRDS("data_new/Pre_Exercise/all_prexercise_metadata.RDS")





# Load Splicing data that excludes the itrons with SE of 1 across all samples
all_pre_splice <- readRDS("data_new/Pre_Exercise/non_perfect_Pre_splicing_data.RDS")


# select only samples present in the metadata
all_pre_metadata <- all_pre_metadata %>%
  filter((seq_sample_id %in% colnames(all_pre_splice[,-1]))) 


# reorder the sample ids to match how they occur in the metadata
all_pre_splice_reordered <- all_pre_splice[ , c("transcript_ID",all_pre_metadata$seq_sample_id)]


# Check if everything matches except the transcript_id
match(colnames(all_pre_splice_reordered), all_pre_metadata$seq_sample_id)



# convert the 1.0 to 0.999. This is because beta-model accepts only values between 0 and one
all_pre_splice_reordered[all_pre_splice_reordered == 1 ] <- 0.999

# This argument would estimate the intercept, and the slope separately
# with uncorrelated random intercept and random slope within each study
<<<<<<< HEAD
arg_1<- list(formula = y ~  scaled_age  + (1|study) + (scaled_age+0|study) +(1|participant), 
=======
arg_1<- list(formula = y ~  scaled_age + (1|study) + (scaled_age+0|study) +(1|participant), 
>>>>>>> a5265ac924e41ee4df97d3e02786f00a60fd0943
             family = glmmTMB::beta_family())


model_1 <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                   arguments = arg_1,
                   data = all_pre_splice_reordered,
                   metadata = all_pre_metadata,
                   samplename = "seq_sample_id",
                   summary_fun = sum_fun,
                   eval_fun = eval_mod,
                   exported = list(),
                   save_models = F,
                   return_models = T,
<<<<<<< HEAD
                   #   subset = 1:100,
                   cores = ncores-2)
model_1$summaries$
=======
                   # subset = 1:10,
                   cores = ncores)

model_1$summaries$ENST00000001008.6_5_12
>>>>>>> a5265ac924e41ee4df97d3e02786f00a60fd0943


names(model_1$models)
# plot(effect_plot)

#exclude those whose summaries are not Null
model_list <- model_1$models[which(model_1$summaries != "NULL")]

baseline_predictions <- data.frame(scaled_age = numeric(),
                                   target = character(), type = character(), stringsAsFactors = FALSE)

for (i in 1:length(model_list)) {
  model_name <- names(model_list)[i]
  
  # Extract the effect of predictors and catch those with null
  effect_plot <- # tryCatch({
    allEffects(model_list[[i]], xlevels=list(scaled_age=seq(from = 0, to = 1, by = 0.1)))

    effect_df <- as.data.frame(effect_plot$scaled_age)
    
    # Add a column for the model name
    effect_df$target <- model_name

    
    # Append to the baseline_predictions dataframe
    baseline_predictions <- rbind(baseline_predictions, effect_df)
}


<<<<<<< HEAD


saveRDS(baseline_predictions, "data_new/baseline_model_predictions.RDS")
=======
saveRDS(model_1, "data_new/baseline_model.RDS")
>>>>>>> a5265ac924e41ee4df97d3e02786f00a60fd0943

saveRDS(baseline_predictions, "data_new/baseline_predictions.RDS")

#Remove all that have output NULL
excl_1<- names(which(model_1$summaries == "NULL"))
geneids_1 <- names(which(model_1$summaries != "NULL"))


# Collect all model summaries
mod_sum_1 <- bind_rows(within(model_1$summaries, rm(excl_1))) %>%
  mutate(target = rep(geneids_1, each = 2)) %>%
#  subset(coef != "(Intercept)")  %>%
  mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
         log2fc = Estimate/log(2),
         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns"))


<<<<<<< HEAD
# Bind all model evaluations
# mod_eval_1 <- bind_rows(within(model_1$evaluations, rm(excl_1)))%>%
#   mutate(target = geneids_1)
# 
# 
# # Merge the evaluations and summaries
# model_cont_1 <- mod_sum_1 %>%
#   inner_join(mod_eval_1, by = "target") # %>%
# filter(adj.p <= 0.05)


saveRDS(mod_sum_1, "data_new/baseline_model_summary.RDS")
=======

saveRDS(mod_sum_1 , "data_new/baseline_model_summary.RDS")
>>>>>>> a5265ac924e41ee4df97d3e02786f00a60fd0943


