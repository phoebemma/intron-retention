# This script builds and saves the baseline data
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





# Load Splicing data
all_pre_splice <- readRDS("data_new/Pre_Exercise/all_pre_Exc_splicing_data.RDS")



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
arg_1<- list(formula = y ~  scaled_age + (1|study) + (scaled_age+0|study) +(1|participant), 
            family = glmmTMB::beta_family())


model_1 <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                 arguments = arg_1,
                 data = all_pre_splice_reordered,
                 metadata = all_pre_metadata,
                 samplename = "seq_sample_id",
                 summary_fun = sum_fun2,
                 eval_fun = eval_mod,
                 exported = list(),
                 save_models = FALSE,
                 return_models = FALSE,
                 # subset = 1:10,
                 cores = ncores-2)


model_1$summaries

model_1$summaries$ENST00000023939.8_6_20


#Remove all that have output NULL
excl_1<- names(which(model_1$summaries == "NULL"))
geneids_1 <- names(which(model_1$summaries != "NULL"))


# Collect all model summaries
mod_sum_1 <- bind_rows(within(model_1$summaries, rm(excl_1))) %>%
  mutate(target = rep(geneids_1, each = 2)) %>%
   subset(coef != "(Intercept)")  %>%
  mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
         log2fc = Estimate/log(2),
         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns"))


# Bind all model evaluations
mod_eval_1 <- bind_rows(within(model_1$evaluations, rm(excl_1)))%>%
  mutate(target = geneids_1)


# Merge the evaluations and summaries
model_cont_1 <- mod_sum_1 %>%
  inner_join(mod_eval_1, by = "target") # %>%
 # filter(adj.p <= 0.05)


saveRDS(model_cont_1, "data_new/models/scaled_age_seperate_slope_intercept_model.RDS")










# The second model will group the participants into age-groups

hist(all_pre_metadata$age)

all_pre_metadata <- all_pre_metadata %>%
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


ggplot(all_pre_metadata, aes(group, fill = group )) +
  geom_bar()+
  ggtitle("Distribution of baseline data")+
  theme(plot.title = element_text(hjust = 0.5))

# Build the third model

# Initialize argument

arg_3 <- list(formula = y ~ group + (1|study)+(1|participant),
              family = glmmTMB::beta_family())


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
                   # subset = 1:10,
                   cores = ncores-2)


model_3$summaries

model_3$summaries$ENST00000023939.8_6_20


#Remove all that have output NULL
excl_3<- names(which(model_3$summaries == "NULL"))
geneids_3 <- names(which(model_3$summaries != "NULL"))


# Collect all model summaries
mod_sum_3 <- bind_rows(within(model_3$summaries, rm(excl_3))) %>%
  mutate(target = rep(geneids_3, each = 8)) %>%
  subset(coef != "(Intercept)")  %>%
  mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
         log2fc = Estimate/log(2),
         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns"))


# Bind all model evaluations
mod_eval_3 <- bind_rows(within(model_3$evaluations, rm(excl_3)))%>%
  mutate(target = geneids_3)


# Merge the evaluations and summaries
model_cont_3 <- mod_sum_3 %>%
  inner_join(mod_eval_3, by = "target") # %>%
# filter(Pr...z.. <= 0.05)

saveRDS(model_cont_3, "data_new/models/grouped_model.RDS")
