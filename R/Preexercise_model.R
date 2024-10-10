# Load the file with libraries
# source("R/libraries.R")
library(dplyr)
library(tidyverse)
library(seqwrap)
library(gridExtra)
library(ggpubr)


# Load the Trainome functions
# Contains functions needed for model building
source("R/Trainome_functions.R")

#Load the metadata of the the datasets

#Copd 
copd_metadata <- readRDS("data/preexercise_data/copd_preExc_metadata.RDS") 
unique(copd_metadata$time)
colnames(copd_metadata)



#Volume

volume_metadata <- readRDS("data/preexercise_data/vol_preExc_metadata.RDS")

colnames(volume_metadata)
unique(volume_metadata$time)


# Contratratrain
Contratrain_metadata <- readRDS("data/preexercise_data/ct_PreExc_metadata.RDS")

# SRP102542
SRP102542_metadata <- readRDS("data/preexercise_data/SRP102542_preExc_metadata.RDS")

#
# Merge all in one
all_pre_metadata <- rbind(copd_metadata, volume_metadata)%>%
  rbind(Contratrain_metadata) %>%

mutate(age_group = case_when(study == "copd" ~ "Old", 
                               study == "vol" ~ "Young",
                               study == "ct" ~ "Young")) %>%
  rbind(SRP102542_metadata) %>%
  
  #Copd and volume and SRP102542 have age data in decimal
  mutate(across(c("age"), round, 0)) %>%
  mutate(sex = factor(sex, levels = c("female", "male")),
         age_group = factor(age_group, levels= c("Young", "Old"))) 

 unique(all_pre_metadata$time)
 
 
ggplot(all_pre_metadata, aes(age, fill = age_group)) +
  geom_bar()+
  ggtitle("Distribution of baseline data")+
  theme(plot.title = element_text(hjust = 0.5))

ggplot(all_pre_metadata, aes(age_group, fill = age_group)) +
  geom_bar()+
  ggtitle("Distribution of baseline data")+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = 1.5)


copd_data <- readRDS("data/preexercise_data/copd_preExc_splicing_data.RDS")


volume_data <- readRDS("data/preexercise_data/vol_preExc_splicing_data.RDS")

contratrain_data <- readRDS("data/preexercise_data/ct_PreExc_splicing_data.RDS")
SRP102542_data <- readRDS("data/preexercise_data/SRP102542_preExc_splicing_data.RDS")


all_pre_splice_cont <- copd_data%>%
  inner_join(volume_data, by = "transcript_ID") %>%
  inner_join(contratrain_data, by = "transcript_ID")%>%
  inner_join(SRP102542_data, by = "transcript_ID") %>%
  drop_na()

long_df <- all_pre_splice_cont %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID) )%>%
  inner_join(all_pre_metadata, by = "seq_sample_id")



# Plot the relationship between age and splicing efficiency
long_df %>%
  group_by( age_group)%>%
  summarise(avg = mean(SE)) %>%
  ggplot(aes(avg, age_group))+
  geom_point(mapping = aes(colour = age_group))+ 
  geom_smooth()+
  ggtitle("Relationship between age and splicing efficiency") +
  xlab(" Average splicing efficiency")+
  theme(plot.title = element_text(hjust = 0.5))

#reorder the column name to match how they occur in the metadata
#Not certain it has an impact though

all_pre_splice_reordered <- all_pre_splice_cont[,  c("transcript_ID",all_pre_metadata$seq_sample_id)]

colnames(all_pre_splice_reordered)

# Check if everything matches except the transcript_id
match(colnames(all_pre_splice_reordered), all_pre_metadata$seq_sample_id)

# invert the data to create a dataframe suitable for zero inflated analyses
# splice_df_inverted <- all_pre_splice_reordered %>%
#   mutate(across(X102PreExcVLR12:X134.subj8sample4, function(x)1-x))



#convert the 1.0 to 0.999. This is becasue beta-model accepts only values between 0 and one
all_pre_splice_reordered[all_pre_splice_reordered == 1 ] <- 0.999


args<- list(formula = y ~  age*sex + (1|study) +(1|participant), 
            family = glmmTMB::beta_family())

SE_model <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
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

SE_model$summaries$ENST00000366528.3_3_1


excl <- names(which(SE_model$summaries == "NULL"))
geneids <- names(which(SE_model$summaries != "NULL"))

#Remove all that have output NULL



mod_sum <- bind_rows(within(SE_model$summaries, rm(excl))) %>%
  mutate(target = rep(geneids, each = 4)) %>%
  subset(coef != "(Intercept)") # %>%
#   mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
# log2fc = Estimate/log(2),
# fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) 



mod_eval <- bind_rows(within(SE_model$evaluations, rm(excl)))%>%
  mutate(target = geneids)


model_cont <- mod_sum %>%
  inner_join(mod_eval, by = "target")

hist(model_cont$Pr...z..)

 saveRDS(model_cont, "data/re_models/primary_preExc_count_interaction_model.RDS")




# Model using age_group


args<- list(formula = y ~  age_group*sex + (1|study) +(1|participant), 
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
  mutate(target = rep(geneids, each = 4))  %>%
  subset(coef != "(Intercept)") #%>%
  # mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
  #        log2fc = Estimate/log(2),
  #        fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) 



mod_eval_group <- bind_rows(within(SE_count_model$evaluations, rm(excl)))%>%
  mutate(target = geneids)


model_cont_group <- mod_sum_group %>%
  inner_join(mod_eval_group, by = "target")

# This contains the interaction between age-group and sex
saveRDS(model_cont_group, "data/re_models/primary_model_extracts/primary_preExc_interaction_group_model.RDS")


# This one contained only age_group and sex, without interaction
#saveRDS(model_cont_group, "data/re_models/primary_model_extracts/primary_preExc_group_model.RDS")
