library(dplyr)
library(tidyverse)
library(seqwrap)


source("R/Trainome_functions.R")
#This is the script where the model will be built for both pre and postexercise.
#The impact of age and of RT will be modeled

#COPD metadata
# copd_metadata <- readRDS("data/processed_data/copd_metadata.RDS")%>%
#   dplyr::select(study, seq_sample_id, participant, time, age)%>%
#   subset(time == "PreExc" | time == "PostExc")%>%
#   drop_na()
#   
# hist(copd_metadata$age)
# max(copd_metadata$age)
# unique(copd_metadata$time)
# 
# 
# #Volume_data
# Vol_metadata <- readRDS("data/processed_data/volume_metadata.RDS")%>%
#   dplyr::select(study, seq_sample_id, participant, time, age)%>%
#   #pre exercise and post exercise are refered to as w0 and w12 respectively
#   subset(time == "w0"| time == "w12") %>%
#   drop_na()
# max(Vol_metadata$age)
# #change the time variables to PreExc and PostExc
# Vol_metadata["time"][Vol_metadata["time"] == "w0"] <- "PreExc"
# Vol_metadata["time"][Vol_metadata["time"] == "w12"] <- "PostExc"
# 
# unique(Vol_metadata$time)



# SRP280348
# To be extracted from this metadata are the old participants who recieved placebo
# and the young participants
SRP280348_metadata <- readRDS("data/processed_data/SRP280348_metadata.RDS") %>%
  #select only the placebo adults and the young adults
  subset(study_arm == "plaPRT" | study_arm  == "Young") %>%
  select(study, seq_sample_id, participant, biopsy, age_group, time)



#SRP102542
SRP102542_metadata <- readRDS("data/processed_data/SRP102542_metadata.RDS")%>%
  subset(exercise_type == "Resistance") %>%
  select(study, seq_sample_id, participant, time, age_group)

# SRP102542_metadata["time"][SRP102542_metadata["time"] == "PreTraining" ] <- "PreExc"
# SRP102542_metadata["time"][SRP102542_metadata["time"] == "PostTraining" ] <- "PostExc"


unique(SRP102542_metadata$time)


#SRP043368
SRP043368_metadata <- readRDS("data/processed_data/SRP043368_metadata.RDS")

SRP043368_metadata  <- SRP043368_metadata %>%
  select(study, seq_sample_id, participant, time)
  


# merge the 5 datasets

full_metadata <- rbind(copd_metadata, Vol_metadata)%>%
  rbind(SRP043368_metadata) %>%
  #create the age_group_column
  mutate(age_group = case_when(study == "copd" ~ "Old", 
                               study == "vol" ~ "Young",
                               study == "SRP043368" ~ "Young" ))%>%
  #add the dataframes that originally have the age_group information
  rbind(SRP102542_metadata)%>%
  rbind(SRP280348_metadata) %>%
  mutate(age_group = factor(age_group, levels = c("Young", "Old")))%>%
  mutate(time = factor(time, levels = c("PreExc", "PostExc")))

#create a new column to extract the age_group data
unique(full_metadata$age_group)


#saveRDS(full_metadata, "data/model/full_metadata.RDS")

#full_metadata <- readRDS("data/model/full_metadata.RDS")




ggplot(full_metadata, aes(time)) +
  geom_bar()+
  labs(x = "time")+
  ggtitle("Distribution of biopsy timepoints")+
  theme(plot.title = element_text(hjust = 0.5))


#LOAD ALL THE SPLICING DATA

copd_splice_df <- readRDS("data/processed_data/copd_splicing_data.RDS")

vol_splice_df <- readRDS("data/processed_data/volume_splicing_data.RDS")

SRP043368_splice_df <- readRDS("data/processed_data/SRP043368_splicing_data.RDS")

SRP102542_splice_df <- readRDS("data/processed_data/SRP102542_splicing_data.RDS")

SRP280348_splice_df <- readRDS("data/processed_data/SRP280348_splicing_data.RDS")


full_splice_df <- copd_splice_df%>%
  inner_join(vol_splice_df, by = "transcript_ID") %>%
  inner_join(SRP043368_splice_df, by = "transcript_ID") %>%
  inner_join(SRP102542_splice_df, by = "transcript_ID")%>%
  inner_join(SRP280348_splice_df, by = "transcript_ID") %>%
  drop_na()



#select only the splicing samples captured in the metadata
splice_intersect <- intersect(colnames(full_splice_df), full_metadata$seq_sample_id)

full_splice_df <- full_splice_df %>%
  subset(select = c("transcript_ID", splice_intersect))
#saveRDS(full_splice_df, "data/model/full_splice_data.RDS")


#Relace all the 1s in the dataframe to 0.99

full_splice_df[full_splice_df == 1 ] <- 0.999


  
args<- list(formula = y ~  age_group* time + (1|study) +(1|participant), 
            #ziformula = ~1,
            family = glmmTMB::beta_family())


full_model <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                      arguments = args,
                      data = full_splice_df,
                      metadata = full_metadata,
                      samplename = "seq_sample_id",
                      summary_fun = sum_fun,
                      eval_fun = eval_mod,
                      exported = list(),
                      save_models = FALSE,
                      return_models = FALSE,
                      cores = ncores-2)
full_model$summaries$ENST00000164139.4_11_11

#The model with ziformula =1
#saveRDS(full_model, "data/model/full_model_with_zero_inflation.RDS")



#The model without ziformula
#saveRDS(full_model, "data/model/full_model.RDS")

#determine the models that raised errors
# temp <-  full_model$errors %>%
#   
#   mutate(err = unlist(errors_fit)) %>%
#   
#   
#   
#   pivot_longer(cols = errors_fit:warn_eval) %>%
#   
#   filter(name == "err_sum") %>%
#   print()
# 
# unlist(temp$value)
# 
# 
# x = temp %>%
#   subset(value != "NULL")
# b <- as.list(x$target)
# 
#   select(target) %>%
#   tolist()


  

#x = as.list(x)

excl <- names(which(full_model$summaries == "NULL"))


#Remove all that have output NULL
full_model$summaries[which(names(full_model$summaries) %in% (excl))] <- NULL

#Remove all that have output NULL
full_model$evaluations[which(names(full_model$evaluations) %in% (excl))] <- NULL

mod_sum <- bind_rows(full_model$summaries) %>%
   mutate(target = rep(names(full_model$summaries), each = 4)) %>%
  subset(coef != "(Intercept)")%>%
  mutate(adj.p = p.adjust(Pr...z.., method = "fdr"))



mod_eval <- bind_rows(full_model$evaluations)%>%
  mutate(target = names(full_model$evaluations))


#merge the model evaluation and summary dataframes

model_full <- mod_sum %>%
  inner_join(mod_eval, by = "target")%>%
  filter(Pr...z.. <= 0.05 )%>%
  filter(Estimate >= 0.1 | Estimate <= -0.1)


hist(model_full$pval.disp)
hist(model_full$Estimate)

#saveRDS(model_full, "data/model/filtered_zero_inflated_model.RDS")
 