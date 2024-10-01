library(dplyr)
library(tidyverse)
library(seqwrap)


source("R/Trainome_functions.R")
#This is the script where the model will be built for both pre and postexercise.
#The impact of age and of RT will be modeled

#COPD metadata
copd_metadata <- readRDS("data/processed_data/copd_metadata.RDS")%>%
  dplyr::select(study, participant,time, seq_sample_id)

# hist(copd_metadata$age)
# max(copd_metadata$age)
# unique(copd_metadata$time)
# 
# 
# Volume_data
Vol_metadata <- readRDS("data/processed_data/volume_metadata.RDS")%>%
  dplyr::select(study, participant,time, seq_sample_id)

  
# unique(Vol_metadata$time)


ct_metadata <- readRDS("data/processed_data/contratrain_metadata.RDS")%>%
  dplyr::select(study, participant,time, seq_sample_id)

# SRP280348
# To be extracted from this metadata are the old participants who recieved placebo
# and the young participants
SRP280348_metadata <- readRDS("data/processed_data/SRP280348_metadata.RDS") %>%
  #select only the placebo adults and the young adults
  subset(study_arm == "plaPRT" | study_arm  == "Young") %>%
  select(study, participant,time, seq_sample_id, age_group)

unique(SRP280348_metadata$time)




#SRP102542

SRP102542_metadata <- readRDS("data/processed_data/SRP102542_metadata.RDS")%>%
  subset(exercise_type == "Resistance") %>%
  select(study, participant,time, seq_sample_id, age_group)



unique(SRP102542_metadata$participant)


#SRP043368
# As it is only baseline, it is contained in the preexercise data
SRP043368_metadata <- readRDS("data/preexercise_data/SRP043368_metadata.RDS")

SRP043368_metadata  <- SRP043368_metadata %>%
  select(study, participant,time, seq_sample_id)
  


# merge the 5 datasets

full_metadata <- rbind(copd_metadata, Vol_metadata)%>%
  rbind(ct_metadata) %>%
  rbind(SRP043368_metadata) %>%
  #create the age_group_column
  mutate(age_group = case_when(study == "copd" ~ "Old", 
                               study == "vol" ~ "Young",
                               study == "ct" ~ "Young",
                               study == "SRP043368" ~ "Young" ))%>%
  #add the dataframes that originally have the age_group information
  rbind(SRP102542_metadata)%>%
  rbind(SRP280348_metadata) %>%
  mutate(age_group = factor(age_group, levels = c("Young", "Old")))%>%
  mutate(time = factor(time, levels = c("PreExc", "PostExc")))

#create a new column to extract the age_group data
unique(full_metadata$age_group)
unique(full_metadata$study)

#saveRDS(full_metadata, "data/model/full_metadata.RDS")

#full_metadata <- readRDS("data/model/full_metadata.RDS")




ggplot(full_metadata, aes(time, fill = time)) +
  geom_bar()+
  labs(x = "time")+
  ggtitle("Distribution of biopsy timepoints for full model")+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = -0.5)

#LOAD ALL THE SPLICING DATA

copd_splice_df <- readRDS("data/processed_data/copd_splicing_data.RDS")

vol_splice_df <- readRDS("data/processed_data/volume_splicing_data.RDS")
ct_splice_df <- readRDS("data/processed_data/contratrain_splicing_data.RDS")

SRP043368_splice_df <- readRDS("data/preexercise_data/SRP043368_splicing_data.RDS")

SRP102542_splice_df <- readRDS("data/processed_data/SRP102542_splicing_data.RDS")

SRP280348_splice_df <- readRDS("data/processed_data/SRP280348_splicing_data.RDS")


full_splice_df <- copd_splice_df%>%
  inner_join(vol_splice_df, by = "transcript_ID") %>%
  inner_join(ct_splice_df, by = "transcript_ID") %>%
  inner_join(SRP043368_splice_df, by = "transcript_ID") %>%
  inner_join(SRP102542_splice_df, by = "transcript_ID")%>%
  inner_join(SRP280348_splice_df, by = "transcript_ID") %>%
  drop_na()

# Visualisation

vis_df <- full_splice_df %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID) )%>%
  inner_join(full_metadata, by = "seq_sample_id")



# Plot the relationship between age and splicing efficiency
vis_df %>%
  group_by( time, age_group)%>%
  summarise(avg = mean(SE)) %>%
  ggplot(aes(avg, age_group))+
  geom_point(mapping = aes(colour = age_group, shape = time))+ 
  #geom_smooth()+
  ggtitle("Relationship between age and splicing efficiency") +
  xlab(" Average splicing efficiency")+
  theme(plot.title = element_text(hjust = 0.5))






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





excl <- names(which(full_model$summaries == "NULL"))

geneids <- names(which(full_model$summaries != "NULL"))


mod_sum <- bind_rows(within(full_model$summaries, rm(excl))) %>%
  mutate(target = rep(geneids, each = 4)) %>%
  subset(coef != "(Intercept)")%>%
  mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
         log2fc = Estimate/log(2),
         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) 



mod_eval <- bind_rows(within(full_model$evaluations, rm(excl)))%>%
  mutate(target = geneids)

#merge the model evaluation and summary dataframes

model_full <- mod_sum %>%
  inner_join(mod_eval, by = "target")
hist(model_full$pval.disp)
hist(model_full$Estimate)

saveRDS(model_full, "data/re_models/all_full_group_model.RDS")
 