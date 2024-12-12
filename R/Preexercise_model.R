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
copd_metadata <- readRDS("data/preexercise_data/copd_preExc_metadata.RDS") %>%
  select(study, participant, sex, time, seq_sample_id, age)
unique(copd_metadata$time)
colnames(copd_metadata)



#Volume

volume_metadata <- readRDS("data/preexercise_data/vol_preExc_metadata.RDS")%>%
  select(study, participant, sex, time, seq_sample_id, age)

colnames(volume_metadata)
unique(volume_metadata$time)


# Contratratrain
Contratrain_metadata <- readRDS("data/preexercise_data/ct_PreExc_metadata.RDS")%>%
  select(study, participant, sex, time, seq_sample_id, age)
colnames(Contratrain_metadata)
unique(Contratrain_metadata$sex)

# SRP102542
SRP102542_metadata <- readRDS("data/preexercise_data/SRP102542_preExc_metadata.RDS")%>%
  select(study, participant, sex, time, seq_sample_id, age, age_group)
colnames(SRP102542_metadata)
unique(SRP102542_metadata$sex)

range(SRP102542_metadata$age)

SRP280348_metadata <- readRDS("data/preexercise_data/SRP280348_preExc_metadata.RDS")%>%
  select(study, participant, sex, time, seq_sample_id, age, age_group)
colnames(SRP280348_metadata)
unique(SRP280348_metadata$sex)
#
# Merge all in one
all_pre_metadata <- rbind(copd_metadata, volume_metadata)%>%
  rbind(Contratrain_metadata) %>%

mutate(age_group = case_when(study == "copd" ~ "Old", 
                               study == "vol" ~ "Young",
                               study == "ct" ~ "Young")) %>%
  rbind(SRP102542_metadata) %>%
  rbind(SRP280348_metadata) %>%
 
  #Copd and volume and SRP102542 have age data in decimal
  mutate(across(c("age"), round, 0)) %>%
  # Create youps where group 1 = those below 30
  #group 2 is those above 30 but below 51
  # group 3 those above 50 but below 71
  # group 4 is those above 70
  
  mutate(group = case_when(age <=30 ~ "<=30" ,
                           age > 30 & age <= 40 ~ ">30 & <=40", 
                           age > 40 & age <= 50 ~ ">40 & <=50",
                           age > 50 & age <= 60 ~ ">50 & <=60",
                           age > 60 & age <= 70 ~ ">60 & <=70",
                           age > 70 ~ ">70")) %>%
  mutate(sex = factor(sex, levels = c("female", "male")),
         age_group = factor(age_group, levels= c("Young", "Old")),
         group = factor(group, levels = c("<=30" , ">30 & <=40", ">40 & <=50",
                                          ">50 & <=60",">60 & <=70", ">70"))) 

 unique(all_pre_metadata$sex)
 
 
ggplot(all_pre_metadata, aes(age, fill = age_group )) +
  geom_bar()+
  ggtitle("Distribution of baseline data")+
  theme(plot.title = element_text(hjust = 0.5))


ggplot(all_pre_metadata, aes(group, fill = group)) +
  geom_bar()+
  ggtitle("Distribution of baseline data")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("Age range of participants")+
  stat_count(geom = "Text", aes(label = ..count..), vjust = 1.5)


ggplot(all_pre_metadata, aes(age_group, fill = age_group)) +
  geom_bar()+
  ggtitle("Distribution of baseline data")+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = 1.5)+
  xlab("Age group of participants")



ggplot(all_pre_metadata, aes(study, fill = group)) +
  geom_bar()+
  ggtitle("Distribution of baseline data")+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = 1.5)


ggplot(all_pre_metadata, aes(sex, fill = age_group)) +
  geom_bar()+
  ggtitle("Distribution of baseline data")+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = 1.5)

copd_data <- readRDS("data/preexercise_data/copd_preExc_splicing_data.RDS")


volume_data <- readRDS("data/preexercise_data/vol_preExc_splicing_data.RDS")

contratrain_data <- readRDS("data/preexercise_data/ct_PreExc_splicing_data.RDS")
SRP102542_data <- readRDS("data/preexercise_data/TSD_SRP102542_preExc_splicing_data.RDS")
SRP280348_data <- readRDS("data/preexercise_data/SRP280348_preExc_splicing_data.RDS")


all_pre_splice_cont <- copd_data%>%
  inner_join(volume_data, by = "transcript_ID") %>%
  inner_join(contratrain_data, by = "transcript_ID")%>%
  inner_join(SRP102542_data, by = "transcript_ID") %>%
   inner_join(SRP280348_data, by = "transcript_ID") %>%
  drop_na()

long_df <- all_pre_splice_cont %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID) )%>%
  inner_join(all_pre_metadata, by = "seq_sample_id")



# Plot the relationship between age and splicing efficiency
x <- long_df %>%
  group_by(age_group)%>%
  summarise(avg = mean(SE)) %>%
  ggplot(aes(avg, age_group))+
  geom_point(mapping = aes(colour = age_group, size = 10))+ 
  geom_smooth()+
  ggtitle("Relationship between age group and splicing efficiency") +
  xlab(" Average splicing efficiency")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_cowplot()


z <- long_df %>%
  group_by(age, sex)%>%
  summarise(avg = mean(SE)) %>%
  ggplot(aes(avg, age))+
  geom_point(mapping = aes(colour = sex ,  size = 5))+ 
  geom_smooth()+
  ggtitle("Relationship between age, gender and splicing efficiency") +
  xlab(" Average splicing efficiency")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_cowplot()
  

y <- long_df %>%
  group_by(age, age_group)%>%
  summarise(avg = mean(SE)) %>%
  ggplot(aes(avg, age))+
  geom_point(mapping = aes(colour = age_group , size = 5))+ 
  geom_smooth()+
  ggtitle("Relationship between age and splicing efficiency") +
  xlab(" Average splicing efficiency")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_cowplot()


ggarrange(y, z)



# reorder the column name to match how they occur in the metadata
# Not certain it has an impact though

# Filter to remove missing values
all_pre_metadata <- all_pre_metadata %>%
  filter((seq_sample_id %in% colnames(all_pre_splice_cont[,-1]))) %>%
  print()


all_pre_splice_reordered <- all_pre_splice_cont[ , c("transcript_ID",all_pre_metadata$seq_sample_id)]

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



# Model for group
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
  mutate(target = rep(geneids, each = 8))  %>%
  subset(coef != "(Intercept)") #%>%
# mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
#        log2fc = Estimate/log(2),
#        fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) 



mod_eval_group <- bind_rows(within(SE_count_model$evaluations, rm(excl)))%>%
  mutate(target = geneids)




model_cont_group <- mod_sum_group %>%
  inner_join(mod_eval_group, by = "target")


# This contains the interaction between age-group and sex
saveRDS(model_cont_group, "data/re_models/primary_model_extracts/primary_preExc_groups1_to_4.RDS")



x <- model_cont_group %>%
  subset(coef == "groupgroup_3")%>%
  mutate(adj.p = p.adjust(Pr...z.., method = "fdr") ,
         log2fc = Estimate/log(2),
         fcthreshold = if_else(abs(log2fc) > 1, "s", "ns")) %>%
  filter(adj.p <= 0.05 &  fcthreshold == "s" )

hist(x$Estimate)
