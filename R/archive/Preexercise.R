#This script contains the analyses and modeling of the prexercise data

#Load the file with libraries
#source("R/libraries.R")
library(dplyr)
library(tidyverse)
library(seqwrap)
library(gridExtra)
library(ggpubr)

#Load the trainome functions file
source("R/Trainome_functions.R")


#Load the metadata of the the datasets

#Copd 
copd_metadata <- readRDS("data/preexercise_data/copd_preExc_metadata.RDS") %>%
  select(study, participant,time, seq_sample_id )
# 
unique(copd_metadata$time)




#Volume

volume_metadata <- readRDS("data/preexercise_data/vol_preExc_metadata.RDS")%>%
  select(study, participant, time, seq_sample_id )
# 


unique(volume_metadata$time)


# Contratratrain
Contratrain_metadata <- readRDS("data/preexercise_data/ct_PreExc_metadata.RDS")%>%
  select(study, participant, time, seq_sample_id )
# 



# #SRP043368
# 
 SRP043368_metadata <- readRDS("data/preexercise_data/SRP043368_metadata.RDS")%>%
   select(study, participant, time, seq_sample_id )
# 
# 
# colnames(SRP043368_metadata)
# 
# #SRP280348
# 
 SRP280348_metadata <- readRDS("data/preexercise_data/SRP280348_preExc_metadata.RDS")%>%
   dplyr::select(study, participant, time, seq_sample_id,  age_group)
# 
 unique(SRP280348_metadata$age_group)
# 
# 
# 
# 
# #SRP102542
# 
 SRP102542_metadata <- readRDS("data/preexercise_data/SRP102542_preExc_metadata.RDS")%>%
   select(study, participant, time, seq_sample_id, age_group )
# 
# unique(SRP102542_metadata$time)

SRP102542_metadata$age_group
#combine all metadata into one


#The first three datasets do not have a column called "age_group".
#They will be merged and the column created


all_pre_metadata <- rbind(copd_metadata, volume_metadata)%>%
  rbind(Contratrain_metadata) %>%
  rbind(SRP043368_metadata) %>%
  #create the age_group_column
  mutate(age_group = case_when(study == "copd" ~ "Old", 
                               study == "vol" ~ "Young",
                               study == "ct" ~ "Young",
                               study == "SRP043368" ~ "Young"))  %>%
  #add the dataframes that originally have the age_group information
   rbind(SRP102542_metadata)  %>%
   rbind(SRP280348_metadata) #%>%
   mutate(age_group = factor(age_group, levels = c("Young", "Old")),
          sex = factor(sex, levels = c("female", "male")))
  
#create a new column to extract the age_group data
unique(all_pre_metadata$age_group)

#saveRDS(all_pre_metadata, "data/preexercise_data/all_preExc_metadata.RDS")

#Plot distribution of age_groups
ggplot(all_pre_metadata, aes(age_group, fill = study)) +
  geom_bar()+
  ggtitle("Distribution of baseline data")+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = 1.5)



#Load each of the splicing data
#The pre-exercise splicing data has been extracted and saved in the preexercise folder

copd_data <- readRDS("data/preexercise_data/copd_preExc_splicing_data.RDS")


volume_data <- readRDS("data/preexercise_data/vol_preExc_splicing_data.RDS")

contratrain_data <- readRDS("data/preexercise_data/ct_PreExc_splicing_data.RDS")


 SRP043368_data <- readRDS("data/preexercise_data/SRP043368_splicing_data.RDS")
# 
# 
 SRP102542_data <- readRDS("data/preexercise_data/SRP102542_preExc_splicing_data.RDS")
# 
# 
 SRP280348_data <- readRDS("data/preexercise_data/SRP280348_preExc_splicing_data.RDS")



#merge all the 5 splicing data
all_pre_splice <- copd_data%>%
  inner_join(volume_data, by = "transcript_ID") %>%
  inner_join(contratrain_data, by = "transcript_ID") %>%
   inner_join(SRP043368_data, by = "transcript_ID") %>%
   inner_join(SRP102542_data, by = "transcript_ID") %>%
   inner_join(SRP280348_data, by = "transcript_ID")%>%
  drop_na()



# saveRDS(all_pre_splice, "data/preexercise_data/all_preExc_splice_data.RDS")


#Visualization
long_df <- all_pre_splice %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID) )%>%
  inner_join(all_pre_metadata, by = "seq_sample_id") # %>%
  #create a new column that states whether intron is retained or not
  # mutate(status = if_else(SE == 1, "none_retained", "retained"),
  #        reverse_status = case_when(SE == 0.0 ~ "completely_retained", 
  #                                   SE > 0.0 & SE < 1.0 ~ "partially_retained", 
  #                                   SE == 1.0 ~ "completely_spliced")) %>%
  # # copd and volume ages are in decimal places, remove
  # mutate(across(c("age"), round, 0))
  # 






# x <- SRP280348_data %>%
#   pivot_longer(names_to = "seq_sample_id",
#                values_to = "SE",
#                cols = -(transcript_ID) )%>%
#   inner_join(SRP280348_metadata, by = "seq_sample_id") 
# 
# unique(x$SE)
# x %>%
#   group_by( age_group)%>%
#   summarise(avg = mean(SE, na.rm = T)) %>%
#   ggplot(aes(avg, age_group))+
#   geom_point(mapping = aes(colour = age_group))+ 
#   geom_smooth()



long_df %>%
   group_by( age_group)%>%
   summarise(avg = mean(SE, na.rm = T)) %>%
  ggplot(aes(avg, age_group))+
  geom_point(mapping = aes(colour = age_group))+ 
  geom_smooth()

cor(long_df$SE, long_df$age)


cor(long_df$SE, long_df$age, method = "spearman")


cor(long_df$SE, long_df$age, method =  "kendall")

long_df %>%
  ggplot(aes(SE, age_group))+
  geom_boxplot()


plot(long_df$age,long_df$SE )





old_pre_df <- long_df %>%
  subset(age_group == "Old") %>%
  group_by(reverse_status)%>%
  count()%>% #count the occurence
  ungroup()%>%
  mutate(perc = n/sum(n))%>% #extract the percentage
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))


young_pre_df <- long_df %>%
  subset(age_group == "Young") %>%
  group_by(reverse_status)%>%
  count()%>% #count the occurence
  ungroup()%>%
  mutate(perc = n/sum(n))%>% #extract the percentage
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))


old_pre <- ggplot(old_pre_df, aes(x = "", y = perc, fill = reverse_status)) +
  geom_col()+
  geom_text(aes(label = labels, x = 1.6),
            position = position_stack(vjust = 0.8)) +
  coord_polar(theta = "y")+
  theme_void()+
  ggtitle("percentage of  introns retained in old participants at baseline")

young_pre <- ggplot(young_pre_df, aes(x = "", y = perc, fill = reverse_status)) +
  geom_col()+
  geom_text(aes(label = labels, x = 1.6),
            position = position_stack(vjust = 0.8)) +
  coord_polar(theta = "y")+
  theme_void()+
  ggtitle("percentage of  introns retained in young participants at baseline")




ggarrange(old_pre, young_pre, #you can also specify rremove("x.text) to remove texts on the x axis 
          labels = c("A", "B") #Add labels like done in publications
        ) 

colnames(all_pre_splice)

order
#reorder the column name to match how they occur in the metadata
#Not certain it has an impact though

all_pre_splice_reordered <- all_pre_splice[,  c("transcript_ID",all_pre_metadata$seq_sample_id)]

colnames(all_pre_splice_reordered)

#Check if everything matches except the transcript_id
match(colnames(all_pre_splice_reordered), all_pre_metadata$seq_sample_id)

#invert the data to create a dataframe suitable for zero inflated analyses
 splice_df_inverted <- all_pre_splice_reordered %>%
   mutate(across(X102PreExcVLR12:SRR12604227, function(x)1-x))



#convert the 1.0 to 0.999. This is becasue beta-model accepts only values between 0 and one
all_pre_splice_reordered[all_pre_splice_reordered == 1 ] <- 0.999







args<- list(formula = y ~  age_group + (1|study) +(1|participant), 
            family = glmmTMB::beta_family())
#ziformula of 1 suggests there is zero inflation
# args_2<- list(formula = y ~  age_group + (1|study) +(1|participant), 
#             ziformula = ~1,
#             family = glmmTMB::beta_family())


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

# 
# zero_inf_mod <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
#                         arguments = args_2,
#                         data = splice_df_inverted,
#                         metadata = all_pre_metadata,
#                         samplename = "seq_sample_id",
#                         summary_fun = sum_fun,
#                         eval_fun = eval_mod,
#                         exported = list(),
#                         save_models = FALSE,
#                         return_models = FALSE,
#                         cores = ncores-2)





SE_model$summaries$ENST00000164139.4_12_11
# zero_inf_mod$summaries$ENST00000164139.4_12_11

excl <- names(which(SE_model$summaries == "NULL"))
geneids <- names(which(SE_model$summaries != "NULL"))



mod_sum <- bind_rows(within(SE_model$summaries, rm(excl))) %>%
  mutate(target = rep(geneids, each = 2)) %>%
  subset(coef != "(Intercept)")%>%
  mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
         log2fc = Estimate/log(2),
         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) 



mod_eval <- bind_rows(within(SE_model$evaluations, rm(excl)))%>%
  mutate(target = geneids)


#merge the model evaluation and summary dataframes

model_full <- mod_sum %>%
  inner_join(mod_eval, by = "target")

saveRDS(model_full, "data/re_models/all_preExc_group_model,RDS")








