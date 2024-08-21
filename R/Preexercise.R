#Load the file with libraries
#source("R/libraries.R")
library(dplyr)
library(tidyverse)
library(seqwrap)

#Load the trainome functions file
source("R/Trainome_functions.R")


#Load the metadata of the the datasets

#Copd 
copd_metadata <- readRDS("data/processed_data/copd_metadata.RDS")

unique(copd_metadata$time)

#Prexercise in COPD is named "PreExc"

PreEXC_copd_meta <- copd_metadata %>%
  subset(time == "PreExc") %>%
  select(study, seq_sample_id, participant) %>%
  drop_na()

#saveRDS(PreEXC_copd_meta, "data/preexercise_data/preexc_copd_metadata.RDS")



#Volume

volume_metadata <- readRDS("data/processed_data/volume_metadata.RDS")

unique(volume_metadata$time)

#Prexercise data in volume is "w0" 

PreExc_vol_meta <- volume_metadata %>%
  subset(time == "w0") %>%
  select(study, seq_sample_id, participant)%>%
  drop_na()

#saveRDS(PreExc_vol_meta, "data/preexercise_data/preexc_volume_metadata.RDS")

#Contratrain
ct_metadata <- readRDS("data/processed_data/contratrain_metadata.RDS")

unique(ct_metadata$time)

#Subset to the Pre-exercise data

PreExc_ct_meta <- ct_metadata %>%
  subset(time == "t1") %>%
  select(study, seq_sample_id, participant) %>%
  drop_na()

#saveRDS(PreExc_ct_meta, "data/preexercise_data/preexc_contratrain_metadata.RDS")
#SRP280348 

SRP280348_meta <- readRDS("data/processed_data/SRP280348_metadata.RDS")

colnames(SRP280348_meta)
#Subset the preexercise  data
#This is the biopsy data 1 for old participants, and 0 for young participants

PreExc_SRP280348_meta <- SRP280348_meta %>%
  subset(biopsy == 1 | biopsy == 0) %>%
  select(study, seq_sample_id, participant, age_group)

#saveRDS(PreExc_SRP280348_meta, "data/preexercise_data/preexc_SRP280348_metadata.RDS")



#SRP043368 is baseline data alone, would be loaded from the processed data folder

SRP043368 <- readRDS("data/processed_data/SRP043368_metadata.RDS")%>%
  select(study, seq_sample_id, participant)


colnames(SRP043368)

#SRP102542

SRP102542_meta <- readRDS("data/processed_data/SRP102542_metadata.RDS")

PreExc_SRP102542_meta <- SRP102542_meta %>%
  subset(biopsy_timepoint == "PreTraining")%>%
  select(study, seq_sample_id, participant, age_group)

#saveRDS(PreExc_SRP102542_meta, "data/preexercise_data/preexc_SRP102542_metadata.RDS")


#combine all metadata into one


#The first three datasets do not have a column called "age_group".
#They will be merged and the column created


all_metadata <- rbind(PreEXC_copd_meta, PreExc_vol_meta)%>%
  rbind(SRP043368) %>%
  #create the age_group_column
  mutate(age_group = case_when(study == "copd" ~ "Old", 
                               study == "vol" ~ "Young",
                               study == "SRP043368" ~ "Young" ))%>%
  #add the dataframes that originally have the age_group information
  rbind(PreExc_SRP102542_meta)%>%
  rbind(PreExc_SRP280348_meta) %>%
  mutate(age_group = factor(age_group, levels = c("Young", "Old")))
  
#create a new column to extract the age_group data
unique(all_metadata$age_group)

#saveRDS(all_metadata, "data/preexercise_data/all_metadata.RDS")

# Post_ex_df <- fullage_group# Post_ex_df <- full_splice_df %>%
#   subset(time == "PostExc")
ggplot(all_metadata, aes(age_group)) +
  geom_bar()


#Load each of the splicing data

copd_data <- readRDS("data/processed_data//copd_splicing_data.RDS")

#get the preexercise data
splice_intersect <- (intersect(colnames(copd_data),
                               PreEXC_copd_meta$seq_sample_id))

#subset the splicing data to only include the intersects
pre_copd_data <- copd_data %>%
  subset( select = c("transcript_ID", splice_intersect))

#saveRDS(pre_copd_data, "data/preexercise_data/preexc_copd_data.RDS")


volume_data <- readRDS("data/processed_data/volume_splicing_data.RDS")

splice_intersect <- (intersect(colnames(volume_data),
                               PreExc_vol_meta$seq_sample_id))

#subset the splicing data to only include the intersects
pre_vol_data <- volume_data %>%
  subset( select = c("transcript_ID", splice_intersect))

#saveRDS(pre_vol_data, "data/preexercise_data/prexc_vol_data.RDS")


SRP043368_data <- readRDS("data/processed_data/SRP043368_splicing_data.RDS")
SRP102542_data <- readRDS("data/processed_data/SRP102542_splicing_data.RDS")


splice_intersect <- (intersect(colnames(SRP102542_data),
                               PreExc_SRP102542_meta$seq_sample_id))

pre_SRP102542_data <- SRP102542_data %>%
  subset( select = c("transcript_ID", splice_intersect))

#saveRDS(pre_SRP102542_data, "data/preexercise_data/pre_SRP102542_data.RDS")


SRP280348_data <- readRDS("data/processed_data/SRP280348_splicing_data.RDS")
splice_intersect <- (intersect(colnames(SRP280348_data),
                               PreExc_SRP280348_meta$seq_sample_id))

pre_SRP280348_data <- SRP280348_data %>%
  subset( select = c("transcript_ID", splice_intersect))

#saveRDS(pre_SRP280348_data, "data/preexercise_data/pre_SRP280348_data.RDS")

#merge all the 5 splicing data
all_splice_df <- pre_copd_data%>%
  inner_join(pre_vol_data, by = "transcript_ID") %>%
  inner_join(pre_SRP102542_data, by = "transcript_ID") %>%
  inner_join(pre_SRP280348_data, by = "transcript_ID")%>%
  inner_join(SRP043368_data, by = "transcript_ID") %>%
  drop_na()


#select only the splicing samples captured in the metadata
splice_intersect <- intersect(colnames(all_splice_df), all_metadata$seq_sample_id)

all_splice_df <- all_splice_df %>%
  subset(select = c("transcript_ID", splice_intersect))
#saveRDS(all_splice_df, "data/preexercise_data/all_splice_data.RDS")


all_splice_df[all_splice_df == 1 ] <- 0.999






#Subset only the introns wih SE 0
colnames(all_splice_df)
#invert the data to make it suitable for zero inflated analyses
# splice_df_inverted <- all_splice_df %>%
#   mutate(across(X144PreExcVLL198:SRR1424756, function(x)1-x))

#Relace all the 1s in the dataframe to 0.99




#ziformula of 1 suggests there is zero inflation

args<- list(formula = y ~  age_group + (1|study) +(1|participant), 
            #ziformula = ~1,
            family = glmmTMB::beta_family())




SE_model <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                         arguments = args,
                         data = all_splice_df,
                         metadata = all_metadata,
                         samplename = "seq_sample_id",
                         summary_fun = sum_fun,
                         eval_fun = eval_mod,
                         exported = list(),
                         save_models = FALSE,
                         return_models = FALSE,
                         cores = ncores-2)


#saveRDS(SE_model, "data/model/preexercise_model.RDS")

SE_model$summaries$ENST00000164139.4_11_11

excl <- names(which(SE_model$summaries == "NULL"))


#Remove all that have output NULL
SE_model$summaries[which(names(SE_model$summaries) %in% (excl))] <- NULL

#Remove all that have output NULL
SE_model$evaluations[which(names(SE_model$evaluations) %in% (excl))] <- NULL

mod_sum <- bind_rows(SE_model$summaries) %>%
  mutate(target = rep(names(SE_model$summaries), each = 2)) %>%
  subset(coef != "(Intercept)")%>%
  mutate(adj.p = p.adjust(Pr...z.., method = "fdr"))



mod_eval <- bind_rows(SE_model$evaluations)%>%
  mutate(target = names(SE_model$evaluations))


#merge the model evaluation and summary dataframes

model_full <- mod_sum %>%
  inner_join(mod_eval, by = "target")%>%
  filter(Pr...z.. <= 0.05 )%>%
  filter(Estimate >= 0.1 | Estimate <= - 0.1)



hist(model_full$pval.disp)
hist(model_full$Estimate)


unique(model_full$Estimate)

#saveRDS(model_full, "data/model/filtered_presercise_model.RDS")






