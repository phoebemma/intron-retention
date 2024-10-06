library(dplyr)
library(tidyverse)
library(seqwrap)


source("R/Trainome_functions.R")

# This contains the full model on data from Trainome group.
# It models both the model built using age as continous variable
# And that using age as a grouped variable. 

#COPD metadata
copd_metadata <- readRDS("data/processed_data/copd_metadata.RDS")
# hist(copd_metadata$age)
# max(copd_metadata$age)
# unique(copd_metadata$time)
# 
# 
# Volume_data
Vol_metadata <- readRDS("data/processed_data/volume_metadata.RDS")

# unique(Vol_metadata$time)


ct_metadata <- readRDS("data/processed_data/contratrain_metadata.RDS")

all_full_metadata <- rbind(copd_metadata, Vol_metadata)%>%
  rbind(ct_metadata) %>%
  #Copd and volume has data in decimal
  mutate(across(c("age"), round, 0)) %>%
  mutate(age_group = case_when(study == "copd" ~ "Old", 
                               study == "vol" ~ "Young",
                               study == "ct" ~ "Young")) %>%
  mutate(sex = factor(sex, levels = c("female", "male")),
         age_group = factor(age_group, levels= c("Young", "Old")),
         time = factor(time, levels = c("PreExc", "PostExc"))) 


ggplot(all_full_metadata, aes(age)) +
  geom_bar()+
  ggtitle("Distribution of full trainome data")+
  theme(plot.title = element_text(hjust = 0.5))

ggplot(all_full_metadata, aes(age_group, fill = study, colour = time)) +
  geom_bar()+
  ggtitle("Distribution of full trainome data")+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = 1.5)


# Load the splicing data

copd_splice_df <- readRDS("data/processed_data/copd_splicing_data.RDS")

vol_splice_df <- readRDS("data/processed_data/volume_splicing_data.RDS")
ct_splice_df <- readRDS("data/processed_data/contratrain_splicing_data.RDS")


all_splice_df <- copd_splice_df%>%
  inner_join(vol_splice_df, by = "transcript_ID") %>%
  inner_join(ct_splice_df, by = "transcript_ID") %>%
  drop_na()



# Visualisation

full_vis_df <- all_splice_df %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID) )%>%
  inner_join(all_full_metadata, by = "seq_sample_id")



# Plot the relationship between age and splicing efficiency
full_vis_df %>%
  group_by( time, age_group)%>%
  summarise(avg = mean(SE)) %>%
  ggplot(aes(avg, age_group))+
  geom_point(mapping = aes(colour = time, shape = age_group))+ 
  geom_smooth()+
  ggtitle("Relationship between age and splicing efficiency") +
  xlab(" Average splicing efficiency")+
  theme(plot.title = element_text(hjust = 0.5))


#select only the splicing samples captured in the metadata
splice_intersect <- intersect(colnames(all_splice_df), all_full_metadata$seq_sample_id)

all_splice_df <- all_splice_df %>%
  subset(select = c("transcript_ID", splice_intersect))



#Relace all the 1s in the dataframe to 0.99

all_splice_df[all_splice_df == 1 ] <- 0.999



args<- list(formula = y ~  age*time + sex + (1|study) +(1|participant), 
            #ziformula = ~1,
            family = glmmTMB::beta_family())


full_model <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                      arguments = args,
                      data = all_splice_df,
                      metadata = all_full_metadata,
                      samplename = "seq_sample_id",
                      summary_fun = sum_fun,
                      eval_fun = eval_mod,
                      exported = list(),
                      save_models = FALSE,
                      return_models = FALSE,
                      cores = ncores-2)
full_model$summaries$ENST00000342232.5_5_8

excl <- names(which(full_model$summaries == "NULL"))

geneids <- names(which(full_model$summaries != "NULL"))


mod_sum <- bind_rows(within(full_model$summaries, rm(excl))) %>%
  mutate(target = rep(geneids, each = 5)) %>%
  subset(coef != "(Intercept)")%>%
  mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
         log2fc = Estimate/log(2),
         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) 



mod_eval <- bind_rows(within(full_model$evaluations, rm(excl)))%>%
  mutate(target = geneids)

#merge the model evaluation and summary dataframes

model_full <- mod_sum %>%
  inner_join(mod_eval, by = "target")
saveRDS(model_full, "data/re_models/primary_full_count_model.RDS")



# Build the same model but use age as a categorical variable

args2<- list(formula = y ~  age_group*time + sex + (1|study) +(1|participant), 
            #ziformula = ~1,
            family = glmmTMB::beta_family())


full_group_model <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                      arguments = args2,
                      data = all_splice_df,
                      metadata = all_full_metadata,
                      samplename = "seq_sample_id",
                      summary_fun = sum_fun,
                      eval_fun = eval_mod,
                      exported = list(),
                      save_models = FALSE,
                      return_models = FALSE,
                      cores = ncores-2)
full_group_model$summaries$ENST00000309881.11_6_7

excl_group <- names(which(full_group_model$summaries == "NULL"))

geneids_group <- names(which(full_group_model$summaries != "NULL"))


mod_sum_group <- bind_rows(within(full_group_model$summaries, rm(excl_group))) %>%
  mutate(target = rep(geneids_group, each = 5)) # %>%
  # subset(coef != "(Intercept)")%>%
  # mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
  #        log2fc = Estimate/log(2),
  #        fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) 



mod_eval_group <- bind_rows(within(full_group_model$evaluations, rm(excl)))%>%
  mutate(target = geneids_group)

#merge the model evaluation and summary dataframes

model_full_group <- mod_sum_group %>%
  inner_join(mod_eval_group, by = "target")

saveRDS(model_full_group, "data/re_models/primary_full_group_unfiltered_model.RDS")
