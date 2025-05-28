library(dplyr)
library(tidyverse)
library(seqwrap)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(effects)
library(scales)
library(glmmTMB)


# Load the Trainome functions
# Contains functions needed for model building
source("R/Trainome_functions.R")



all_splice_df <- readRDS("data/all_splice.RDS") %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  drop_na() 


# Extract the introns that are perfectly spliced across all samples

# Get the perfectly spliced introns across all datasets
# High_SE_df <- all_splice_df %>%
#   pivot_longer(names_to = "seq_sample_id",
#                values_to = "SE",
#                cols = -(transcript_ID)) %>%
#   summarise(.by = transcript_ID, 
#             mode = getmode(SE),
#             min = min(SE), 
#             max = max(SE)) %>%
#   filter(min == 1) 

all_full_metadata <- readRDS("data_new/processed_data/all_full_metadata.RDS")
colnames(all_full_metadata)


# REORDER THE SEQUENCE ID TO MATCH BOTH DATAFRAMMES
all_splice_reordered <- all_splice_df[, c("transcript_ID",all_full_metadata$seq_sample_id)]

# Check if everything matches except the transcript_id
match(colnames(all_splice_reordered), all_full_metadata$seq_sample_id)


# convert the 1.0 to 0.999. This is becasue beta-model accepts only values between 0 and one
all_splice_reordered[all_splice_reordered == 1 ] <- 0.999



# model using Seqwrap


args_full <-list(formula = y ~  scaled_age * time + sex + (1|study) +(1|participant), 
                 family = glmmTMB::beta_family())


container <- seqwrap_compose(data = all_splice_reordered,
                             metadata = all_full_metadata,
                             samplename = "seq_sample_id",
                             modelfun = glmmTMB::glmmTMB,
                             arguments = args_full,
                             summary_fun = sum_with_pred,
                             eval_fun = eval_mod )



model <- seqwrap(container,
                 summary_fun = sum_with_pred,
                 eval_fun = eval_mod,
                 return_models = F,
                 # subset = 1:15,
                 cores = ncores-2)

model_summary <- seqwrap_summarise(model)
head(model_summary)
str(model)
saveRDS(model_summary, "data/RT_model_summary.RDS")
# extract the model summaries
mod_sum <- model_summary$summaries

RT_model_summary <- mod_sum %>%
  dplyr::select(coef, target,Estimate, Pr...z..) %>%
  # select the model results of those affected by RT alone or interaction with aging
  filter(coef == "timePostExc" | coef == "scaled_age:timePostExc" | coef == "scaled_age") %>%
  mutate(adj.p = p.adjust( Pr...z.., method = "fdr"),
         log2fc = Estimate/log(2),
         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns"),
         odds_ratio = exp(Estimate))

RT_predictions <- mod_sum %>%
  dplyr::select(scaled_age, time, fit, target) %>%
  drop_na() %>%
  # Extract the transcript_id from the target
  mutate(transcript_ID = str_split(target, "_",simplify= T) [,1])


# Load the gene annotation file
gene_annotation <- readRDS("data_new/ensembl_gene_annotation.RDS")


# merge dataset
RT_merged <- RT_predictions %>%
  inner_join(RT_model_summary, by = "target") %>%
  drop_na() %>%
  mutate(effect = case_when(Estimate > 0 & Pr...z.. <= 0.05 ~ "Improved SE",
                            Estimate < 0 & Pr...z.. <= 0.05 ~ "Reduced SE" ,
                            Estimate < 0 & Pr...z.. > 0.05 ~ "No effect",
                            Estimate > 0 & Pr...z.. > 0.05 ~ "No effect")) %>%
  inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version"))


saveRDS(RT_merged, "data/RT_model_df.RDS")




# plot the distribution of gene biotypes
RT_merged %>%
  ggplot(aes(transcript_biotype, fill = transcript_biotype))+
  geom_bar()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.title = element_text(hjust = 0.5))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = -0.5) +
  ggtitle("Distribution of gene biotypes in our dataset") +
  ylab("Number of introns per biotype") +
  facet_grid(~effect)




# Explore tyhe effect of aging alone
Aging_effect <- RT_merged %>%
  filter(coef == "scaled_age" & time == "PreExc")


# extract the disicnt effect types
summary_aging <- Aging_effect %>%
  group_by(effect) %>%
  summarize(num_targets = n_distinct(target))

# get a dataset that includes a legen that displays the number of introns in each group
Aging_effect <- Aging_effect %>%
  left_join(summary_aging, by = "effect")%>%
  mutate(effect_label = paste(effect, "(No of introns:", num_targets, ")"))



aging_plot <- ggplot(Aging_effect, aes(x = scaled_age, y = fit, color = effect_label)) +
  geom_smooth(method = "lm", se = FALSE) + # se = FALSE to remove confidence intervals
  theme_minimal() +
  labs(title = "Relationship between aging and Splicing Efficiency at baseline", 
       x = "Scaled Age", 
       y = "Splicing Efficiency") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 10, face = "italic"),
        legend.title = element_text(size = 10, face = "bold", hjust = 0.5))



# Investigate chromosome level infoprmation
chr_RT <- RT_model_summary %>%
  filter(coef == "scaled_age" ) %>%
  drop_na() %>%
  mutate(effect = case_when(Estimate > 0 & Pr...z.. <= 0.05 ~ "Improved SE",
                            Estimate < 0 & Pr...z.. <= 0.05 ~ "Reduced SE" ,
                            Estimate < 0 & Pr...z.. > 0.05 ~ "No effect",
                            Estimate > 0 & Pr...z.. > 0.05 ~ "No effect"),
         chromosome = str_extract(target, "(?<=_)[^_]+$"),
         chromosome = factor(chromosome, levels = c("X", "1", "2", "3", "4", "5",
                                                    "6", "7", "8", "9", "10", "11",
                                                    "12", "13", "14", "15", "16", "17",
                                                    "18", "19", "20", "21", "22")))


chromosome_distribution <- chr_RT %>%
  # filter(effect != "No effect" ) %>% 
  group_by(target) %>%
  ggplot(aes(chromosome))+
  geom_bar()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.title = element_text(hjust = 0.5))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = -0.5) +
  ggtitle("Distribution of introns across chromosomes") +
  ylab("Number of introns per chromosome") +
  facet_grid(~effect)


# Extract the most affected introns using  fit values
doubled_aging <- Aging_effect %>%
  # select for those significantly affected by aging
  filter(fcthreshold == "s") %>%
  arrange(fit) #%>%
#   slice_head(n = 200) %>%
#   dplyr::select(target, external_gene_name, time, fit, Estimate,transcript_biotype, scaled_age, external_transcript_name) 
length(unique(doubled_aging$external_gene_name))   