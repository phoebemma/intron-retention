library(dplyr)
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(ggridges)
library(UpSetR)
library(ggplot2)
library(ggplotify)
library(clusterProfiler)
library(org.Hs.eg.db)

# Load the postexercise model summary

RT_model_summary <- readRDS("data_new/simpler_RT_model_summary.RDS") %>%
  # select the model results of those affected by RT alone or interaction with aging
  filter(coef == "timePostExc" | coef == "scaled_age:timePostExc" | coef == "scaled_age")%>%
  # add the odds ratio for each Estimate
  mutate(odds_ratio = exp(Estimate),
         type = "model coefficient")


# load the predictions 

RT_predictions <- readRDS("data_new/predictions_simpler_RT_model.RDS")%>%
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

saveRDS(RT_merged, "data_new/RT_model_outputs/RT_model_df.RDS")



# Explore tyhe effect of aging alone
Aging_effect <- RT_merged %>%
  filter(coef == "scaled_age")

# extract the disicnt effect types
summary_aging <- Aging_effect %>%
  group_by(effect) %>%
  summarize(num_targets = n_distinct(target))

Aging_effect <- Aging_effect %>%
  left_join(summary_aging, by = "effect")%>%
  mutate(effect_label = paste(effect, "(No of introns:", num_targets, ")"))



ggplot(Aging_effect %>%
         filter(time == "PreExc"), aes(x = scaled_age, y = fit, color = effect_label)) +
  geom_smooth(method = "lm", se = FALSE) + # se = FALSE to remove confidence intervals
  theme_minimal() +
  labs(title = "Relationship between Age and Splicing Efficiency at baseline", 
       x = "Scaled Age", 
       y = "Splicing Efficiency") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 10, face = "italic"),
        legend.title = element_text(size = 12, face = "bold", hjust = 0.5))




# Investigate chromosome level infoprmation
chr_RT <- RT_merged %>%
  filter(coef == "scaled_age" & time == "PreExc") %>%
  mutate(chromosome = str_extract(target, "(?<=_)[^_]+$"))

chr_RT %>%
  filter(effect != "No effect" & coef == "scaled_age") %>% ggplot(aes(chromosome))+
  geom_bar()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = -0.5) +
  ggtitle("Distribution of introns across chromosomes") +
  ylab("Number of chromosome") +
  facet_grid(~effect)



# Extract the highest ranked introns by absolute Estimate
effect_ranked <- chr_RT %>%
  arrange(fit) %>%
  slice_head(n = 100) %>%
  group_by(target) %>%
  # inner_join(RT_predictions, by = "target") %>%
  # inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version")) %>%
  dplyr::select(target, external_gene_name,  fit, Estimate,transcript_biotype, scaled_age, external_transcript_name) 

  length(unique(effect_ranked$target))


ggplot(effect_ranked, aes(x = scaled_age, y = fit, colour = target, group = target)) +
  geom_line(aes(alpha = 0.5, colour = "grey"), show.legend = F) + 
  theme_minimal()+
  scale_alpha_identity()+
  scale_color_manual(values = c("grey"), guide = "none")+
  labs(title = "Relationship between Age and Splicing Efficiency in the top 20 affected introns", 
       x = "Scaled Age", 
       y = "Splicing Efficiency") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 10, face = "italic"),
        legend.title = element_text(size = 12, face = "bold", hjust = 0.5))



length(unique(effect_ranked$external_gene_name))


# Extract those affected by RT 
RT_effect_alone <- RT_merged %>%
  filter(coef == "timePostExc")

# extract the disicnt effect types
summary_RT <- RT_effect_alone %>%
  group_by(effect) %>%
  summarize(num_targets = n_distinct(target))


RT_effect_alone <- RT_effect_alone %>%
  left_join(summary_RT, by = "effect")%>%
  mutate(effect_label = paste(effect, "(No of introns:", num_targets, ")"))

ggplot(RT_effect_alone, aes(x = scaled_age, y = fit, color = effect_label, linetype = time)) +
  geom_smooth(method = "lm", se = FALSE) + # se = FALSE to remove confidence intervals
  theme_minimal() +
  labs(title = "Relationship between RT and Splicing Efficiency", 
       x = "Scaled Age", 
       y = "Splicing Efficiency") +
  theme(plot.title = element_text(hjust = 0.5))





# Extract those with stronger effect and rank them

RT_effect_ranked <- RT_effect_alone %>%
  filter( effect != "No effect" & time == "PostExc") %>%
  arrange(fit) %>%
  slice_head(n = 200) 

RT_effect_ranked <- RT_effect_alone %>%
  filter(target %in% RT_effect_ranked$target)

ggplot(RT_effect_ranked %>%
         filter(target == "ENST00000380384.5_4_9" ), aes(x = scaled_age, y = fit, linetype = time)) +
  geom_line(aes(alpha = 0.5, colour = "grey"), show.legend = F) + 
  theme_minimal()+
  scale_alpha_identity()+
  scale_color_manual(values = c("grey"), guide = "none")+
  labs(title = "Effect of RT on top 20 affected introns", 
       x = "Scaled Age", 
       y = "Splicing Efficiency") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 10, face = "italic"),
        legend.title = element_text(size = 12, face = "bold", hjust = 0.5))

length(unique(RT_effect_df$external_transcript_name))

# hOW RT afects introns affected by agin

#get the intercept between RT and aging

#GET THOSE AFECTED BY AGINg

RT_alone <- RT_effect_alone %>%
  filter(effect != "No effect")

Aging_alone <- Aging_effect %>%
  filter(effect != "No effect")

int_aging_RT <- intersect(RT_alone$target, Aging_alone$target)



int_df <- RT_effect_alone %>%
  filter( target %in% int_aging_RT)

summary_int <- int_df %>%
  group_by(effect) %>%
  summarize(num_targets = n_distinct(target))


ggplot(int_df, aes(x = scaled_age, y = fit, color = time)) +
  geom_smooth(method = "lm", se = FALSE) + # se = FALSE to remove confidence intervals
  theme_minimal() +
  labs(title = "Relationship between RT and Splicing Efficiency in age_affected introns", 
       x = "Scaled Age", 
       y = "Splicing Efficiency") +
  geom_text(data = summary_int, aes(x = Inf, y = Inf, label = paste("Number of introns:", num_targets)), 
            hjust = 1.1, vjust = 1.1, size = 3, color = "black", inherit.aes = FALSE) +
  theme(plot.title = element_text(hjust = 0.5))+
  facet_wrap(~effect)


ggplot(int_df, aes(x = scaled_age, y = fit,  group = target, colour = time)) +
  geom_line(aes(alpha = 0.5, colour = "grey"), show.legend = F) + 
  theme_minimal()+
  labs(title = "Relationship between age and splicing efficiency", x = "Scaled Age", y = "Splicing efficiency") +
  facet_wrap(~time) +
  scale_alpha_identity()+
  scale_color_manual(values = c("grey"), guide = "none")+
  geom_text(data = summary_int, aes(x = Inf, y = Inf, label = paste("Number of introns:", num_targets)), 
            hjust = 1.1, vjust = 1.1, size = 3, color = "black", inherit.aes = FALSE) +
  theme(plot.title = element_text(hjust = 0.5))




# select only those affected by RT
RT_effect <- RT_merged %>%
  filter(effect != "No effect" ) %>%
  filter(coef == "scaled_age" )

length(unique(RT_effect$target))

# Extract the number of the various effect groups
summary_df <- RT_effect %>%
  group_by(effect) %>%
  summarize(num_targets = n_distinct(target))


ggplot(RT_effect, aes(x = scaled_age, y = fit,  group = target)) +
  geom_line(alpha = 0.5) + 
  theme_minimal()+
  labs(title = "Relationship between RT and splicing efficiency", x = "Scaled Age", y = "Splicing efficiency") +
  facet_wrap(~time) +
  scale_alpha_identity()+
  # scale_color_manual(values = c("grey"), guide = "none")+
  # geom_text(data = summary_df, aes(x = Inf, y = Inf, label = paste("Number of introns:", num_targets)), 
  #           hjust = 1.1, vjust = 1.1, size = 3, color = "black", inherit.aes = FALSE) +
  theme(plot.title = element_text(hjust = 0.5))


ggplot(RT_effect, aes(x = scaled_age, y =fit, group = coef, colour = effect )) +
  geom_smooth(method = "lm", se = FALSE) + # se = FALSE to remove confidence intervals
  theme_minimal() +
  labs(title = "Relationship between Age and Splicing Efficiency", 
       x = "Scaled Age", 
       y = "Splicing Efficiency") +
  facet_wrap(~time)
#  geom_text(data = summary_df, aes(x = max(baseline_merged$scaled_age), y = predict(lm(fit ~ scaled_age, data = baseline_merged)), label = paste("Number of introns:", num_targets)), 
#            hjust = 1, vjust = 1, size = 3, color = "black", inherit.aes = FALSE) +
  theme(plot.title = element_text(hjust = 0.5))





RT_alone <- RT_effect %>%
  filter(coef == "timePostExc") 
# RT_alone %>%
#   group_by(target) %>%
#   ggplot(aes(y= Estimate, x = odds_ratio, colour = time, group = target))+
#   geom_line() +
#   labs(title = "RT Effects Plot",
#        x = "scaled age",
#        y = "splicing efficiency") +
# #  theme_minimal() +
#   facet_wrap(~effect)
length(unique(RT_alone$target))


RT_alone_double <- RT_alone %>%
  filter(abs(log2fc) >= 0.5)
length(unique(RT_alone_double$target))

# ggplot(RT_alone , aes(x = Estimate, y = odds_ratio, colour = time)) +
#   geom_point(size = 3, alpha = 0.5) +
#   theme_minimal() +
#   labs(title = "Relationship between Fit and Scaled Age", x = "Scaled Age", y = "Fit") +
#   scale_color_manual(values = c("PreExc" = "gray20", "PostExc" = "gray70"))   # Custom color mapping
#  geom_smooth(method = "lm", col = "red")  # Adding a regression line

ggplot(RT_alone, aes(x = scaled_age, y = fit,  colour = time, group = target)) +
  geom_line(alpha = 0.1) +
  theme_minimal() +
  labs(title = "SE of introns affected by RT alone", x = "Scaled Age", y = "Splicing efficiency")+
  facet_wrap(~time)#+
# scale_color_manual(values = c("PreExc" = "black", "PostExc" = "gray70"))   # Custom color mapping
#  geom_smooth(method = "lm", col = "red")  # Adding a regression line




# bin <- RT_alone[1:20,]
# 
# write.csv(bin, "x.csv")

RT_with_age  <- RT_effect %>%
  filter(coef == "scaled_age:timePostExc")
length(unique(RT_with_age$target))

ggplot(RT_with_age, aes(x = scaled_age, y = fit, group = target, colour = time)) +
  geom_line(alpha = 0.1) +
  theme_minimal() +
  labs(title = "SE of introns affected by RT in interaction with aging", x = "Scaled Age", y = "Splicing efficiency")+
  facet_wrap(~time)

RT_age_double <- RT_with_age %>%
  filter(abs(log2fc) >= 0.5)



length(unique(RT_age_double$target))




# Load the baseline model output to query intersects between that and RT
baseline_df <- readRDS("data_new/baeline_merged_model_output.RDS") %>%
  filter(effect != "No effect")
# 

# filter targets at baseline that are in the RT model

int_df <- RT_alone %>%
  filter( target %in% baseline_df$target)


length(unique(int_df$target))

ggplot(int_df, aes(x = scaled_age, y = fit,  group = target, colour = time )) +
 
  geom_line(alpha = 0.2) + 
  theme_minimal()+
  labs(title = "Relationship between RT and sE in aging affected introns", x = "Scaled Age", y = "Splicing efficiency") +
  facet_wrap(~time) +
  scale_alpha_identity()+
 # scale_color_manual(values = c("grey"), guide = "none")+
  # geom_text(data = summary_df, aes(x = Inf, y = Inf, label = paste("Number of introns:", num_targets)), 
  #           hjust = 1.1, vjust = 1.1, size = 3, color = "black", inherit.aes = FALSE) +
  theme(plot.title = element_text(hjust = 0.5))









int_base <-  RT_effect %>%
  filter( target %in% baseline_df$target)




ggplot(int_base, aes(x = scaled_age, y = fit,  group = target, colour = time )) +
  geom_line(alpha = 0.4) + 
  theme_minimal()+
  labs(title = "Effect of RT in introns with age-related decline", x = "Scaled Age", y = "Splicing efficiency") +
  facet_wrap(~time) +
  scale_alpha_identity()+
  # scale_color_manual(values = c("grey"), guide = "none")+
  # geom_text(data = summary_df, aes(x = Inf, y = Inf, label = paste("Number of introns:", num_targets)), 
  #           hjust = 1.1, vjust = 1.1, size = 3, color = "black", inherit.aes = FALSE) +
  theme(plot.title = element_text(hjust = 0.5))


length(unique(int_base$target))



Improved_RT <- RT_effect %>%
  filter(effect == "Improved SE")


Reduced_RT <- RT_effect %>%
  filter(effect == "Reduced SE")


Inte_both <- intersect(Improved_RT$transcript_ID, Reduced_RT$transcript_ID)

# Load the batch-corrected gene expression data
gene_exp <- readRDS("data_new/gene_counts/batch_corrected_genecounts.RDS")


ref_gene <-  gene_exp[gene_exp$gene_id %in% RT_merged$ensembl_gene_id_version,] %>%
  # remove the version number to match gene_id
  mutate(gene_id = gsub("\\..*", "",  gene_id))



# Functional annotation of the genes affected
ego_df_effect <- enrichGO(gene = unique(RT_effect$ensembl_gene_id) ,
                          keyType = "ENSEMBL",
                          universe = ref_gene$gene_id,
                          OrgDb = org.Hs.eg.db, 
                          ont = "BP", 
                          pAdjustMethod = "BH", 
                          qvalueCutoff = 0.05, 
                          readable = T)



## Output results from GO analysis to a table
cluster_summary_effect <- data.frame(ego_df_effect)

dotplot(ego_df_effect,
        
        font.size = 8, title = "Enriched biological processes in genes affected RT") +
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))



# Functional annotation of the genes affected
ego_df_Improved_effect <- enrichGO(gene = unique(Improved_RT$ensembl_gene_id) ,
                          keyType = "ENSEMBL",
                          universe = ref_gene$gene_id,
                          OrgDb = org.Hs.eg.db, 
                          ont = "BP", 
                          pAdjustMethod = "BH", 
                          qvalueCutoff = 0.05, 
                          readable = T)



## Output results from GO analysis to a table
cluster_summary_Improved <- data.frame(ego_df_Improved_effect)

dotplot(ego_df_Improved_effect,
        
        font.size = 8, title = "Enriched biological processes in genes improved by RT") +
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))





# Reduced effect
# Functional annotation of the genes affected
ego_df_Reduced_effect <- enrichGO(gene = unique(Reduced_RT$ensembl_gene_id) ,
                                   keyType = "ENSEMBL",
                                   universe = ref_gene$gene_id,
                                   OrgDb = org.Hs.eg.db, 
                                   ont = "BP", 
                                   pAdjustMethod = "BH", 
                                   qvalueCutoff = 0.05, 
                                   readable = T)



## Output results from GO analysis to a table
cluster_summary_Reduced <- data.frame(ego_df_Reduced_effect)

dotplot(ego_df_Reduced_effect,
        
        font.size = 8, title = "Enriched biological processes in genes improved by RT") +
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))


