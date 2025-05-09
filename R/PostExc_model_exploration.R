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

# Load the model summary

full_model_summary <- readRDS("data_new/simpler_model_summary.RDS")

# divid tnto two dataframes
# bind rows was used to merge the prediction and model fit. Thus leaving many NAs
# Divining itno two dataframes and merging again would solve the issue


RT_model_summary <- full_model_summary %>%
  dplyr::select(coef, target,Estimate, Pr...z.., adj.p, log2fc, fcthreshold , odds_ratio, type) %>%
  # select the model results of those affected by RT alone or interaction with aging
  filter(coef == "timePostExc" | coef == "scaled_age:timePostExc" | coef == "scaled_age")


# load the predictions 

RT_predictions <- full_model_summary %>%
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

saveRDS(RT_merged, "data_new/RT_model_outputs/RT_model_df.RDS")



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


# Load original splicing data 
all_splicing <- readRDS("data_new/processed_data/all_splice_data.RDS")
all_metadata <- readRDS("data_new/processed_data/all_full_metadata.RDS")

long_df <- all_splicing %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID) )%>%
  inner_join(all_metadata, by = "seq_sample_id")




  # based on the ranked list, select the disticnt ones for plotting
ranked <- Aging_effect %>%
 filter(target %in% c("ENST00000526182.1_1_1", "ENST00000258888.6_9_15",
                        "ENST00000504055.1_2_6", "ENST00000299432.7_6_10",
                      "ENST00000290219.11_5_21", "ENST00000359596.8_70_19",
                      "ENST00000409751.1_4_2", "ENST00000511188.2_9_4",
                      "ENST00000228641.4_2_12", "ENST00000420686.5_2_2",
                      "ENST00000375985.5_6_20", "ENST00000599111.5_10_19")) %>%
  dplyr::select(target, scaled_age, fit)
 
ranked_original  <- long_df %>%
   filter(transcript_ID %in% c("ENST00000526182.1_1_1", "ENST00000258888.6_9_15",
                               "ENST00000504055.1_2_6", "ENST00000299432.7_6_10",
                               "ENST00000290219.11_5_21", "ENST00000359596.8_70_19", 
                               "ENST00000409751.1_4_2", "ENST00000511188.2_9_4",
                               "ENST00000228641.4_2_12", "ENST00000420686.5_2_2",
                               "ENST00000375985.5_6_20", "ENST00000599111.5_10_19") &
            time == "PreExc") %>%
  group_by(scaled_age, transcript_ID) %>%
  summarize(SE = mean(SE)) %>%
  mutate(target = transcript_ID,
         fit = SE) %>%
  dplyr::select(target, scaled_age, fit)
 ranked$type <- "prediction"
 ranked_original$type <- "SE"
 
 
 ranked_df <- rbind(ranked, ranked_original) %>%
   mutate(transcript_ID = str_split(target, "_",simplify= T) [,1])%>%
   inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version"))

top_decline <- ggplot(ranked_df, aes(x = scaled_age, y = fit,  linetype = type)) +
  geom_line(aes(alpha = 1)) + 
  theme_minimal()+
  scale_alpha_identity()+
  facet_wrap(~target, scales = "free")+
  labs(title = "Some of the introns with most age-associated decline in SE", 
       x = "Scaled Age", 
       y = "Splicing Efficiency") +
      geom_text(data = ranked_df, aes( x = Inf, y = Inf, label =paste("Transcript name:", external_transcript_name)),
                vjust = 1.1, hjust = 1.1, size = 2.5, colour = "darkgrey")+
   theme(plot.title = element_text(hjust = 0.5))+
  # show only the type legend
  guides(colour = "none", linetype = guide_legend(title = "type"))

 


# select for those that appeared to have improved with aging
improved_with_aging <- Aging_effect %>%
  filter(effect == "Improved SE" & fcthreshold == "s") %>%
  arrange(fit) %>%
  dplyr::select(target, external_gene_name, time, fit, Estimate,transcript_biotype, scaled_age, external_transcript_name) 


ranked_decline <- Aging_effect %>%
  filter(target %in% c("ENST00000252992.8_4_1", "ENST00000304992.11_38_17",
                       "ENST00000242728.5_4_12", "ENST00000412401.3_1_2",
                       "ENST00000529464.5_4_8", "ENST00000471642.6_8_1",
                       "ENST00000475184.6_4_10", "ENST00000544216.8_4_19",
                       "ENST00000268220.12_3_15", "ENST00000398884.7_2_6",
                       "ENST00000505523.1_2_4", "ENST00000371083.4_2_1")) %>%
  dplyr::select(target, scaled_age, fit)
# extract those improved ones from splicing data

decline_original <- long_df %>%
  filter(transcript_ID %in% c("ENST00000252992.8_4_1", "ENST00000304992.11_38_17",
                              "ENST00000242728.5_4_12", "ENST00000412401.3_1_2",
                              "ENST00000529464.5_4_8", "ENST00000471642.6_8_1",
                              "ENST00000475184.6_4_10", "ENST00000544216.8_4_19",
                              "ENST00000268220.12_3_15", "ENST00000398884.7_2_6",
                              "ENST00000505523.1_2_4", "ENST00000371083.4_2_1")) %>%
  group_by(scaled_age, transcript_ID) %>%
  summarize(SE = mean(SE)) %>%
  mutate(target = transcript_ID,
         fit = SE) %>%
  dplyr::select(target, scaled_age, fit)

ranked_decline$type <- "prediction"
decline_original$type <- "SE"


decline_df <- rbind(ranked_decline, decline_original) %>%
  mutate(transcript_ID = str_split(target, "_",simplify= T) [,1])%>%
  inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version"))




# Access those that improved with aging
top_improved <- ggplot(decline_df, aes(x = scaled_age, y = fit,  linetype = type)) +
  geom_line() + 
  theme_minimal()+
  scale_alpha_identity()+
  facet_wrap(~target)+
  labs(title = "Some of the introns with most age-associated improvement in SE", 
       x = "Scaled Age", 
       y = "Splicing Efficiency") +
  geom_text(data = decline_df, aes( x = Inf, y = Inf, label =paste("Transcript name:", external_transcript_name)),
            vjust = 1.1, hjust = 1.1, size = 3, colour = "darkgrey")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(colour = "none", linetype = guide_legend(title = "type"))




# explore the GO of genes containing the differentailly spliced introns by aging



# Load the batch-corrected gene expression data
gene_exp_df <- readRDS("data_new/gene_counts/batch_corrected_genecounts.RDS") %>%
  # remove the version number to match gene_id
  mutate(gene_id = gsub("\\..*", "",  gene_id))


# select only genes present in the splicing data
gene_exp <- gene_exp_df[gene_exp_df$gene_id %in% RT_merged$ensembl_gene_id, ]

Aging_df <- Aging_effect  %>%
  filter(effect != "No effect")

# Functional annotation of the genes affected
ego_aging <- enrichGO(gene = unique(Aging_df$ensembl_gene_id) ,
                          keyType = "ENSEMBL",
                          universe = gene_exp$gene_id,
                          OrgDb = org.Hs.eg.db, 
                          ont = "BP", 
                          pAdjustMethod = "BH", 
                          qvalueCutoff = 0.05, 
                          readable = T)


## Output results from GO analysis to a table
cluster_aging <- data.frame(ego_aging)

go_aging <- dotplot(ego_aging,
                      
                      font.size = 8, title = "Enriched biological processes in genes containing introns with aging-associated SE") +
  theme(axis.text = element_text(size = 10), axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 10),
        plot.title = element_text(hjust = 0) )




# those not influenced 
# no_aging <- Aging_effect %>%
#   filter(effect == "No effect")


# Functional annotation of the genes affected
ego_aging_CC <- enrichGO(gene = unique(Aging_df$ensembl_gene_id) ,
                      keyType = "ENSEMBL",
                      universe = gene_exp$gene_id,
                      OrgDb = org.Hs.eg.db, 
                      ont = "CC", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05, 
                      readable = T)


cluster_aging_cc <- data.frame(ego_aging_CC)


go_aging_CC <- dotplot(ego_aging_CC,
                    
                    font.size = 8, title = "Cellular compartment of genes containing introns with aging-associated SE") +
  theme(axis.text = element_text(size = 10), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 10),
        plot.title = element_text(hjust = 0))










# subset these effects into different dataframes# scolour = # subset these effects into different dataframes# subset theseeffect effects into different dataframes
no_effect <- RT_merged %>% 
  filter(effect == "No effect")
Improved_Se <- RT_merged %>%
  filter(effect == "Improved SE")
Reduced_se <- RT_merged %>%
  filter(effect == "Reduced SE")

# Get the intersect based on transcript IDs. This is to determine if there are transcripts
# across the three groups

intersect_df <- Reduce(intersect, list(no_effect$external_gene_name, 
                                       Improved_Se$external_gene_name, Reduced_se$external_gene_name))

intersect_efect <- intersect(no_effect$transcript_ID, Reduced_se$transcript_ID)


# create a Venn list
int_data <- list("No_effect" = no_effect$external_gene_name,
                 "Improved SE" = Improved_Se$external_gene_name,
                 "Reduced SE" = Reduced_se$external_gene_name)

# Plot the UpSet plot
 upset(fromList(int_data), order.by = "freq",
                       text.scale = 1.5, 
                       mainbar.y.label = "Number of intersecting genes", 
                       sets.x.label = "Number of genes", 
                       sets.bar.color = c("grey30", "grey50", "grey60"),
                       number.angles = 45, 
                       point.size = 4.5, 
                       matrix.color = "black",
                       line.size = 2)



# plot figure 1
ggarrange(aging_plot, chromosome_distribution, 
          go_aging, go_aging_CC,
           ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"),
 #         align = "hv",
#          axis = "tblr",
          label.x = 0.05,
 #         label.y = 0.05,
          font.label = list(size = 11),
          heights = c(1.3, 1),
          widths = c(1, 1.2))


ggsave("Figures/newest_version/Figure1.png", bg = "white",scale = 2.5, dpi = 400)
# there was no obvious molecular function







# Explore perfectly spliced introns across all samples
High_SE_df <- readRDS("data_new/processed_data/perfectly_spliced_across_all.RDS") %>%
  separate("transcript_ID", c("transcript_ID", "intron_ID", "chr"), sep = "_") %>%
  inner_join(gene_annotation, by = c("transcript_ID" = "ensembl_transcript_id_version"), copy = T) 




ego_perfect <- enrichGO(gene = unique(High_SE_df$ensembl_gene_id) ,
                         keyType = "ENSEMBL",
                         universe = gene_exp_df$gene_id,
                         OrgDb = org.Hs.eg.db, 
                         ont = "BP", 
                         pAdjustMethod = "BH", 
                         qvalueCutoff = 0.05, 
                         readable = T)


cluster_perfect <- data.frame(ego_perfect)


perfect_bp <- dotplot(ego_perfect,
        
        font.size = 8, title = "Biological processes of genes containing perfectly spliced introns") +
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20),
        plot.title = element_text(hjust = 0))


ego_perfect_cc <- enrichGO(gene = unique(High_SE_df$ensembl_gene_id) ,
                        keyType = "ENSEMBL",
                        universe = gene_exp_df$gene_id,
                        OrgDb = org.Hs.eg.db, 
                        ont = "cc", 
                        pAdjustMethod = "BH", 
                        qvalueCutoff = 0.05, 
                        readable = T)


cluster_perfect_cc <- data.frame(ego_perfect_cc)


perfect_cc <- dotplot(ego_perfect_cc,
        
        font.size = 8, title = "Cellular compartment of genes containing perfectly spliced introns") +
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20),
        plot.title = element_text(hjust = 0))


ggarrange(perfect_bp, perfect_cc)







# explore the effects of resistance training

RT_effect_alone <- RT_merged %>%
  filter(coef == "timePostExc" )

# extract the disicnt effect types
summary_RT <- RT_effect_alone %>%
  group_by(effect) %>%
  summarize(num_targets = n_distinct(target))


RT_effect_alone <- RT_effect_alone %>%
  left_join(summary_RT, by = "effect")%>%
  mutate(effect_label = paste(effect, "(No of introns:", num_targets, ")"))

RT_plot <- ggplot(RT_effect_alone, aes(x = scaled_age, y = fit, color = effect_label, linetype = time)) +
  geom_smooth(method = "lm", se = FALSE) + # se = FALSE to remove confidence intervals
  theme_minimal() +
  labs(title = "Relationship between RT and Splicing Efficiency", 
       x = "Scaled Age", 
       y = "Splicing Efficiency") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 10, face = "italic"),
        legend.title = element_text(size = 12, face = "bold", hjust = 0.5))



# explore the age-dependent effects
Age_dependent_RT <- RT_merged %>%
  filter(coef == "scaled_age:timePostExc")

summary_age_rt <- Age_dependent_RT %>%
  group_by(effect) %>%
  summarize(num_targets = n_distinct(target))

Age_dependent_RT <- Age_dependent_RT %>%
  left_join(summary_age_rt, by = "effect")%>%
  mutate(effect_label = paste(effect, "(No of introns:", num_targets, ")"))

age_rt_plot <- ggplot(Age_dependent_RT, aes(x = scaled_age, y = fit, color = effect_label, linetype = time)) +
  geom_smooth(method = "lm", se = FALSE) + # se = FALSE to remove confidence intervals
  theme_minimal() +
  labs(title = "Relationship between aging, RT and Splicing Efficiency", 
       x = "Scaled Age", 
       y = "Splicing Efficiency") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 10, face = "italic"),
        legend.title = element_text(size = 12, face = "bold", hjust = 0.5))


ggarrange(RT_plot, age_rt_plot,
          labels = c("A", "B"),
          align = "hv")


RT_ranked <- RT_effect_alone %>%
  filter(effect != "No effect" & time == "PostExc") %>%
  arrange(fit) #%>%
 # dplyr::select(target, external_gene_name, time, fit, Estimate,transcript_biotype, scaled_age, external_transcript_name, transcript_length) 

 # filter(fcthreshold == "s")
length(unique(RT_ranked$target))


RT_ranked_list <- RT_effect_alone %>%
  filter(target %in% c("ENST00000504055.1_2_6", "ENST00000299432.7_6_10",
                        "ENST00000307259.9_5_5","ENST00000301364.10_13_17",
                        "ENST00000295314.9_9_1","ENST00000258888.6_9_15",
                       "ENST00000566354.1_1_16", "ENST00000553489.1_1_12", 
                       "ENST00000228641.4_2_12", "ENST00000306434.8_7_2",
                      "ENST00000432688.5_14_2", "ENST00000310417.9_12_3" ))


Top_affected_introns <- ggplot(RT_ranked_list, aes(x = scaled_age, y = fit, color = time, linetype = time)) +
  geom_smooth(method = "lm", se = FALSE) + # se = FALSE to remove confidence intervals
  theme_minimal() +
  labs(title = "Relationship between RT and Splicing Efficiency in  top  affected introns", 
       x = "Scaled Age", 
       y = "Splicing Efficiency") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 10, face = "italic"),
        legend.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  facet_wrap(~target,scale = "free",  ncol = 3) +
  geom_text(data = RT_ranked_list, aes( x = Inf, y = Inf, label =paste("Transcript name:", external_transcript_name)),
            vjust = 1.1, hjust = 1.1, size = 3, colour = "darkgrey")



# explore on gene to check for transcript specificity

VP <- RT_effect_alone %>%
  filter(external_gene_name == "VPS4A")
ggplot(VP, aes(x = scaled_age, y = fit, colour = time, linetype = time)) +
  geom_line() + 
  theme_minimal()+
  scale_alpha_identity()+
  facet_wrap(~target, scale = "free") +
  geom_text(data = VP, aes( x = Inf, y = Inf, label =external_transcript_name),
            vjust = 1.1, hjust = 1.1, size = 3, colour = "darkgrey")

# Check GO for those affected by RT


# Functional annotation of the genes affected
ego_RT <- enrichGO(gene = unique(RT_ranked$ensembl_gene_id),
                      keyType = "ENSEMBL",
                      universe = gene_exp_df$gene_id,
                      OrgDb = org.Hs.eg.db, 
                      ont = "BP", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05, 
                      readable = T)


## Output results from GO analysis to a table
cluster_RT <- data.frame(ego_RT)

go_RT <- dotplot(ego_RT,
                    
                    font.size = 8, title = "Enriched biological processes in genes containing introns with RT-associated SE") +
  theme(axis.text = element_text(size = 10), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 10))







# Functional annotation of the genes affected
ego_RT_cc <- enrichGO(gene = unique(RT_ranked$ensembl_gene_id),
                   keyType = "ENSEMBL",
                   universe = gene_exp_df$gene_id,
                   OrgDb = org.Hs.eg.db, 
                   ont = "cc", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = T)


## Output results from GO analysis to a table
cluster_RT_cc <- data.frame(ego_RT_cc)

go_RT_cc <- dotplot(ego_RT_cc,
                 
                 font.size = 8, title = "Enriched cellular compartments in genes containing introns with RT-associated SE") +
  theme(axis.text = element_text(size = 10), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 10))







# add the gene annotation to the long data
# to visualise relationship between SE and transcript length
long_df <- long_df %>%
  mutate(transcript = str_split(transcript_ID, "_",simplify= T) [,1]) %>%
  inner_join(gene_annotation, by= c("transcript" = "ensembl_transcript_id_version"))




summary_df <- long_df %>%
  group_by(transcript_ID, transcript_length) %>%
  summarize(average_SE = round(mean(SE), digits = 2))

# Plot
ggplot(summary_df, aes(x = average_SE, y = transcript_length)) +
  geom_line() + 
  theme_minimal() +
  scale_alpha_identity()


# get those affected by the interaction of aging and RT

# Inter_df <- RT_merged %>%
#   filter(coef == "scaled_age:timePostExc" & effect != "No effect") %>%
#   arrange(fit) #%>%
# #  filter(fcthreshold == "s")
# 
# 
# 
# inter_ranked_list <- RT_effect_alone %>%
#   filter(target %in% c("ENST00000264775.9_5_5", "ENST00000380384.5_4_9",
#                        "ENST00000369538.4_8_1", "ENST00000586316.5_1_19",
#                        "ENST00000287777.5_2_3", "ENST00000274897.9_11_6",
#                        "ENST00000674475.1_35_10", "ENST00000295314.9_8_1"))
# 
# 
# inter_plot <- ggplot(inter_ranked_list, aes(x = scaled_age, y = fit, color = time, linetype = time)) +
#   geom_smooth(method = "lm", se = FALSE) + # se = FALSE to remove confidence intervals
#   theme_minimal() +
#   labs(title = "Some Introns with SE associated to the interaction of aging and RT", 
#        x = "Scaled Age", 
#        y = "Splicing Efficiency") +
#   theme(plot.title = element_text(hjust = 0.5),
#         legend.text = element_text(size = 10, face = "italic"),
#         legend.title = element_text(size = 12, face = "bold", hjust = 0.5))+
#   facet_wrap(~target, scales = "free", ncol = 3) +
#   geom_text(data = inter_ranked_list, aes( x = Inf, y = Inf, label =paste("Transcript name:", external_transcript_name)),
#             vjust = 1.1, hjust = 1.1, size = 3, colour = "darkgrey")




# Effect of RT on introns with age-associated SE

RT_ranked <- RT_effect_alone %>%
  filter(effect != "No effect" )
 age_rt_int <- RT_ranked[RT_ranked$target %in% Aging_df$target,] 
   
 length(unique(age_rt_int$target))
 length(unique(age_rt_int$ensembl_gene_id))
 

RT_on_age_specific_plot <-  ggplot(age_rt_int, aes(x = scaled_age, y = fit,  linetype = time)) +
   geom_smooth(method = "lm", se = FALSE) + # se = FALSE to remove confidence intervals
   theme_minimal() +
   labs(title = "Effect of RT on introns with age-associated effects on RT", 
        x = "Scaled Age", 
        y = "Splicing Efficiency") +
   theme(plot.title = element_text(hjust = 0.5))+
   theme(plot.title = element_text(hjust = 0.5),
         legend.text = element_text(size = 10, face = "italic"),
         legend.title = element_text(size = 12, face = "bold", hjust = 0.5))
 
 age_rt_int_list <- RT_ranked %>%
   filter(target %in% c("ENST00000258888.6_9_15", "ENST00000487126.5_6_10", 
                        "ENST00000424301.6_10_5", "ENST00000536007.5_27_12", 
                        "ENST00000456057.5_9_6", "ENST00000555869.5_9_14",
                        "ENST00000502372.1_3_3","ENST00000397183.6_22_17"))
 
 
 
 Topmost_aging_Rt <- ggplot( age_rt_int_list, aes(x = scaled_age, y = fit, color = time, linetype = time)) +
   geom_smooth(method = "lm", se = FALSE) + # se = FALSE to remove confidence intervals
   theme_minimal() +
   labs(title = "Effect of RT on topmost age-affected introns", 
        x = "Scaled Age", 
        y = "Splicing Efficiency") +
   theme(plot.title = element_text(hjust = 0.5),
         legend.text = element_text(size = 10, face = "italic"),
         legend.title = element_text(size = 12, face = "bold", hjust = 0.5))+
   facet_wrap(~target, scales = "free", ncol = 3) +
   geom_text(data =  age_rt_int_list, aes( x = Inf, y = Inf, label =paste("Transcript name:", external_transcript_name)),
             vjust = 1.1, hjust = 1.1, size = 3, colour = "darkgrey")
 
 
 
 ggarrange(RT_plot, Top_affected_introns,
           RT_on_age_specific_plot, inter_plot, 
            go_RT , Topmost_aging_Rt,
           nrow = 3, ncol = 2,
           labels = c("A", "B", "C", "D", "E", "F"),
           label.x = 0.05,
           font.label = list(size = 11),
           heights = c(1.3, 1, 1),
           widths = c(1, 1.2))
 
 
 ggsave("Figures/Figure2.png", bg = "white",scale = 2.5, dpi = 400)
