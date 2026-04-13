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
library(magick)

# Load the model summary
model_summary <- readRDS("data/RT_model_summary.RDS")

# Load original splicing data 
all_splicing <- readRDS("data/all_splice.RDS") %>%
  drop_na()
all_metadata <- readRDS("data/all_full_metadata.RDS")


# Load the batch-corrected gene expression data
gene_exp_df <- readRDS("data_new/gene_counts/batch_corrected_genecounts.RDS") #%>%
# remove the version number to match gene_id
# mutate(gene_id = gsub("\\..*", "",  gene_id))



# Load the gene annotation file
gene_annotation <- readRDS("data/ensembl_gene_annotation.RDS")

# Load one file from which we will extract intron length
# This is valid as only introns quantified in all samples were included in the analyses
intron_length <- readr::read_tsv("data_new/Alpha_Omega_SpliceQ_outputs/A_102.tsv") %>%
  
  distinct(across(6:ncol(.)), .keep_all = T) %>% # Removes duplicates based on columns 6 to end
  mutate(transcript_ID = paste0(transcript_ID, "_", intron_ID, "_", chr),
         intron_length = abs((sj3start - sj5end) + 1) ) %>% # Ensures positive length regardless of strand
  dplyr::select(transcript_ID, intron_length)





# Extract the model estimate and associated variables
RT_model_summary<- model_summary$summaries %>%
  # The summaries contain both the model summaries and the predictions
  # select only the summary outputs
  dplyr::select(coef, target,Estimate, Pr...z..) %>%
  # select the model results of those affected by RT alone or interaction with aging
  filter(coef!= "(Intercept)") %>%
 # filter(coef == "timePostExc" | coef == "scaled_age:timePostExc" | coef == "scaled_age") %>%
  drop_na() %>%
  group_by(coef) %>% 
  mutate(adj.p = p.adjust( Pr...z.., method = "fdr"),
         log2fc = Estimate/log(2),
         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns"),
         odds_ratio = exp(Estimate)
       ) %>%
  ungroup() %>%
  
  
  inner_join(intron_length, by = c("target" = "transcript_ID"))




# Extract the prdictions

# From the model summary, extract those that are predictions
RT_predictions <- model_summary$summaries%>%
  dplyr::select(scaled_age, time, fit, target) %>%
  drop_na() %>%
  # Extract the transcript_id from the target
  mutate(transcript_ID = str_split(target, "_",simplify= T) [,1])



# Merge the summary and predictions into one clearner dataframe
# merge dataset
RT_merged <- RT_predictions %>%
  inner_join(RT_model_summary, by = "target") %>%
 # drop_na() %>%
  # based on model estimate and p values, group the introns into how they are affected by aging
  mutate(effect = case_when(Estimate > 0 & adj.p <= 0.05 ~ "Improved SE", 
                            Estimate < 0 & adj.p <= 0.05 ~ "Reduced SE" ,
                            Estimate < 0 & adj.p > 0.05 ~ "No effect",
                            Estimate > 0 & adj.p > 0.05 ~ "No effect")) %>%
  inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version"))

# length(unique(RT_merged$target))
#saveRDS(RT_merged, "data/RT_model_df.RDS")



# plot the distribution of gene biotypes of the genec containing the introns
intron_distribution <- RT_merged %>%
  distinct(target, transcript_biotype, .keep_all = T) %>%
  ggplot(aes(transcript_biotype, fill = transcript_biotype))+
  geom_bar()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")+ # no legend
  stat_count(geom = "Text", aes(label = ..count..), vjust = -0.1) +
  ggtitle("Distribution of biotypes of genes containing the introns") +
  ylab("Number of introns per biotype") 
  







# Plot the gene biotype to check if the distribution of biotypes also fit that of those containing introns
gene_distribution <- gene_annotation %>%
  distinct(external_gene_name, .keep_all = T) %>%
  inner_join((gene_exp_df %>%
                dplyr::select(gene_name) ), by = c("external_gene_name" = "gene_name")) %>%
#  distinct(gene_id, transcript_biotype, .keep_all = T)%>%
  ggplot(aes(transcript_biotype, fill = transcript_biotype))+
  geom_bar()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")+
  stat_count(geom = "Text", aes(label = ..count..), vjust = -0.1) +
  ggtitle("Distribution of gene biotypes in the gene expression data") +
  ylab("Number of genes per biotype")



ggarrange(intron_distribution, NULL,  gene_distribution,
          ncol = 3,
          widths = c(1,0.1,1),
          labels = c("A", "B"))

ggsave("plots/FigureEV1.png", bg = "white",scale = 2.5, dpi = 400)





# Explore the effect of aging alone and extract the number of introns per category
Aging_effect <- RT_merged %>%
  filter(coef == "scaled_age" & time == "PreExc") %>%
  group_by(effect) %>%
  mutate(num_targets = n_distinct(target)) %>%
  ungroup() %>%
  mutate(effect_label = paste(effect, "(No of introns:", num_targets, ")"))



# Visualise the plot
aging_plot <- ggplot(Aging_effect, aes(x = scaled_age, y = fit, color = effect_label)) +
  geom_smooth(method = "lm", se = FALSE) + # se = FALSE to remove confidence intervals
  theme_minimal() +
  labs(title = "Relationship between aging and Splicing Efficiency at baseline", 
       x = "Scaled Age", 
       y = "Splicing Efficiency") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 10, face = "italic"),
        legend.title = element_text(size = 10, face = "bold", hjust = 0.5))
print(aging_plot)

# # Investigate chromosome level infoprmation
# chr_RT <- RT_model_summary %>%
#   filter(coef == "scaled_age" ) %>%
#   drop_na() %>%
#   mutate(effect = case_when(Estimate > 0 & Pr...z.. <= 0.05 ~ "Improved SE",
#                             Estimate < 0 & Pr...z.. <= 0.05 ~ "Reduced SE" ,
#                             Estimate < 0 & Pr...z.. > 0.05 ~ "No effect",
#                             Estimate > 0 & Pr...z.. > 0.05 ~ "No effect"),
#          chromosome = str_extract(target, "(?<=_)[^_]+$"),
#          chromosome = factor(chromosome, levels = c("X", "1", "2", "3", "4", "5",
#                                                     "6", "7", "8", "9", "10", "11",
#                                                     "12", "13", "14", "15", "16", "17",
#                                                     "18", "19", "20", "21", "22")))
# 
# 
# chromosome_distribution <- chr_RT %>%
#   # filter(effect != "No effect" ) %>% 
#   group_by(target) %>%
#   ggplot(aes(chromosome))+
#   geom_bar()+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#         plot.title = element_text(hjust = 0.5))+
#   stat_count(geom = "Text", aes(label = ..count..), vjust = -0.5) +
#   ggtitle("Distribution of introns across chromosomes") +
#   ylab("Number of introns per chromosome") +
#   facet_grid(~effect)



# plot the length of the introns


length_introns <- Aging_effect %>%
  dplyr::select(target, intron_length, effect) %>%
  distinct()

# Create main plot
main_plot <- length_introns %>%
  ggplot(aes(x = intron_length)) +
  geom_histogram(fill = "steelblue", color = "white") +
  labs(
    title = "Distribution of Intron Lengths",
    x = "Intron Length",
    y = "Number of Introns"
  ) +
  scale_x_continuous(limits = c(70, 30000)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Create inset plot for Reduced SE
Red_plot <- length_introns %>%
  filter(effect == "Reduced SE") %>%
  ggplot(aes(x = intron_length)) +
  geom_histogram(fill = "tomato", color = "white") +
  labs(title = "Declining SE with Aging", x = "Intron Length", y = "Number of Introns") +
  scale_x_continuous(limits = c(70, NA)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Create inset plot for Improved SE
Imp_plot <- length_introns %>%
  filter(effect == "Improved SE") %>%
  ggplot(aes(x = intron_length)) +
  geom_histogram(fill = "seagreen", color = "white") +
  labs(title = "Improving SE with Aging", x = "Intron Length", y = "Number of Introns") +
  scale_x_continuous(limits = c(70, NA)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Combine plots
length_plot <- ggdraw() +
  draw_plot(main_plot) +
  draw_plot(Red_plot, x = 0.65, y = 0.65, width = 0.3, height = 0.3) +
  draw_plot(Imp_plot, x = 0.65, y = 0.25, width = 0.3, height = 0.3)

# Display
print(length_plot)






# Extract the most affected introns using  fit values
doubled_aging <- Aging_effect %>%
  filter(fcthreshold == "s") %>%
  arrange(fit) %>%
  dplyr::select(target, external_gene_name, time, fit, Estimate,adj.p, transcript_biotype, scaled_age, external_transcript_name)


# 2. Reshape splicing data and join with metadata
long_df <- all_splicing %>%
  pivot_longer(cols = -transcript_ID, names_to = "seq_sample_id", values_to = "SE") %>%
  inner_join(all_metadata, by = "seq_sample_id")

# 3. Select top introns for plotting
top_introns <- c(
  "ENST00000526182.1_1_1", "ENST00000343257.7_20_7", "ENST00000689936.2_70_19",
  "ENST00000258888.6_9_15", "ENST00000504055.1_2_6", "ENST00000487126.5_6_10",
  "ENST00000290219.11_5_21", "ENST00000272167.10_7_1", "ENST00000372642.5_2_9",
  "ENST00000649427.1_5_2", "ENST00000228641.4_2_12", "ENST00000306336.6_3_2"
)

# Predicted fit values
ranked <- Aging_effect %>%
  filter(target %in% top_introns) %>%
  dplyr::select(target, scaled_age, fit) %>%
  mutate(type = "prediction")

# Observed SE values
ranked_original <- long_df %>%
  filter(transcript_ID %in% top_introns, time == "PreExc") %>%
  group_by(scaled_age, transcript_ID) %>%
  summarize(SE = mean(SE), .groups = "drop") %>%
  mutate(target = transcript_ID, fit = SE, type = "SE") %>%
  dplyr::select(target, scaled_age, fit, type)

# 4. Combine and annotate
ranked_df <- bind_rows(ranked, ranked_original) %>%
  mutate(transcript_ID = str_split(target, "_", simplify = TRUE)[, 1]) %>%
  inner_join(gene_annotation, by = c("transcript_ID" = "ensembl_transcript_id_version"))

# 5. Plot
top_decline <- ggplot(ranked_df, aes(x = scaled_age, y = fit, shape = type, colour = type)) +
  geom_point(aes(alpha = 1)) +
  facet_wrap(~target, scales = "free") +
  theme_minimal() +
  scale_alpha_identity() +
  labs(
    title = "Introns with most age-associated decline in SE",
    x = "Scaled Age",
    y = "Splicing Efficiency"
  ) +
  geom_text(
    aes(x = Inf, y = Inf, label = paste("Transcript name:", external_transcript_name)),
    vjust = 1.1, hjust = 1.1, size = 2.5, colour = "darkgrey"
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides( shape = guide_legend(title = "type"))

print(top_decline)






# select for those that appeared to have improved with aging
improved_with_aging <- Aging_effect %>%
  filter(effect == "Improved SE", fcthreshold == "s") %>%
  arrange(fit) %>%
  dplyr::select(target, external_gene_name, time, fit, Estimate,adj.p, transcript_biotype, scaled_age, external_transcript_name)

# 2. Define top improved introns
top_improved_ids <- c( "ENST00000650546.1_2_1")

# 3. Extract predicted fit values
ranked_improved <- Aging_effect %>%
  filter(target %in% top_improved_ids) %>%
  dplyr::select(target, scaled_age, fit) %>%
  mutate(type = "prediction")

# 4. Extract observed SE values
improved_original <- long_df %>%
  filter(transcript_ID %in% top_improved_ids) %>%
  group_by(scaled_age, transcript_ID) %>%
  summarize(SE = mean(SE), .groups = "drop") %>%
  mutate(target = transcript_ID, fit = SE, type = "SE") %>%
  dplyr::select(target, scaled_age, fit, type)

# 5. Combine and annotate
improved_df <- bind_rows(ranked_improved, improved_original) %>%
  mutate(transcript_ID = str_split(target, "_", simplify = TRUE)[, 1]) %>%
  inner_join(gene_annotation, by = c("transcript_ID" = "ensembl_transcript_id_version"))

# 6. Plot SE improvement with aging
top_improved <- ggplot(improved_df, aes(x = scaled_age, y = fit, shape = type, color = type)) +
  geom_point(alpha = 1) +
  # geom_smooth(method = "lm", se = FALSE, aes(group = type), linetype = "dashed", size = 0.6) +
  facet_wrap(~target, scales = "free") +
  theme_minimal() +
  scale_color_manual(values = c("prediction" = "blue", "SE" = "green")) +
  labs(
    title = "Intron with age-associated Improvement in Splicing Efficiency",
    x = "Scaled Age",
    y = "Splicing Efficiency"
  ) +
  geom_text(
    aes(x = Inf, y = Inf, label = paste("Transcript name:", external_transcript_name)),
    vjust = 1.1, hjust = 1.1, size = 3, colour = "darkgrey"
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(shape = guide_legend(title = "Type"), color = guide_legend(title = "Type"))




ggarrange(top_decline,NULL,  top_improved,
          nrow = 3,
          heights = c(1,0.1,1),
          labels = c("A", "", "B"),
          align = "hv")


ggsave("plots/FigureEV2.png", bg = "white",scale = 2.5, dpi = 400)

# explore the GO of genes containing the differentailly spliced introns by aging




# select only genes present in the splicing data
gene_exp <- gene_exp_df[gene_exp_df$gene_name %in% RT_merged$external_gene_name, ]

Aging_df <- Aging_effect  %>%
  filter(effect != "No effect")

saveRDS(Aging_df, "data/Aging_affected_introns.RDS")

# Functional annotation of the genes affected
ego_aging <- enrichGO(gene = unique(Aging_df$external_gene_name) ,
                      keyType = "SYMBOL",
                      universe = gene_exp$gene_name,
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
ego_aging_CC <- enrichGO(gene = unique(Aging_df$external_gene_name) ,
                         keyType = "SYMBOL",
                         universe = gene_exp$gene_name,
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




# explore those not affected by aging

no_aging_effect <- Aging_effect  %>%
  filter(effect == "No effect")


ego_no_aging <- enrichGO(gene = unique(no_aging_effect$external_gene_name) ,
                      keyType = "SYMBOL",
                      universe = gene_exp$gene_name,
                      OrgDb = org.Hs.eg.db, 
                      ont = "BP", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05, 
                      readable = T)



cluster_no_aging <- data.frame(ego_no_aging)

go_no_aging <- dotplot(ego_no_aging,
                    
                    font.size = 8, title = "Enriched biological processes in genes containing introns not affected by aging") +
  theme(axis.text = element_text(size = 10), axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 10),
        plot.title = element_text(hjust = 0) )




ggarrange(aging_plot, length_plot, go_aging,
          ncol = 2, nrow = 2, labels = c("A", "B", "C"),
         #align = "v",
          #          axis = "tblr",
          label.x = 0.05,
          #         label.y = 0.05,
          font.label = list(size = 11),
          heights = c(1.3, 1),
          widths = c(1, 1.2))


ggsave("plots/Figure1.png", bg = "white",scale = 2.5, dpi = 400)
# there was no obvious molecular function




# subset these effects into different dataframes
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
# ggarrange(aging_plot, length_plot, 
#           go_aging, go_aging_CC,
#           ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"),
#           #         align = "hv",
#           #          axis = "tblr",
#           label.x = 0.05,
#           #         label.y = 0.05,
#           font.label = list(size = 11),
#           heights = c(1.3, 1),
#           widths = c(1, 1.2))
# 
# 
# ggsave("Figures/newest_version/Figure1.png", bg = "white",scale = 2.5, dpi = 400)
# # there was no obvious molecular function











# Explore the effect of aging and exercisealone and extract the number of introns per category
Age_dependent_RT <- RT_merged %>%
  filter(coef == "scaled_age:timePostExc") %>%
  group_by(effect) %>%
  mutate(num_targets = n_distinct(target)) %>%
  ungroup() %>%
  mutate(effect_label = paste(effect, "(No of introns:", num_targets, ")"))


# Visualise the plot
age_rt_plot <- ggplot(Age_dependent_RT, aes(x = scaled_age, y = fit, color = effect_label, linetype = time)) +
  geom_smooth(method = "lm", se = FALSE) + # se = FALSE to remove confidence intervals
  theme_minimal() +
  labs(title = "Age-dependent relationship between RT and Splicing Efficiency", 
       x = "Scaled Age", 
       y = "Splicing Efficiency") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 10, face = "italic"),
        legend.title = element_text(size = 12, face = "bold", hjust = 0.5))


# explore the effects of resistance training


# extract only the age-dpendent effects. IE exclude those without effect
Age_dep_RT <- Age_dependent_RT %>%
  filter(effect != "No effect")

saveRDS(Age_dep_RT, "data/RT_and_age_effect.RDS")
# Explore the effect of aging alone and extract the number of introns per category
RT_effect <- RT_merged %>%
  filter(coef == "timePostExc" & effect != "No effect") %>%
  group_by(effect) %>%
  mutate(num_targets = n_distinct(target)) %>%
  ungroup() %>%
  mutate(effect_label = paste(effect, "(No of introns:", num_targets, ")")) #%>%
  # To make it exclusive to those not affected by an interaction with age,
  # We exclude those in the age-dependent RT effect group
  #filter(!target %in% Age_dep_RT$target)


saveRDS(RT_effect, "data/RT_alone_effect.RDS")


RT_plot  <- ggplot(RT_effect, aes(x = scaled_age, y = fit, color = effect_label, linetype = time)) +
  geom_smooth(method = "lm", se = FALSE) + # se = FALSE to remove confidence intervals
  theme_minimal() +
  labs(title = "Relationship between  RT and Splicing Efficiency", 
       x = "Scaled Age", 
       y = "Splicing Efficiency") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 10, face = "italic"),
        legend.title = element_text(size = 12, face = "bold", hjust = 0.5))




# How did RT affect the aging-affected introns

age_variables <- Aging_df %>%
  dplyr::select(target, effect, effect_label) %>%
  dplyr::rename(age_effect = effect)


Rt_on_age <- RT_effect %>%
  filter(target %in% Aging_df$target) %>%
  # drop the effect label %>%
  dplyr::select(-effect_label) %>%
  inner_join(age_variables, by = "target") %>%
  # extract its own effect label
  group_by(effect) %>%
  mutate(num_targets = n_distinct(target)) %>%
  ungroup() %>%
  mutate(effect_label = paste(effect, "(No of introns:", num_targets, ")")) 
  

  ggplot(Rt_on_age, aes(x = scaled_age, y = fit, color = effect_label, linetype = time)) +
  geom_smooth(method = "lm", se = FALSE) + # se = FALSE to remove confidence intervals
  theme_minimal() +
  labs(title = "Relationship between  RT and Splicing Efficiency", 
       x = "Scaled Age", 
       y = "Splicing Efficiency") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 10, face = "italic"),
        legend.title = element_text(size = 12, face = "bold", hjust = 0.5))
  










# Get a ranking of the RT introns
  
RT_ranked <- RT_effect %>%
  filter(time == "PostExc" , effect == "Improved SE", fcthreshold == "s") %>%
  arrange(fit) 

RT_ranked_list <- RT_effect %>%
  filter(target %in% c("ENST00000401702.5_2_22"
    # "ENST00000369541.4_1_1", "ENST00000504055.1_2_6",
    #                    "ENST00000487126.5_6_10","ENST00000258888.6_9_15",
    #                    "ENST00000301364.10_13_17","ENST00000714472.1_10_16",
    #                    "ENST00000553489.1_1_12", "ENST00000228641.4_2_12", 
    #                    "ENST00000490003.5_1_3", "ENST00000395905.8_6_3",
    #                    "ENST00000675164.1_16_10", "ENST00000168977.7_6_19"
    ))


Top_affected_introns <- ggplot(RT_ranked_list, aes(x = scaled_age, y = fit, color = time, linetype = time)) +
  geom_smooth(method = "lm", se = FALSE) + # se = FALSE to remove confidence intervals
  theme_minimal() +
  labs(title = "Relationship between RT and Splicing Efficiency in  12  RT_associated introns", 
       x = "Scaled Age", 
       y = "Splicing Efficiency") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 10, face = "italic"),
        legend.title = element_text(size = 12, face = "bold", hjust = 0.5))+
 # facet_wrap(~target,scale = "free",  ncol = 3) +
  geom_text(data = RT_ranked_list, aes( x = Inf, y = Inf, label =paste("Transcript name:", external_transcript_name)),
            vjust = 1.1, hjust = 1.1, size = 3, colour = "darkgrey")
# explore on gene to check for transcript specificity and intron_specificity


# Functional annotation of the genes affected
ego_RT <- enrichGO(gene = unique(RT_effect$external_gene_name),
                   keyType = "SYMBOL",
                   universe = gene_exp_df$gene_name,
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = T)



go_RT <- dotplot(ego_RT,
                 
                 font.size = 8, title = "Enriched biological processes in genes containing introns with RT-associated SE") +
  theme(axis.text = element_text(size = 10), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 10))
go_RT




# Functional annotation of the genes affected
ego_age_RT <- enrichGO(gene = unique(Age_dep_RT$external_gene_name),
                       keyType = "SYMBOL",
                       universe = gene_exp_df$gene_name,
                       OrgDb = org.Hs.eg.db, 
                       ont = "BP", 
                       pAdjustMethod = "BH", 
                       qvalueCutoff = 0.05, 
                       readable = T)



go_age_RT <- dotplot(ego_age_RT,
                     
                     font.size = 8, title = "Enriched biological processes in genes containing introns with age-dependent RT associations on SE") +
  theme(axis.text = element_text(size = 10), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 10))







# Functional annotation of the genes affected
ego_RT_cc <- enrichGO(gene = unique(RT_effect$ensembl_gene_id),
                      keyType = "ENSEMBL",
                      universe = gene_exp_df$gene_id,
                      OrgDb = org.Hs.eg.db, 
                      ont = "cc", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05, 
                      readable = T)



go_RT_cc <- dotplot(ego_RT_cc,
                    
                    font.size = 8, title = "Enriched cellular compartments in genes containing introns with RT-associated SE") +
  theme(axis.text = element_text(size = 10), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 10))


# check if there are ontological differences between genes containing the reduced and imrpoved introns

# # Get a ranking of the RT introns
# RT_improved <- RT_effect %>%
#   filter(effect == "Improved SE" )
# 
# 
# # Functional annotation of the genes affected
# ego_RT_reduced <- enrichGO(gene = unique(RT_improved$ensembl_gene_id),
#                    keyType = "ENSEMBL",
#                    universe = gene_exp_df$gene_id,
#                    OrgDb = org.Hs.eg.db, 
#                    ont = "BP", 
#                    pAdjustMethod = "BH", 
#                    qvalueCutoff = 0.05, 
#                    readable = T)
# 
# 
# ## Output results from GO analysis to a table
# cluster_RT_reduced <- data.frame(ego_RT_reduced)
# 
# dotplot(ego_RT_reduced,
#         
#         font.size = 8, title = "Enriched biological processes in genes containing introns with RT-associated decline in SE") +
#   theme(axis.text = element_text(size = 10), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 10))





# extract introns affected by aging, and by RT
# 
# age_rt_int <- RT_ranked[RT_ranked$target %in% Aging_df$target,] 
# 
# ranked_age_rt <- age_rt_int%>%
#   group_by(time) %>%
#   mutate(new_fit = fit - lag(fit)) %>%
#   ungroup()%>%
#   filter(time == "PostExc") %>%
#   arrange(new_fit) %>%
#   dplyr::select(target, external_gene_name,  new_fit, Estimate,transcript_biotype, scaled_age, external_transcript_name, transcript_length) 
# 
# 
# length(unique(age_rt_int$target))
# length(unique(age_rt_int$ensembl_gene_id))
# 
# 
# RT_on_age_specific_plot <-  ggplot(age_rt_int, aes(x = scaled_age, y = fit,  linetype = time)) +
#   geom_smooth(method = "lm", se = FALSE) + # se = FALSE to remove confidence intervals
#   theme_minimal() +
#   labs(title = "Effect of RT on introns with age-associated effects on RT", 
#        x = "Scaled Age", 
#        y = "Splicing Efficiency") +
#   theme(plot.title = element_text(hjust = 0.5))+
#   theme(plot.title = element_text(hjust = 0.5),
#         legend.text = element_text(size = 10, face = "italic"),
#         legend.title = element_text(size = 12, face = "bold", hjust = 0.5))
# 
# age_rt_int_list <- RT_ranked %>%
#   # filter(target %in% c("ENST00000504055.1_2_6", "ENST00000409751.1_4_2", 
#   #                      "ENST00000258888.6_9_15", "ENST00000553294.1_2_12", 
#   #                      "ENST00000228641.4_2_12", "ENST00000487126.5_5_10",
#   #                      "ENST00000567071.5_1_15","ENST00000437669.5_9_5",
#   #                      "ENST00000361466.7_30_12", "ENST00000640575.2_1_2",
#   #                     "ENST00000268661.8_7_16", "ENST00000262113.9_14_8" ))
#   filter(target %in% c("ENST00000369541.4_1_1", "ENST00000504055.1_2_6",
#                        "ENST00000258888.6_9_15", "ENST00000487126.5_6_10",
#                        "ENST00000649427.1_5_2", "ENST00000483208.5_16_15",
#                        "ENST00000553489.1_1_12","ENST00000228641.4_2_12",
#                        "ENST00000547405.5_26_12", "ENST00000508045.5_8_17",
#                        "ENST00000228641.4_1_12", "ENST00000262113.9_14_8" ))
# 
# 
# 
# Topmost_aging_Rt <- ggplot( age_rt_int_list, aes(x = scaled_age, y = fit, color = time, linetype = time)) +
#   geom_smooth(method = "lm", se = FALSE) + # se = FALSE to remove confidence intervals
#   theme_minimal() +
#   labs(title = "Effect of RT on 9 topmost age-affected introns", 
#        x = "Scaled Age", 
#        y = "Splicing Efficiency") +
#   theme(plot.title = element_text(hjust = 0.5),
#         legend.text = element_text(size = 10, face = "italic"),
#         legend.title = element_text(size = 12, face = "bold", hjust = 0.5))+
#   facet_wrap(~target, scales = "free", ncol = 3) +
#   geom_text(data =  age_rt_int_list, aes( x = Inf, y = Inf, label =paste("Transcript name:", external_transcript_name)),
#             vjust = 1.1, hjust = 1.1, size = 3, colour = "darkgrey")
# 
# 
# 
# ggarrange(RT_plot, Top_affected_introns,
#           RT_on_age_specific_plot, inter_plot, 
#           go_RT , Topmost_aging_Rt,
#           nrow = 3, ncol = 2,
#           labels = c("A", "B", "C", "D", "E", "F"),
#           label.x = 0.05,
#           font.label = list(size = 11),
#           heights = c(1.3, 1, 1),
#           widths = c(1, 1.2))
# 
