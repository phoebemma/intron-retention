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



# extract the model summaries
mod_sum <- model_summary$summaries

# To calculate the length of each intron, we would extract one sample and generate the intron length
 # This is valid as introns included in the analysis were only those quantified across all samples
intron_length <- readr::read_tsv("data_new/Alpha_Omega_SpliceQ_outputs/A_102.tsv")
intron_length <- intron_length[!duplicated(intron_length[, 6:ncol(intron_length)]),] # removes all duplicates from the 6 column upwards


intron_length <- intron_length%>%
  mutate(transcript_ID = paste0(transcript_ID,"_", intron_ID,"_", chr),
         intron_length = abs((sj3start - sj5end) + 1) )  %>% # abs to remove effect of strand orientation
  dplyr::select(transcript_ID, intron_length)
colnames(intron_length)



# Extract the model estimate and associated variables
RT_model_summary <- mod_sum %>%
  dplyr::select(coef, target,Estimate, Pr...z..) %>%
  # select the model results of those affected by RT alone or interaction with aging
  filter(coef == "timePostExc" | coef == "scaled_age:timePostExc" | coef == "scaled_age") %>%
  mutate(adj.p = p.adjust( Pr...z.., method = "fdr"),
         log2fc = Estimate/log(2),
         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns"),
         odds_ratio = exp(Estimate)) %>%
  inner_join(intron_length, by = c("target" = "transcript_ID"))







# From the model summary, extract those that are predictions
RT_predictions <- mod_sum %>%
  dplyr::select(scaled_age, time, fit, target) %>%
  drop_na() %>%
  # Extract the transcript_id from the target
  mutate(transcript_ID = str_split(target, "_",simplify= T) [,1])







# merge dataset
RT_merged <- RT_predictions %>%
  inner_join(RT_model_summary, by = "target") %>%
  drop_na() %>%
  # based on model estimate and p values, group the introns into how they are affected by aging
  mutate(effect = case_when(Estimate > 0 & Pr...z.. <= 0.05 ~ "Improved SE", 
                            Estimate < 0 & Pr...z.. <= 0.05 ~ "Reduced SE" ,
                            Estimate < 0 & Pr...z.. > 0.05 ~ "No effect",
                            Estimate > 0 & Pr...z.. > 0.05 ~ "No effect")) %>%
  inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version"))

# length(unique(RT_merged$target))
saveRDS(RT_merged, "data/RT_model_df.RDS")






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



ggarrange(intron_distribution, gene_distribution,
          ncol = 2,
          labels = c("A", "B"))

# Explore tyhe effect of aging alone
Aging_effect <- RT_merged %>%
  filter(coef == "scaled_age" & time == "PreExc")


# extract the disicnt effect types
summary_aging <- Aging_effect %>%
  group_by(effect) %>%
  summarize(num_targets = n_distinct(target))

# get a dataset that includes a legend that displays the number of introns in each group
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


# 
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
  
# range(length_introns$intron_length)
# subset the introns with improved Se 
Imp_length <- length_introns %>%
  filter(effect == "Improved SE") 
# range(Imp_length$intron_length)

# subset those with reduced SE
Red_length <- length_introns %>%
  filter(effect == "Reduced SE") 

# range(Red_length$intron_length)
# plot the data
main_plot <- length_introns %>%
  ggplot(aes(intron_length)) +
  geom_histogram(fill = "steelblue", colour = "white")+
  ggtitle(" Distribution of intron lengths") +
  xlab("Intron length")+
  ylab("Number of introns") +
  scale_x_continuous(limits = c(70, 30000)) +
  # annotate("text", x = Inf, y = Inf, 
  #          label = paste0("Range: ", 
  #           round(min(length_introns$intron_length)), " - ", round(max(length_introns$intron_length)), " base pairs"),
  #          hjust = 1, vjust = 69, size = 4, colour = "darkgrey") +
  theme(plot.title = element_text(hjust = 0.5))


Red_plot <- Red_length %>%
  ggplot(aes(intron_length)) +
  geom_histogram(fill = "tomato", colour = "white")+
  ggtitle("Declining SE with aging ") +
  xlab(" Intron length")+
  ylab("Number of introns")+
  scale_x_continuous(limits = c(70, NA)) +
  theme(plot.title = element_text(hjust = 0.5))#+
  # annotate("text", x = Inf, y = Inf, 
  #          label = paste0("Range: ", 
  #                         round(min(Red_length$intron_length)), " - ", round(max(Red_length$intron_length)),  " base pairs"),
  #          hjust = 1.1, vjust = 29, size = 4, colour = "darkred" )




Imp_plot <- Imp_length %>%
  ggplot(aes(intron_length)) +
  geom_histogram(fill = "seagreen", colour = "white")+
  ggtitle("Improving SE with aging ") +
  xlab("Intron length")+
  ylab("Number of introns")+
  scale_x_continuous(limits = c(70, NA))  +
  theme(plot.title = element_text(hjust = 0.5)) #+
  # annotate("text", x = Inf, y = Inf, 
  #          label = paste0("Range: ", 
  #                         round(min(Imp_length$intron_length)), " - ", round(max(Imp_length$intron_length)),  " base pairs"),
  #          hjust = 1.1, vjust = 29, size = 4, colour = "darkgreen" )
  # 



# combine the three plots
length_plot <- ggdraw()+
  draw_plot(main_plot) +
  draw_plot(Red_plot, x = 0.65, y = 0.65, width = 0.3, height = 0.3) +
  draw_plot(Imp_plot , x = 0.65, y = 0.25, width = 0.3, height = 0.3)

print(length_plot)









# Plot Figure 1


# plot figure 1




# Extract the most affected introns using  fit values
doubled_aging <- Aging_effect %>%
  # select for those significantly affected by aging
  filter(fcthreshold == "s") %>%
  arrange(fit) %>%
#   slice_head(n = 200) %>%
   dplyr::select(target, external_gene_name, time, fit, Estimate,transcript_biotype, scaled_age, external_transcript_name) 
length(unique(doubled_aging$external_gene_name))   




long_df <- all_splicing %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID) )%>%
  inner_join(all_metadata, by = "seq_sample_id")




# based on the ranked list, select the disticnt ones for plotting
ranked <- Aging_effect %>%
  filter(target %in% c("ENST00000526182.1_1_1", "ENST00000343257.7_20_7",
                       "ENST00000689936.2_70_19", "ENST00000258888.6_9_15",
                       "ENST00000504055.1_2_6", "ENST00000487126.5_6_10",
                       "ENST00000290219.11_5_21", "ENST00000272167.10_7_1",
                       "ENST00000372642.5_2_9", "ENST00000649427.1_5_2",
                       "ENST00000228641.4_2_12", "ENST00000306336.6_3_2")) %>%
  dplyr::select(target, scaled_age, fit)

ranked_original  <- long_df %>%
  filter(transcript_ID %in% c("ENST00000526182.1_1_1", "ENST00000343257.7_20_7",
                              "ENST00000689936.2_70_19", "ENST00000258888.6_9_15",
                              "ENST00000504055.1_2_6", "ENST00000487126.5_6_10",
                              "ENST00000290219.11_5_21", "ENST00000272167.10_7_1",
                              "ENST00000372642.5_2_9", "ENST00000649427.1_5_2",
                              "ENST00000228641.4_2_12", "ENST00000306336.6_3_2") &
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


ranked_improved <- Aging_effect %>%
  filter(target %in% c("ENST00000369541.4_1_1", "ENST00000304992.11_38_17",
                       "ENST00000640292.2_3_1", "ENST00000242728.5_4_12",
                       "ENST00000640575.2_1_2", "ENST00000615631.5_6_8",
                       "ENST00000471642.6_8_1", "ENST00000544216.8_4_19",
                       "ENST00000609742.3_4_10", "ENST00000558134.5_3_15",
                       "ENST00000456057.5_9_6", "ENST00000650546.1_2_1")) %>%
  dplyr::select(target, scaled_age, fit)
# extract those improved ones from splicing data

improved_original <- long_df %>%
  filter(transcript_ID %in% c("ENST00000369541.4_1_1", "ENST00000304992.11_38_17",
                              "ENST00000640292.2_3_1", "ENST00000242728.5_4_12",
                              "ENST00000640575.2_1_2", "ENST00000615631.5_6_8",
                              "ENST00000471642.6_8_1", "ENST00000544216.8_4_19",
                              "ENST00000609742.3_4_10", "ENST00000558134.5_3_15",
                              "ENST00000456057.5_9_6", "ENST00000650546.1_2_1")) %>%
  group_by(scaled_age, transcript_ID) %>%
  summarize(SE = mean(SE)) %>%
  mutate(target = transcript_ID,
         fit = SE) %>%
  dplyr::select(target, scaled_age, fit)

ranked_improved$type <- "prediction"
improved_original$type <- "SE"


improved_df <- rbind(ranked_improved, improved_original) %>%
  mutate(transcript_ID = str_split(target, "_",simplify= T) [,1])%>%
  inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version"))




# Access those that improved with aging
top_improved <- ggplot(improved_df, aes(x = scaled_age, y = fit,  linetype = type)) +
  geom_line() + 
  theme_minimal()+
  scale_alpha_identity()+
  facet_wrap(~target, scales = "free")+
  labs(title = "Some of the introns with most age-associated improvement in SE", 
       x = "Scaled Age", 
       y = "Splicing Efficiency") +
  geom_text(data = improved_df, aes( x = Inf, y = Inf, label =paste("Transcript name:", external_transcript_name)),
            vjust = 1.1, hjust = 1.1, size = 3, colour = "darkgrey")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(colour = "none", linetype = guide_legend(title = "type"))




ggarrange(top_decline,NULL,  top_improved,
          nrow = 3,
          heights = c(1,0.1,1),
          labels = c("A", "", "B"),
          align = "hv")


ggsave("plots/FigureEV2.png", bg = "white",scale = 2.5, dpi = 400)

# explore the GO of genes containing the differentailly spliced introns by aging




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





ggarrange(intron_distribution, gene_distribution, aging_plot, length_plot, 
          ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"),
          align = "v",
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
ggarrange(aging_plot, length_plot, 
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















# explore the age-dependent effects
Age_dependent_RT <- RT_merged %>%
  filter(coef == "scaled_age:timePostExc")

summary_age_rt <- Age_dependent_RT %>%
  group_by(effect) %>%
  summarize(num_targets = n_distinct(target))

Age_dependent_RT <- Age_dependent_RT %>%
  left_join(summary_age_rt, by = "effect")%>%
  mutate(effect_label = paste(effect, "(No of introns:", num_targets, ")")) %>%
   filter(effect != "No effect")

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


length(unique(Age_dependent_RT$target))


# explore the effects of resistance training


RT_effect <- RT_merged %>%
  filter(coef == "timePostExc" & effect != "No effect")


# extract the disicnt effect types
summary_RT <- RT_effect %>%
  group_by(effect) %>%
  summarize(num_targets = n_distinct(target))


RT_effect <- RT_effect %>%
  left_join(summary_RT, by = "effect")%>%
  mutate(effect_label = paste(effect, "(No of introns:", num_targets, ")"))



# Extract those not in the interaction set
RT_alone_df <- RT_effect %>%
  filter(!target %in% Age_dependent_RT$target)

length(unique(RT_effect$target))

interaction_alone <- Age_dependent_RT %>%
  filter(!target %in% RT_effect$target)

length(unique(RT_alone_df$target))


length(unique(interaction_alone$target))

RT_plot  <- ggplot(RT_alone_df, aes(x = scaled_age, y = fit, color = effect_label, linetype = time)) +
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
RT_ranked <- RT_alone_df %>%
#  group_by(time) %>%
#  mutate(new_fit = fit - lag(fit)) %>%
#  ungroup()%>%
  filter(time == "PostExc") %>%
  
  arrange(fit) #%>%
 # dplyr::select(target, external_gene_name, external_transcript_name, fit, Estimate,transcript_biotype, scaled_age, external_transcript_name, transcript_length) 

# filter(fcthreshold == "s")
length(unique(RT_ranked$target))


RT_ranked_list <- RT_alone_df %>%
  # filter(target %in% c("ENST00000504055.1_2_6", "ENST00000369541.4_1_1",
  #                      "ENST00000714472.1_10_16","ENST00000301364.10_13_17",
  #                      "ENST00000258888.6_9_15","ENST00000307259.9_5_5",
  #                      "ENST00000380394.9_5_9", "ENST00000306434.8_7_2", 
  #                      "ENST00000553489.1_1_12", "ENST00000649427.1_5_2",
  #                      "ENST00000375337.4_4_9", "ENST00000168977.7_5_19"))
  filter(target %in% c("ENST00000369541.4_1_1", "ENST00000504055.1_2_6",
                       "ENST00000487126.5_6_10","ENST00000258888.6_9_15",
                       "ENST00000301364.10_13_17","ENST00000714472.1_10_16",
                       "ENST00000553489.1_1_12", "ENST00000228641.4_2_12", 
                       "ENST00000490003.5_1_3", "ENST00000395905.8_6_3",
                       "ENST00000675164.1_16_10", "ENST00000168977.7_6_19"))


Top_affected_introns <- ggplot(RT_ranked_list, aes(x = scaled_age, y = fit, color = time, linetype = time)) +
  geom_smooth(method = "lm", se = FALSE) + # se = FALSE to remove confidence intervals
  theme_minimal() +
  labs(title = "Relationship between RT and Splicing Efficiency in  12  RT_associated introns", 
       x = "Scaled Age", 
       y = "Splicing Efficiency") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 10, face = "italic"),
        legend.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  facet_wrap(~target,scale = "free",  ncol = 3) +
  geom_text(data = RT_ranked_list, aes( x = Inf, y = Inf, label =paste("Transcript name:", external_transcript_name)),
            vjust = 1.1, hjust = 1.1, size = 3, colour = "darkgrey")



# explore on gene to check for transcript specificity and intron_specificity


RT_effect_alone <- RT_merged %>%
  filter(coef == "timePostExc")

VP <- RT_effect_alone %>%
  filter(external_gene_name == "TMOD4")
ggplot(VP, aes(x = scaled_age, y = fit, colour = time, linetype = time)) +
  geom_line() + 
  theme_minimal()+
  scale_alpha_identity()+
  facet_wrap(~target + external_transcript_name, scale = "free") # +
# geom_text(data = VP, aes( x = Inf, y = Inf, label =external_transcript_name),
#           vjust = 1.1, hjust = 1.1, size = 3, colour = "darkgrey")

# load the saved upset plot
upset_plot <- image_read("Figures/Fig EV1.png")

# load the saved VP plot
int_spec_plot <- image_read("Figures/Rplot.png")
combined <- image_append(c(upset_plot, int_spec_plot)) %>%
#  image_scale("000") %>%
  image_scale("700%")

print(combined)

image_write(combined, "Figures/EV1.png")
# Check GO for those affected by RT

# Functional annotation of the genes affected
ego_RT <- enrichGO(gene = unique(RT_effect$ensembl_gene_id),
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
ego_RT_cc <- enrichGO(gene = unique(RT_effect$ensembl_gene_id),
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



ggarrange(RT_plot, age_rt_plot, go_RT, go_RT_cc,
          labels = c("A", "B", "C", "D"),
          font.label = list(size = 12),
          align = "v")


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
