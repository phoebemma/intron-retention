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

RT_model_summary <- readRDS("data_new/processed_data/simpler_model_summary.RDS") %>%
  # select the model results of those affected by RT alone or interaction with aging
  filter(coef == "timePostExc" | coef == "scaled_age:timePostExc" | coef == "scaled_age")%>%
  # add the odds ratio for each Estimate
  mutate(odds_ratio = exp(Estimate),
         type = "model coefficient")


# load the predictions 

RT_predictions <- readRDS("data_new/processed_data/predictions_simpler_RT_model.RDS")%>%
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
         chromosome = str_extract(target, "(?<=_)[^_]+$"))

 chromosome_distribution <- chr_RT %>%
 # filter(effect != "No effect" ) %>% 
  group_by(target) %>%
  ggplot(aes(chromosome))+
  geom_bar()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.title = element_text(hjust = 0.5))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = -0.5) +
  ggtitle("Distribution of introns across chromosomes") +
  ylab("Number of chromosome") +
  facet_grid(~effect)


# Extract the most affected introns using  fit values
 doubled_aging <- Aging_effect %>%
   filter(fcthreshold == "s") %>%
   arrange(fit)%>%
   slice_head(n = 200) %>%
   dplyr::select(target, external_gene_name, time, fit, Estimate,transcript_biotype, scaled_age, external_transcript_name) 
   

   # based on the ranked list, select the disticnt ones for plotting
ranked <- Aging_effect %>%
 filter(target %in% c("ENST00000526182.1_1_1", "ENST00000258888.6_9_15",
                        "ENST00000487126.5_6_10", "ENST00000228641.4_2_12",
                      "ENST00000392676.8_34_17"))
 
 ranked_improved <- Aging_effect %>%
   filter(target %in% c( "ENST00000575116.1_1_17",
                        "ENST00000220584.9_5_8", "ENST00000252898.11_6_11"))
 


top_decline <- ggplot(ranked, aes(x = scaled_age, y = fit, colour = target)) +
  geom_line(aes(alpha = 1), show.legend = F) + 
  theme_minimal()+
  scale_alpha_identity()+
  facet_wrap(~target, scales = "free")+
  labs(title = "Some of the introns with most age-associated decline in SE", 
       x = "Scaled Age", 
       y = "Splicing Efficiency") +
    geom_text(data = ranked, aes( x = Inf, y = Inf, label =paste("Transcript name:", external_transcript_name)),
              vjust = 1.1, hjust = 1.1, size = 3, colour = "darkgrey")+
   theme(plot.title = element_text(hjust = 0.5))
 

top_improved <- ggplot(ranked_improved, aes(x = scaled_age, y = fit, colour = target)) +
  geom_line(aes(alpha = 1), show.legend = F) + 
  theme_minimal()+
  scale_alpha_identity()+
  facet_wrap(~target, scales = "free")+
  labs(title = "Some of the introns with most age-associated improvement in SE", 
       x = "Scaled Age", 
       y = "Splicing Efficiency") +
  geom_text(data = ranked_improved, aes( x = Inf, y = Inf, label =paste("Transcript name:", external_transcript_name)),
            vjust = 1.1, hjust = 1.1, size = 3, colour = "darkgrey")+
  theme(plot.title = element_text(hjust = 0.5))




# Extract the introns with Improved SE upon aging

Improved <- Aging_effect %>%
  filter(effect == "Improved SE") 

length(unique(Improved$external_gene_name))


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
  theme(axis.text = element_text(size = 10), axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 10))




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
  theme(axis.text = element_text(size = 10), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 10))










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
          top_decline, top_improved, 
          go_aging, go_aging_CC,
           ncol = 2, nrow = 3, labels = c("A", "B", "C", "D", "E", "F"),
 #         align = "hv",
#          axis = "tblr",
          label.x = 0.05,
 #         label.y = 0.05,
          font.label = list(size = 11),
          heights = c(1.3, 1, 1),
          widths = c(1, 1.2)
          )


ggsave("Figures/Figure1.png", bg = "white",scale = 2.5, dpi = 400)
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


dotplot(ego_perfect,
        
        font.size = 8, title = "Biological processes of genes containing perfectly spliced introns across all samples") +
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))


length(unique(High_SE_df$ensembl_gene_id))







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





RT_ranked <- RT_effect_alone %>%
  filter(effect != "No effect") %>%
  arrange(fit) #%>%
 # filter(fcthreshold == "s")
length(unique(RT_ranked$target))


RT_ranked_list <- RT_effect_alone %>%
  filter(target %in% c("ENST00000258888.6_9_15", "ENST00000487126.5_6_10",
                       , "ENST00000264775.9_5_5",
                        "ENST00000318023.11_12_15",
                       "ENST00000597889.1_4_19", "ENST00000682605.1_6_3", 
                       "ENST00000682605.1_6_3"))


Top_affected_introns <- ggplot(RT_ranked_list, aes(x = scaled_age, y = fit, color = time, linetype = time)) +
  geom_smooth(method = "lm", se = FALSE) + # se = FALSE to remove confidence intervals
  theme_minimal() +
  labs(title = "Relationship between RT and Splicing Efficiency in  top  affected introns", 
       x = "Scaled Age", 
       y = "Splicing Efficiency") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 10, face = "italic"),
        legend.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  facet_wrap(~target, scales = "free", ncol = 3) +
  geom_text(data = RT_ranked_list, aes( x = Inf, y = Inf, label =paste("Transcript name:", external_transcript_name)),
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












# get those affected by the interaction of aging and RT

Inter_df <- RT_merged %>%
  filter(coef == "scaled_age:timePostExc" & effect != "No effect") %>%
  arrange(fit) #%>%
#  filter(fcthreshold == "s")



inter_ranked_list <- RT_effect_alone %>%
  filter(target %in% c("ENST00000264775.9_5_5", "ENST00000380384.5_4_9",
                       "ENST00000369538.4_8_1", "ENST00000586316.5_1_19",
                       "ENST00000287777.5_2_3", "ENST00000274897.9_11_6",
                       "ENST00000674475.1_35_10", "ENST00000295314.9_8_1"))


inter_plot <- ggplot(inter_ranked_list, aes(x = scaled_age, y = fit, color = time, linetype = time)) +
  geom_smooth(method = "lm", se = FALSE) + # se = FALSE to remove confidence intervals
  theme_minimal() +
  labs(title = "Some Introns with SE associated to the interaction of aging and RT", 
       x = "Scaled Age", 
       y = "Splicing Efficiency") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 10, face = "italic"),
        legend.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  facet_wrap(~target, scales = "free", ncol = 3) +
  geom_text(data = inter_ranked_list, aes( x = Inf, y = Inf, label =paste("Transcript name:", external_transcript_name)),
            vjust = 1.1, hjust = 1.1, size = 3, colour = "darkgrey")





# Effect of RT on introns with age-associated SE

 age_rt_int <- RT_ranked[RT_ranked$target %in% Aging_df$target,] 
   
 length(unique(age_rt_int$target))
 length(unique(age_rt_int$ensembl_gene_id))
 

RT_on_age_specific_plot <-  ggplot(age_rt_int, aes(x = scaled_age, y = fit, color = effect_label, linetype = time)) +
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
