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

# load baseline model summary 

baseline_summary <- readRDS("data_new/simpler_baseline_model_summary.RDS")%>%
  # filter only the scaled age
  filter(coef == "scaled_age" ) %>%
  # add the odds ratio for each Estimate
  mutate(odds_ratio = exp(Estimate),
         type = "model coefficient")


# Load the gene annotation file
gene_annotation <- readRDS("data_new/ensembl_gene_annotation.RDS")

# baseline_model predictions
baseline_predictions <- readRDS("data_new/simpler_baseline_predictions.RDS") %>%
  # Extract the transcript_id from the target
  mutate(transcript_ID = str_split(target, "_",simplify= T) [,1])


# merge both dataframes
baseline_merged <- baseline_predictions %>%
  inner_join(baseline_summary, by = "target") %>%
  drop_na() %>%
  # create an effect column that classifies the introns into those that aging has an effect on or not
  mutate(effect = case_when(Estimate > 0 & Pr...z.. <= 0.05 ~ "Improved SE",
                            Estimate < 0 & Pr...z.. <= 0.05 ~ "Reduced SE" ,
                            Estimate < 0 & Pr...z.. > 0.05 ~ "No effect",
                            Estimate > 0 & Pr...z.. > 0.05 ~ "No effect")) %>%
  inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version"))



saveRDS(baseline_merged, "data_new/baeline_merged_model_output.RDS")



# Extract the number of the various effect groups
summary_df <- baseline_merged %>%
  group_by(effect) %>%
  summarize(num_targets = n_distinct(target))


ggplot(baseline_merged, aes(x = scaled_age, y = fit,  group = target)) +
  geom_line(aes(alpha = 0.5, colour = "grey"), show.legend = F) + 
  theme_minimal()+
  labs(title = "Relationship between age and splicing efficiency", x = "Scaled Age", y = "Splicing efficiency") +
  facet_wrap(~effect) +
  scale_alpha_identity()+
  scale_color_manual(values = c("grey"), guide = "none")+
  geom_text(data = summary_df, aes(x = Inf, y = Inf, label = paste("Number of introns:", num_targets)), 
            hjust = 1.1, vjust = 1.1, size = 3, color = "black", inherit.aes = FALSE) +
  theme(plot.title = element_text(hjust = 0.5))

# subset these effects into different dataframes
no_effect <- baseline_merged %>% 
  filter(effect == "No effect")
Improved_Se <- baseline_merged %>%
  filter(effect == "Improved SE")
Reduced_se <- baseline_merged %>%
  filter(effect == "Reduced SE")

# Get the intersect based on transcript IDs. This is to determine if there are transcripts
# across the three groups

intersect_df <- Reduce(intersect, list(no_effect$transcript_ID, 
                                       Improved_Se$transcript_ID, Reduced_se$transcript_ID))

intersect_efect <- intersect(no_effect$transcript_ID, Reduced_se$transcript_ID)


# create a Venn list
int_data <- list("No_effect" = no_effect$transcript_ID,
                  "Improved_SE" = Improved_Se$transcript_ID,
                  "Reduced SE" = Reduced_se$transcript_ID)

# Plot the UpSet plot
 upset(fromList(int_data), order.by = "freq",
      text.scale = 1.5, 
      mainbar.y.label = "Number of intersecting transcripts ", 
      sets.x.label = "Number of transcripts", 
      sets.bar.color = c("grey30", "grey50", "grey60"),
      number.angles = 45, 
      point.size = 4.5, 
      matrix.color = "black",
      line.size = 2)



 
 
 
 Reduced_se%>%
   ggplot(aes(transcript_biotype))+
   geom_bar()+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
   stat_count(geom = "Text", aes(label = ..count..), vjust = -0.5) +
   ggtitle("Annotation of genes containing introns with improved SE upon aging") +
   ylab("Number of genes")
 
 
 
 # Exploring the relationship between the affected genes and those not affected
 
 
 # Load the batch-corrected gene expression data
gene_exp <- readRDS("data_new/gene_counts/batch_corrected_genecounts.RDS")
 
 
# load the baseline metadata
 pre_metadata <- readRDS("data_new/Pre_Exercise/all_prexercise_metadata.RDS")
 
 # Select the pre-exercise splicing data
 pre_intersect <- intersect(colnames(gene_exp), pre_metadata$seq_sample_id)
 
 
 # select only those in the metadata
 # this is because the full expression data contains mid exercise data
 
 pre_counts  <- gene_exp %>%
   subset(select = c("gene_id", pre_intersect))


 
 
# Explore the expression pattern of those at the intersect
 #  get the intersect at gene level
 intersect_df_gene <- Reduce(intersect, list(no_effect$ensembl_gene_id_version, 
                                        Improved_Se$ensembl_gene_id_version , Reduced_se$ensembl_gene_id_version))
 
 
 int_genes <- pre_counts[pre_counts$gene_id %in% intersect_df_gene,]
 
int_genes %>%
  pivot_longer(names_to = "seq_sample_id",
              values_to = "count",
              cols = -(gene_id))%>%
  # group_by(gene_id) %>%
  inner_join(pre_metadata, by = "seq_sample_id")%>%
  ggplot(aes(x = scaled_age, y = count)) +
  geom_point()



# does a relationship exist between those affected and those not

intersect_efect_gene <- intersect(no_effect$ensembl_gene_id_version, Reduced_se$ensembl_gene_id_version)



effect_genes <- pre_counts[pre_counts$gene_id %in% intersect_efect_gene,]

effect_genes %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "count",
               cols = -(gene_id))%>%
  # group_by(gene_id) %>%
  inner_join(pre_metadata, by = "seq_sample_id")%>%
  ggplot(aes(x = scaled_age, y = count)) +
  geom_point()




# investigate if there is a pattern in those not affected by aging
 no_effect_genes <- pre_counts[pre_counts$gene_id %in% no_effect$ensembl_gene_id_version,]
 
 
 no_effect_genes %>%
   pivot_longer(names_to = "seq_sample_id",
                values_to = "count",
                cols = -(gene_id))%>%
   # group_by(gene_id) %>%
   inner_join(pre_metadata, by = "seq_sample_id")%>%
   ggplot(aes(x = scaled_age, y = count)) +
   geom_point()
 


 
 
 
 # For gene ontology analyses, extract the genes included in our baseline model
 
 
 ref_gene <-  pre_counts[pre_counts$gene_id %in% baseline_merged$ensembl_gene_id_version,] %>%
   # remove the version number to match gene_id
   mutate(gene_id = gsub("\\..*", "",  gene_id))
 
effect_genes <- effect_genes %>%
  mutate(gene_id = gsub("\\..*", "",  gene_id))
 
 
 # Functional annotation of the genes affected
 ego_df_effect <- enrichGO(gene = effect_genes$gene_id ,
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
              
              font.size = 8, title = "Enriched biological processes in genes affected by aging at baseline") +
   theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))
 
 
 
 
 
 # Explore those not affected by aging
 
 no_effect_genes <- no_effect_genes %>%
   mutate(gene_id = gsub("\\..*", "",  gene_id))
 
 ego_df_no_effect <- enrichGO(gene = no_effect_genes$gene_id ,
                           keyType = "ENSEMBL",
                           universe = ref_gene$gene_id,
                           OrgDb = org.Hs.eg.db, 
                           ont = "BP", 
                           pAdjustMethod = "BH", 
                           qvalueCutoff = 0.05, 
                           readable = T)
 
 
 
 ## Output results from GO analysis to a table
 cluster_summary_no_effect <- data.frame(ego_df_no_effect)
 
 # dotplot(ego_df_no_effect,
 #         
 #         font.size = 8, title = "Enriched bioloical processes in genes not affected by aging at baseline") +
 #   theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))
 # 

 
 
 
 
 
 # Clear previous plots
# grid.newpage()
# 
# # Plot the Venn diagram
# venn.plot <- venn.diagram(
#   x = venn_data,
#   main = " Distribution of transcripts containing introns in the data",
#   main.cex = 1.5,
#   category.names = c("No effect", "Positive effect", "Negative effect"),
#   filename = NULL,
#   output = T,
#   cat.cex = 1.2, # Label size
#   cat.col = c("red", "blue", "green"), # Label color
#   cat.fontfamily = "serif", # Font family
#   cat.pos = c(-20, 20, 180), # Label position
#   cat.dist = c(0.05, 0.05, 0.05), # Distance of labels from circles
# 
# )
# # Display the Venn diagram
# grid.draw(venn.plot)

