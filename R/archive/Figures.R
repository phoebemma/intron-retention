# Script to reproduce the figures used in the manuscript
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(cowplot)
library(gridExtra)
library(grid)
library(gt)
library(biomaRt)
library(org.Hs.eg.db)
library(clusterProfiler)

source("R/Trainome_functions.R")

# Image of the table containing the six introns that were poorly spliced across all baseline samples
# Load all the baseline splicing data

all_pre_splice <- readRDS("data_new/Pre_Exercise/all_pre_Exc_splicing_data.RDS")

# Figure 1
# Extract the poorly spliced introns

low_SE_df <- all_pre_splice %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID)) %>%
  summarise(.by = transcript_ID, 
             mode = getmode(SE),
            min = min(SE), 
            max = max(SE), 
            mean = mean(SE),
  ) %>%
  filter(max <= 0.2)

low_SE <- low_SE_df %>%
  separate("transcript_ID", c("transcript_ID", "intron_ID", "chr"), sep = "_")

# get the saved gene annotation file
gene_annotation <- readRDS("data_new/ensembl_gene_annotation.RDS") 

# annotation of the poorly spliced introns 
annotation_low_SE <- inner_join(low_SE, gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version"))




# save the poorly spliced introns as png image
table_plot_low <- tableGrob(annotation_low_SE %>%
                              dplyr::select(transcript_ID, intron_ID, transcript_biotype,external_gene_name, mode, min, max, ))

# Open a jpeg device
jpeg("Figures/Table1.jpeg", width = 780, height = 480, quality = 100)

# Draw the table
grid.draw(table_plot_low)

# Close the png device
dev.off()






# Extract the perfectly spliced introns as baseline
# Get the perfectly spliced introns across all datasets
High_SE_df <- all_pre_splice %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID)) %>%
  summarise(.by = transcript_ID, 
            mode = getmode(SE),
            min = min(SE), 
            max = max(SE)) %>%
  filter(min == 1) 


High_SE <- High_SE_df %>%
  separate("transcript_ID", c("transcript_ID", "intron_ID", "chr"), sep = "_")

annotation_high_SE <- inner_join(High_SE, gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version"))



high_distribution <- annotation_high_SE %>%
  ggplot(aes(transcript_biotype, fill = transcript_biotype))+
  geom_bar(width = 1, alpha =.8 )+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = 1)+
  ggtitle("Biotypes of transcripts containing the perfectly spliced introns at baseline")+
  ylab("Number of unique transcripts")





# Gene ontology
# annotation of perfectly spliced introns
ego_df_high <- enrichGO(gene = annotation_high_SE$ensembl_gene_id,
                        keyType = "ENSEMBL",
                        OrgDb = org.Hs.eg.db, 
                        ont = "BP", 
                        pAdjustMethod = "BH", 
                        qvalueCutoff = 0.05, 
                        readable = T)

## Output results from GO analysis to a table
cluster_summary <- ego_df_high
high_ontology <- dotplot(ego_df_high,
                         
                         font.size = 8, title = "Enriched biological processes preexercise") +
  theme(axis.text = element_text(size = 10), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 10))





# Post exercise data

post_splice_df <- readRDS("data_new/processed_data/PostEXc_splice_data.RDS")

High_SE_full <- post_splice_df %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID)) %>%
  summarise(.by = transcript_ID, 
            mode = getmode(SE),
            min = min(SE), 
            max = max(SE)) %>%
  filter(min == 1) %>%
  separate("transcript_ID", c("transcript_ID", "intron_ID", "chr"), sep = "_")

annotation_high_SE_full <- inner_join(High_SE_full, gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version"))

# Visualize
high_post <- annotation_high_SE_full %>%
  ggplot(aes(transcript_biotype, fill = transcript_biotype))+
  geom_bar(width = 1, alpha =.8 )+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = 1)+
  ggtitle("Biotypes of transcripts containing the perfectly spliced introns at postexercise")+
  ylab("Number of unique transcripts")+
  xlab("Transcript bioypes")


ego_df_post <- enrichGO(gene = annotation_high_SE_full$ensembl_gene_id,
                        keyType = "ENSEMBL",
                        OrgDb = org.Hs.eg.db, 
                        ont = "BP", 
                        pAdjustMethod = "BH", 
                        qvalueCutoff = 0.05, 
                        readable = T)

## Output results from GO analysis to a table
cluster_summary_post <- ego_df_post
high_ontology_post <- dotplot(ego_df_post,
                         
                         font.size = 8, title = "Enriched biological processes at postexercise") +
  theme(axis.text = element_text(size = 10), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 10))





ggarrange(high_distribution + rremove("x.title"), high_ontology + rremove("x.title"), high_post, high_ontology_post,
          labels = c("A", "B", "C", "D"),
          align = "h")




ggsave("Images_tables/Figure_1_data_vis.png", dpi = 400, scale = 1.8)








 
 
