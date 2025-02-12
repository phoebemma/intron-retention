library(dplyr)
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(org.Hs.eg.db)
library(clusterProfiler)



# This model analyses the interaction between  the age_groups and RT, and sex
model <- readRDS("data_new/models/full_data_RT_scaled_age_int_model.RDS") %>%
  # filter to only those with adjusted p values at or below 0.05
  filter(Pr...z..<= 0.05)

length(unique(model$target))

hist(model$Estimate)

model %>%
  ggplot(aes(x = coef)) +
  geom_bar()+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
#hist(sev_above$Estimate)
# Load the gene annotation extracted from Ensemble
gene_annotation <- readRDS("data_new/ensembl_gene_annotation.RDS")

RT_effects <- model %>%
  filter(coef == "timePostExc" | coef == "scaled_age:timePostExc") %>%
  #filter(fcthreshold == "s") %>%
  separate("target", c("transcript_ID", "intron_ID", "chr"), sep = "_") %>%
  inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version")) %>%
  dplyr::select(Estimate, coef,  transcript_ID, transcript_biotype, ensembl_gene_id, external_gene_name, transcript_length)



saveRDS(RT_effects, "data_new/processed_data/annotation_RT_effects.RDS")


# How many unique transcripts
length(unique(RT_effects$transcript_ID))

# How many unique genes
length(unique(RT_effects$ensembl_gene_id))

RT_effects %>%
  ggplot(aes(transcript_biotype))+
  geom_bar()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = -0.5) +
  ggtitle("Annotation of genes containing differentially spliced introns in full data") +
  ylab("Number of genes")




# Extract those ds by exercise
filt_exc <- RT_effects %>%
  filter(coef == "timePostExc")
hist(filt_exc$Estimate)


filt_int <- RT_effects %>%
  filter(coef == "scaled_age:timePostExc")
hist(filt_int$Estimate)

# Which introns are captured in both datasets
inter_both <-  filt_exc[filt_exc$transcript_ID %in% filt_int$transcript_ID,]

# Which introns are unique to RT effect alone
exc_not_int <- filt_exc[!filt_exc$transcript_ID %in% filt_int$transcript_ID,]


ego_df_postEXC <- enrichGO(gene = filt_exc$ensembl_gene_id,
                      keyType = "ENSEMBL",
                      OrgDb = org.Hs.eg.db, 
                      ont = "BP", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05, 
                      readable = T)

## Output results from GO analysis to a table
cluster_summary_postExc <- data.frame(ego_df_postEXC)

a <- dotplot(ego_df_postEXC,
             
             font.size = 8, title = "Enriched biological processes due to RT") +
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))





ego_df_postEXC_2 <- enrichGO(gene = filt_int$ensembl_gene_id,
                           keyType = "ENSEMBL",
                           OrgDb = org.Hs.eg.db, 
                           ont = "BP", 
                           pAdjustMethod = "BH", 
                           qvalueCutoff = 0.05, 
                           readable = T)

## Output results from GO analysis to a table
cluster_summary_postExc_2 <- data.frame(ego_df_postEXC_2)

b <- dotplot(ego_df_postEXC_2,
        
        font.size = 8, title = "Enriched biological processes based on interaction between age and RT") +
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))



ego_df_postEXC_3 <- enrichGO(gene = exc_not_int$ensembl_gene_id,
                             keyType = "ENSEMBL",
                             OrgDb = org.Hs.eg.db, 
                             ont = "BP", 
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.05, 
                             readable = T)

## Output results from GO analysis to a table
cluster_summary_postExc_3 <- data.frame(ego_df_postEXC_3)

c <- dotplot(ego_df_postEXC_3,
             
             font.size = 8, title = "Enriched biological processes based not in interaction dataset") +
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))

ggarrange(a, b,
          labels = c("A", "B"))

ggsave("Images_tables/DS_introns_at_PostExc.png", bg = "white" ,  scale = 2)




