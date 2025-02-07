library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(ggpubr)
library(cowplot)
library(gridExtra)
library(grid)

# Load the primary pre-exercise model
# This was built using the formula y ~  scaled_age + (1|study) + (scaled_age+0|study) +(1|participant)
Pre_group <- readRDS("data_new/models/scaled_age_seperate_slope_intercept_model.RDS") # %>%
 # drop_na()
colnames(Pre_group)
unique(Pre_group$coef)

filt_pre_group <- Pre_group %>%
  filter(Pr...z..<= 0.05 )



filt_pre_group %>%
  ggplot(aes(x = coef)) +
  geom_bar()+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))



unique(filt_pre_group$coef)
length(unique(filt_pre_group$target))


 # Appears to be the best model
saveRDS(filt_pre_group, "data_new/models/filt_scaled_age_seperate_slope_intercept_model.RDS")



# explore the annotations of the ds introns

# Load the gene annotation file
 gene_annotation <- readRDS("data_new/ensembl_gene_annotation.RDS")



# Extract the transcript_Id of the ds introns
 
 filt <- filt_pre_group %>%
   separate("target", c("transcript_ID", "intron_ID", "chr"), sep = "_") %>%
    inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version")) %>%
   dplyr::select(Estimate, transcript_ID, transcript_biotype, ensembl_gene_id, external_gene_name, transcript_length)

 hist(filt$Estimate)
 # How many unique transcripts
length(unique(filt$transcript_ID))

# How many unique genes
length(unique(filt$ensembl_gene_id))

filt %>%
   ggplot(aes(transcript_biotype))+
   geom_bar()+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
   stat_count(geom = "Text", aes(label = ..count..), vjust = -0.5) +
  ggtitle("Annotation of genes containing differentially spliced introns at baseline") +
  ylab("Number of genes")

cor(filt$Estimate, filt$transcript_length)

cor.test(filt$Estimate, filt$transcript_length)


ego_df_mf <- enrichGO(gene = filt$ensembl_gene_id,
                   keyType = "ENSEMBL",
                   OrgDb = org.Hs.eg.db, 
                   ont = "MF", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = T)

## Output results from GO analysis to a table
 cluster_summary_mf <- data.frame(ego_df_mf)

a <- dotplot(ego_df_mf,
        
        font.size = 8, title = "Enriched molecular functions in genes containing differentially spliced introns at baseline") +
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))


ego_df_bp <- enrichGO(gene = filt$ensembl_gene_id,
                      keyType = "ENSEMBL",
                      OrgDb = org.Hs.eg.db, 
                      ont = "BP", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05, 
                      readable = T)

cluster_summary_bp <- data.frame(ego_df_bp)
b <- dotplot(ego_df_bp,
        
        font.size = 8, title = "Enriched biological functions in genes containing differentially spliced introns at baseline") +
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))




ego_df_cc <- enrichGO(gene = filt$ensembl_gene_id,
                      keyType = "ENSEMBL",
                      OrgDb = org.Hs.eg.db, 
                      ont = "CC", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05, 
                      readable = T)

cluster_summary_cc <- data.frame(ego_df_cc)
c <- dotplot(ego_df_cc,
        
        font.size = 8, title = "Enriched cellular compartments in genes containing differentially spliced introns at baseline") +
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))


ggarrange(a, b, c, )

ggsave("Figures/DS_introns_at_baseline.png", bg = "white" ,  scale = 2)



# Some of the models had positive estimates , while others had negative
# That is, those that increase with each scaled age, and those that decrease. 

# Apparently, a couple of introns improve with age

improved_SE <- filt %>%
  filter(Estimate > 0)

decreased_SE <- filt %>%
  filter(Estimate < 0)

# investigate the enriched biological functions in the enriched versus decreased 


ego_df_improved <- enrichGO(gene = improved_SE$ensembl_gene_id,
                      keyType = "ENSEMBL",
                      OrgDb = org.Hs.eg.db, 
                      ont = "BP", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05, 
                      readable = T)

cluster_summary_improved <- data.frame(ego_df_improved)
d <- dotplot(ego_df_improved,
             
             font.size = 8, title = "Enriched biological processes in genes with positive estimates at baseline") +
  theme(axis.text = element_text(size = 12), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 12))




ego_df_decreased <- enrichGO(gene = decreased_SE$ensembl_gene_id,
                            keyType = "ENSEMBL",
                            OrgDb = org.Hs.eg.db, 
                            ont = "BP", 
                            pAdjustMethod = "BH", 
                            qvalueCutoff = 0.05, 
                            readable = T)

cluster_summary_decreased <- data.frame(ego_df_decreased)
e <- dotplot(ego_df_decreased,
        
        font.size = 8, title = "Enriched biological processes in genes with negative estimates at baseline") +
  theme(axis.text = element_text(size = 12), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 12))



ego_df_improved <- enrichGO(gene = improved_SE$ensembl_gene_id,
                            keyType = "ENSEMBL",
                            OrgDb = org.Hs.eg.db, 
                            ont = "MF", 
                            pAdjustMethod = "BH", 
                            qvalueCutoff = 0.05, 
                            readable = T)

cluster_summary_improved <- data.frame(ego_df_improved)
f <- dotplot(ego_df_improved,
             
             font.size = 8, title = "Enriched molecular functions in genes with positive estimates at baseline") +
  theme(axis.text = element_text(size = 12), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 12))


ego_df_decreased <- enrichGO(gene = decreased_SE$ensembl_gene_id,
                             keyType = "ENSEMBL",
                             OrgDb = org.Hs.eg.db, 
                             ont = "MF", 
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.05, 
                             readable = T)

cluster_summary_decreased <- data.frame(ego_df_decreased)
g <- dotplot(ego_df_decreased,
             
             font.size = 8, title = "Enriched molecular functions in genes with negative estimates at baseline") +
  theme(axis.text = element_text(size = 12), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 12))







ego_df_improved <- enrichGO(gene = improved_SE$ensembl_gene_id,
                            keyType = "ENSEMBL",
                            OrgDb = org.Hs.eg.db, 
                            ont = "CC", 
                            pAdjustMethod = "BH", 
                            qvalueCutoff = 0.05, 
                            readable = T)

cluster_summary_improved <- data.frame(ego_df_improved)
h <- dotplot(ego_df_improved,
             
             font.size = 8, title = "Enriched cellular_compartments in genes with positive estimates at baseline") +
  theme(axis.text = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12))




ego_df_decreased <- enrichGO(gene = decreased_SE$ensembl_gene_id,
                             keyType = "ENSEMBL",
                             OrgDb = org.Hs.eg.db, 
                             ont = "CC", 
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.05, 
                             readable = T)

cluster_summary_decreased <- data.frame(ego_df_decreased)
i <- dotplot(ego_df_decreased,
             
             font.size = 8, title = "Enriched cellular compartments in genes with negative estimates at baseline") +
  theme(axis.text = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12))



ggarrange(d, e, f, g, h,i,
          ncol = 2,
          nrow = 3, 
          legend = "none")

ggsave("Figures/ds_introns_based_on_estimates.png", bg = "white" ,  scale = 2)

