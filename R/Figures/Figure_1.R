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


# This data contains the introns defined as differentially spliced. That is those with pvalues at or below 0.05
aging_model <- readRDS("data_new/models/filt_scaled_age_seperate_slope_intercept_model.RDS")

# Load the gene annotation data
gene_annotation <- readRDS("data_new/ensembl_gene_annotation.RDS")







# Plot the estimate distribution of the differentially spliced introns
plot_ds_introns <- aging_model %>% 
  ggplot(aes(x=Estimate))+
  geom_histogram( colour= "lightblue") +
  ylab("Number of ds introns per Estimate")+
  xlab("Estimated change in SE upon aging")
ggtitle("Model Estimate of differentially spliced introns")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))+
  theme_cowplot()





# EXtract the introns  with positvie extimates
filt_positive <- aging_model %>%
  filter(Estimate > 0) %>%
  separate("target", c("transcript_ID", "intron_ID", "chr"), sep = "_") %>%
  inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version")) %>%
  dplyr::select(Estimate, transcript_ID, intron_ID,  transcript_biotype, ensembl_gene_id, external_gene_name, transcript_length)%>%
  # merge the intron_id back to the transcript ID so we track specific introns 
  mutate(intron_ID = paste0(transcript_ID, "_", intron_ID))



# Extract those with negative estimates 
filt_negative <- aging_model %>%
  filter(Estimate < 0) %>%
  separate("target", c("transcript_ID", "intron_ID", "chr"), sep = "_") %>%
  inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version")) %>%
  dplyr::select(Estimate, transcript_ID,intron_ID,  transcript_biotype, ensembl_gene_id, external_gene_name, transcript_length)%>%
  mutate(intron_ID = paste0(transcript_ID, "_", intron_ID))


# Explore on average the intron-specificity of effects
# this counts the numbers of times each transcript_Id is recorded in the data
# This denotes the number of ds introns per transcript 
dist_com_pos <- as.data.frame(table(filt_positive$transcript_ID))

dist_pos_trans <- dist_com_pos %>%
  ggplot(aes(Freq))+
  geom_bar()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = -0.1) +
  ylab("Number of transcripts")+
  xlab("Number of introns")+
  scale_x_continuous()+
  theme_cowplot()



dist_comm_neg <- as.data.frame(table(filt_negative$transcript_ID))

dist_neg_trans <- dist_comm_neg %>%
  ggplot(aes(Freq))+
  geom_bar()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = -0.1) +
  ylab("Number of transcripts")+
  xlab("Number of introns")+
  theme_cowplot()




# Get the data above at the gene level too

dist_com_pos_gene <- as.data.frame(table(filt_positive$external_gene_name))

dist_pos_gene <- dist_com_pos_gene %>%
  ggplot(aes(Freq))+
  geom_bar()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = -0.1) +
  ylab("Number of genes")+
  xlab("Number of introns")+
  scale_x_continuous()+
  theme_cowplot()



dist_comm_neg_gene <- as.data.frame(table(filt_negative$external_gene_name))

dist_neg_gene <- dist_comm_neg_gene %>%
  ggplot(aes(Freq))+
  geom_bar()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = -0.1) +
  ylab("Number of genes")+
  xlab("Number of introns")+
  theme_cowplot()




ggarrange(plot_ds_introns,
          arrangeGrob(dist_neg_trans + rremove("xlab"), dist_pos_trans +rremove("ylab") + rremove("xlab"),dist_neg_gene, dist_pos_gene +rremove("ylab"), ncol = 2 , nrow = 2), # put both plots under the first
          
          nrow = 3,
          heights = c(2, 5, 3) )

ggsave("Images_tables/Figure2_aging_effects.png", dpi = 300, scale = 2, bg = "white")



# Gene ontology analyses

# load the batch-corrected gene counts data

gene_counts <- readRDS("data_new/gene_counts/batch_corrected_genecounts.RDS") %>%
  # remove the version number to match gene_id
  mutate(gene_id = gsub("\\..*", "",  gene_id))


# gene ontology (biological process)  analysis of the genes containing introns with improved SE with aging
ego_df_mf <- enrichGO(gene = filt_positive$ensembl_gene_id,
                      keyType = "ENSEMBL",
                      OrgDb = org.Hs.eg.db, 
                      universe = gene_counts$gene_id,
                      ont = "BP", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05, 
                      readable = T)



a<-  dotplot(ego_df_mf,
             
             font.size = 8, title = "Enriched bioloical processes in genes with positive estimate at baseline") +
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))




# for genes containing introns with SE reduced by aging
ego_df_bp <- enrichGO(gene = filt_negative$ensembl_gene_id,
                      keyType = "ENSEMBL",
                      OrgDb = org.Hs.eg.db, 
                      universe = gene_counts$gene_id,
                      ont = "BP", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05, 
                      readable = T)


b <- dotplot(ego_df_bp,
             
             font.size = 8, title = "Enriched biological processes in genes with negative estimates at baseline") +
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))




# Looking at the cellular compartments
ego_df_mf_2 <- enrichGO(gene = filt_positive$ensembl_gene_id,
                        keyType = "ENSEMBL",
                        OrgDb = org.Hs.eg.db,
                        universe = gene_counts$gene_id,
                        ont = "CC", 
                        pAdjustMethod = "BH", 
                        qvalueCutoff = 0.05, 
                        readable = T)



c<-  dotplot(ego_df_mf_2,
             
             font.size = 8, title = "Enriched cellular components in genes with positive estimate at baseline") +
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))



ego_df_bp_2 <- enrichGO(gene = filt_negative$ensembl_gene_id,
                        keyType = "ENSEMBL",
                        OrgDb = org.Hs.eg.db, 
                        universe = gene_counts$gene_id,
                        ont = "CC", 
                        pAdjustMethod = "BH", 
                        qvalueCutoff = 0.05, 
                        readable = T)


d <- dotplot(ego_df_bp_2,
             
             font.size = 8, title = "Enriched cellular components in genes with negative estimates at baseline") +
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))





ggarrange(a + rremove("xlab"), b +rremove("xlab"), c, d ,
          labels = c("A","B", "C", "D"),
          align = "hv")


ggsave("Images_tables/Figure3_introns_at_baseline.png", bg = "white" ,  scale = 3, dpi = 400)







