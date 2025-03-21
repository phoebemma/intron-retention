# This explores the baseline model outputs
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
Pre_group <- readRDS("data_new/models/scaled_age_seperate_slope_intercept_model.RDS") 



# filter based on p values
filt_pre_group <- Pre_group %>%
  filter(Pr...z..<= 0.05  )




# Visualize the EStimates
plot_ds_introns <- filt_pre_group %>% 
  ggplot(aes(x=Estimate))+
  geom_histogram( colour= "lightblue") +
  ylab("Number of ds introns per Estimate")+
  ggtitle("Model Estimate of differentially spliced introns")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_cowplot()



saveRDS(filt_pre_group, "data_new/models/filt_scaled_age_seperate_slope_intercept_model.RDS")


# Load the gene annotation file
 gene_annotation <- readRDS("data_new/ensembl_gene_annotation.RDS")


 
 # Extract the transcript_Id of the ds introns with positive Estimates
 
 filt_positive <- filt_pre_group %>%
   filter(Estimate > 0) %>%
   separate("target", c("transcript_ID", "intron_ID", "chr"), sep = "_") %>%
    inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version")) %>%
   dplyr::select(Estimate, transcript_ID, intron_ID,  transcript_biotype, ensembl_gene_id, external_gene_name, transcript_length)%>%
   # merge the intron_id back to the transcript ID so we track specific introns 
   mutate(intron_ID = paste0(transcript_ID, "_", intron_ID))
 
 
 # Filter those with negative estimates
 filt_negative <- filt_pre_group %>%
   filter(Estimate < 0) %>%
   separate("target", c("transcript_ID", "intron_ID", "chr"), sep = "_") %>%
   inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version")) %>%
   dplyr::select(Estimate, transcript_ID,intron_ID,  transcript_biotype, ensembl_gene_id, external_gene_name, transcript_length)%>%
   mutate(intron_ID = paste0(transcript_ID, "_", intron_ID))


 
 # Visualization 
filt_positive %>%
   ggplot(aes(transcript_biotype))+
   geom_bar()+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
   stat_count(geom = "Text", aes(label = ..count..), vjust = -0.5) +
  ggtitle("Annotation of genes containing introns with improved SE upon aging") +
  ylab("Number of genes")



filt_negative %>%
  ggplot(aes(transcript_biotype))+
  geom_bar()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = -0.5) +
  ggtitle("Annotation of genes containing introns with reduced SE upon aging") +
  ylab("Number of genes")




# Explore if these negativesly impacted introns shared any genes in common
length(intersect(filt_negative$transcript_ID, filt_positive$transcript_ID))

common_genes_neg <- filt_negative[filt_negative$transcript_ID %in% filt_positive$transcript_ID,]

common_genes_pos <- filt_positive[filt_positive$transcript_ID %in% filt_negative$transcript_ID,]



# check distribution of introns with negative estimates
dist_comm_neg <- as.data.frame(table(filt_negative$transcript_ID))


# Load the full pre exercise splicing data
all_pre_splice <- readRDS("data_new/Pre_Exercise/all_pre_Exc_splicing_data.RDS")

df <- all_pre_splice%>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID) )%>%
  group_by(transcript_ID) %>%
  summarize(avg = round(mean(SE), digits = 2)) %>%
  separate("transcript_ID", c("transcript_ID", "intron_ID", "chr"), sep = "_") %>%
  inner_join(annotation, by = c("transcript_ID" = "ensembl_transcript_id_version"), copy = T) %>%
  mutate(.by = external_transcript_name, 
         SE_per_gene = mean(avg))


# How many introns in a transcript of baselne data
dist <- as.data.frame(table(df$transcript_ID))


plot_all_introns <- dist %>%
  ggplot(aes(Freq))+
  geom_bar()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = -0.5) +
  ggtitle("Introns per transcript in the baseline data") +
  ylab("Number of transcripts")+
  xlab("Number of introns")



# Explore on average the intron-specificity of effects
dist_com_pos <- as.data.frame(table(filt_positive$transcript_ID))

dist_1 <- dist_com_pos %>%
  ggplot(aes(Freq))+
  geom_bar()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = -0.5) +
  ggtitle("Introns per transcript with improved SE") +
  ylab("Number of transcripts")+
  xlab("Number of introns")+
  scale_x_continuous()+
  theme_cowplot()
  



dist_2 <- dist_comm_neg %>%
  ggplot(aes(Freq))+
  geom_bar()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = -0.5) +
  ggtitle("Introns per transcript with reduced SE ") +
  ylab("Number of transcripts")+
  xlab("Number of introns")+
  theme_cowplot()



#  
ggarrange(plot_ds_introns,
          arrangeGrob(dist_2, dist_1, ncol = 2), # put both plots under the first
          nrow = 2, labels = c("A", "B"))
 
 ggsave("Images_tables/Improved_reduced_introns.png", dpi = 400, scale = 1.5, bg = "white")


# Load the batch corrected gene expression data
 
 gene_counts <- readRDS("data_new/gene_counts/batch_corrected_genecounts.RDS") %>%
   # remove the version number to match gene_id
   mutate(gene_id = gsub("\\..*", "",  gene_id))

 
 
# Functional annotation of the positive estimates 
ego_df_mf <- enrichGO(gene = filt_positive$ensembl_gene_id,
                   keyType = "ENSEMBL",
                   universe = gene_counts$gene_id,
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = T)


## Output results from GO analysis to a table
 cluster_summary_pos <- data.frame(ego_df_mf)

a<-  dotplot(ego_df_mf,
        
        font.size = 8, title = "Enriched bioloical processes in genes with positive estimate at baseline") +
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))




# Negative estimates
 ego_df_bp <- enrichGO(gene = filt_negative$ensembl_gene_id,
                      keyType = "ENSEMBL",
                      OrgDb = org.Hs.eg.db, 
                      universe = gene_counts$gene_id,
                      ont = "BP", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05, 
                      readable = T)

cluster_summary_neg <- data.frame(ego_df_bp)
b <- dotplot(ego_df_bp,
        
        font.size = 8, title = "Enriched biological processes in genes with negative estimates at baseline") +
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))




# Cellular compartments, introns with positive estimates
ego_df_mf_2 <- enrichGO(gene = filt_positive$ensembl_gene_id,
                      keyType = "ENSEMBL",
                      OrgDb = org.Hs.eg.db, 
                      ont = "CC", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05, 
                      readable = T)

## Output results from GO analysis to a table
cluster_summary_pos_2 <- data.frame(ego_df_mf_2)

c<-  dotplot(ego_df_mf_2,
             
             font.size = 8, title = "Enriched cellular components in genes with positive estimate at baseline") +
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))




# CC of the negative estimates
ego_df_bp_2 <- enrichGO(gene = filt_negative$ensembl_gene_id,
                      keyType = "ENSEMBL",
                      OrgDb = org.Hs.eg.db, 
                      ont = "CC", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05, 
                      readable = T)

cluster_summary_neg_2 <- data.frame(ego_df_bp_2)
d <- dotplot(ego_df_bp_2,
             
             font.size = 8, title = "Enriched cellular components in genes with negative estimates at baseline") +
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))





ggarrange(a, b, c, d ,
          labels = c("A","B", "C", "D"),
          align = "h")


ggsave("Images_tables/DS_introns_at_baseline.png", bg = "white" ,  scale = 2)



