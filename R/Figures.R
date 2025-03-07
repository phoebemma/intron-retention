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







#Figure 2


# Pre exercise model outputs
#this contains the filtered outputs of the aging model
# Contains only models with p valus at or below 0.05
aging_model <- readRDS("data_new/models/filt_scaled_age_seperate_slope_intercept_model.RDS")


plot_ds_introns <- aging_model %>% 
  ggplot(aes(x=Estimate))+
  geom_histogram( colour= "lightblue") +
  ylab("Number of ds introns per Estimate")+
  xlab("Estimated change in SE upon aging")
  ggtitle("Model Estimate of differentially spliced introns")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))+
  theme_cowplot()

# EXtract those with positvie extimates
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
  

  intersect(filt_positive$transcript_ID, filt_negative$transcript_ID)
  

  # Explore on average the intron-specificity of effects
  dist_com_pos <- as.data.frame(table(filt_positive$transcript_ID))
  
  dist_1 <- dist_com_pos %>%
    ggplot(aes(Freq))+
    geom_bar()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
    stat_count(geom = "Text", aes(label = ..count..), vjust = -0.1) +
  #  ggtitle("Introns per transcript with improved SE") +
    ylab("Number of transcripts")+
    xlab("Number of introns")+
    scale_x_continuous()+
    theme_cowplot()
  


  dist_comm_neg <- as.data.frame(table(filt_negative$transcript_ID))
  
  dist_2 <- dist_comm_neg %>%
    ggplot(aes(Freq))+
    geom_bar()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
    stat_count(geom = "Text", aes(label = ..count..), vjust = -0.1) +
 #   ggtitle("Introns per transcript with reduced SE ") +
    ylab("Number of transcripts")+
    xlab("Number of introns")+
    theme_cowplot()
  
  
  
  
  # Get the data above at the gene level too
  
  dist_com_pos_gene <- as.data.frame(table(filt_positive$external_gene_name))
  
  dist_1_b <- dist_com_pos_gene %>%
    ggplot(aes(Freq))+
    geom_bar()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
    stat_count(geom = "Text", aes(label = ..count..), vjust = -0.1) +
  #  ggtitle("Introns per gene with improved SE") +
    ylab("Number of genes")+
    xlab("Number of introns")+
    scale_x_continuous()+
    theme_cowplot()
  
  
  
  dist_comm_neg_gene <- as.data.frame(table(filt_negative$external_gene_name))
  
  dist_2_b <- dist_comm_neg_gene %>%
    ggplot(aes(Freq))+
    geom_bar()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
    stat_count(geom = "Text", aes(label = ..count..), vjust = -0.1) +
   # ggtitle("Introns per gene with reduced SE ") +
    ylab("Number of genes")+
    xlab("Number of introns")+
    theme_cowplot()
  
  
  
  
  ggarrange(plot_ds_introns,
            arrangeGrob(dist_2 + rremove("xlab"), dist_1 +rremove("ylab") + rremove("xlab"),dist_2_b, dist_1_b +rremove("ylab"), ncol = 2 , nrow = 2), # put both plots under the first
            
            nrow = 3,
           heights = c(2, 5, 3) )
  
  ggsave("Images_tables/Figure2_aging_effects.png", dpi = 300, scale = 2, bg = "white")
  
  
  
 
  
  
  
  
  
  # Figure 3 
  
  ego_df_mf <- enrichGO(gene = filt_positive$ensembl_gene_id,
                        keyType = "ENSEMBL",
                        OrgDb = org.Hs.eg.db, 
                        ont = "BP", 
                        pAdjustMethod = "BH", 
                        qvalueCutoff = 0.05, 
                        readable = T)
  

  
  a<-  dotplot(ego_df_mf,
               
               font.size = 8, title = "Enriched bioloical processes in genes with positive estimate at baseline") +
    theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))
  
  
  
  ego_df_bp <- enrichGO(gene = filt_negative$ensembl_gene_id,
                        keyType = "ENSEMBL",
                        OrgDb = org.Hs.eg.db, 
                        ont = "BP", 
                        pAdjustMethod = "BH", 
                        qvalueCutoff = 0.05, 
                        readable = T)
  

  b <- dotplot(ego_df_bp,
               
               font.size = 8, title = "Enriched biological processes in genes with negative estimates at baseline") +
    theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))
  
  
  
  
  
  ego_df_mf_2 <- enrichGO(gene = filt_positive$ensembl_gene_id,
                          keyType = "ENSEMBL",
                          OrgDb = org.Hs.eg.db, 
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
  
  
  
  
  
  
  
  
  
  # Figure 4
  # Visualization of the impact of RT on Se
  RT_effects <- readRDS("data_new/processed_data/annotation_RT_effects.RDS")
  
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
  
  filt_int_neg <- filt_int %>%
    filter(Estimate < 0)
  
 resistance <-  filt_exc %>%
    ggplot(aes(x=Estimate))+
    geom_histogram( colour= "lightblue") +
    ylab("Number of ds introns per Estimate")+
    xlab("Change in SE after RT")
  ggtitle("Model Estimate of differentially spliced introns by RT")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))+
    theme_cowplot()
  
  
  # Load the full model to extract the effects of aging alone
  model <- readRDS("data_new/models/full_data_RT_scaled_age_int_model.RDS") %>%
    # filter to only those with  p values at or below 0.05
    filter(Pr...z..<= 0.05)
  
  p <- model %>%
    filter(coef ==  "scaled_age")
  
  
  o <- p%>%
    #filter(fcthreshold == "s") %>%
    separate("target", c("transcript_ID", "intron_ID", "chr"), sep = "_") %>%
    inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version")) %>%
    dplyr::select(Estimate, coef,  transcript_ID,intron_ID, transcript_biotype, ensembl_gene_id, external_gene_name, transcript_length)%>%
    # merge the intron_id back to the transcript ID so we track specific introns 
    mutate(intron_ID = paste0(transcript_ID, "_", intron_ID))
  
  saveRDS(o, "data_new/models/annotated_aging_RT_model.RDS")
  
 aging <-  o %>%
    ggplot(aes(x=Estimate))+
    geom_histogram( colour= "lightblue") +
    ylab("Number of ds introns per Estimate")+
    xlab("Change in SE due to aging")
  ggtitle("Model Estimate of differentially spliced introns by aging in full dataset")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))+
    theme_cowplot()
  
  
  
  
  # Annotation
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
  
  
  

  
  ggarrange(a, b,
            labels = c("A", "B"))
  
  ggsave("Images_tables/Figure4_annotation_RT.png", bg = "white" ,  scale = 2)
  
  

  
  
 


  y <- model %>%
    filter(coef == "scaled_age:timePostExc")
  x <- model %>%
    filter(coef == "timePostExc")
  intersect(x$target, y$target)
  
  c <- intersect(a$target, x$target)
  
  
  
  aging_in_RT <- p[p$target %in% x$target,]
  
  RT_in_aging <- x[x$target %in% p$target,]
  
 
 q <-   ggplot()+
     geom_histogram(data = aging_in_RT, aes(x = Estimate),  fill= "green") +
    geom_histogram(data = RT_in_aging, aes(x = Estimate),  fill= "darkblue")+
    ylab("Number of ds introns per Estimate")+
    xlab("Change in SE by aging and by RT")+
#  ggtitle("Introns differentially spliced by both aging and RT")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))+
    theme_cowplot()
  
  
  
 r <-  ggplot()+
    geom_point(data = aging_in_RT, aes(x = target,  y = Estimate ), colour = "green") +
    geom_point(data = RT_in_aging, aes(x = target, y = Estimate), colour = "darkblue")+
    theme(axis.text.y = element_text(size = 8))+
    xlab("DS introns by aging and RT")+
    ylab("Change in SE by aging and by RT")+
 #   ggtitle("Introns differentially spliced by both aging and RT")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))

#   
# aging_and_RT <- RT_in_aging %>%
#   inner_join( aging_in_RT, by = "target") %>%
#   dplyr::select(target, Estimate.x, Estimate.y)
# aging_and_RT %>%
#   ggplot(aes(x = Estimate.x, y = Estimate.y)) +
#   geom_point()
# 
# cor.test(aging_and_RT$Estimate.x, aging_and_RT$Estimate.y)
  
  # EXpplore those affected by aging but captured in interaction between aging and RT
  
  int_in_aging <- y[y$target %in% a$target,]
  
  aging_in_int <- a[a$target %in% y$target,]
  
 s <-  ggplot()+
    geom_point(data = int_in_aging, aes(x = target,  y = Estimate ), colour = "red") +
    geom_point(data = aging_in_int, aes(x = target, y = Estimate), colour = "purple")+
    theme(axis.text.y = element_text(size = 8))+
    xlab("DS introns by aging and interaction with RT")+
    ylab("Change in SE by aging and interaction with RT")+
  #  ggtitle("Introns differentially spliced by aging and its interaction with RT")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))
  
  
  
  
 t <-  ggplot()+
    geom_histogram(data = int_in_aging, aes(x = Estimate),  fill= "red") +
    geom_histogram(data = aging_in_int, aes(x = Estimate),  fill= "purple")+
    ylab("Number of ds introns per Estimate")+
    xlab("Change in SE by aging and interaction with RT")+
  #  ggtitle("Introns differentially spliced by aging and its interaction with RT")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))+
    theme_cowplot()
  
  
 
 ggarrange(resistance+ rremove("xlab"), aging + rremove("xlab"), q + rremove("xlab"),r+ rremove("xlab"),t+ rremove("xlab"),s+ rremove("xlab"),
          # labels = c("A", "B", "C", "D", "E", "F"),
           align = "v", ncol = 2, nrow = 3)
 
 
 
 ggsave("Images_tables/Figure5_ds_introns_aging_RT.png", bg = "white" ,  scale = 2.2, dpi = 400)
 
 
 
 
 
 
 # EXplore tyhe intersect annotation
 
 anno <- RT_in_aging %>%
   separate("target", c("transcript_ID", "intron_ID", "chr"), sep = "_") %>%
   inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version")) %>%
   dplyr::select(Estimate, coef,  transcript_ID,intron_ID, transcript_biotype, ensembl_gene_id, external_gene_name, transcript_length)%>%
   # merge the intron_id back to the transcript ID so we track specific introns 
   mutate(intron_ID = paste0(transcript_ID, "_", intron_ID))

 
 
 
 
 
