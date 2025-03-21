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

# Load the gene annotation file
gene_annotation <- readRDS("data_new/ensembl_gene_annotation.RDS")%>%
  # remove the version number to match gene_id
  mutate(gene_id = gsub("\\..*", "",  gene_id))


# The second model looks at the change in Se by age groups
group_model <- readRDS("data_new/models/grouped_model.RDS") %>%
  filter(Pr...z..<= 0.05 ) %>%
  separate("target", c("transcript_ID", "intron_ID", "chr"), sep = "_") %>%
  inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version")) %>%
  dplyr::select(coef, Estimate, transcript_ID, intron_ID,  transcript_biotype, ensembl_gene_id, external_gene_name, transcript_length)%>%
  # merge the intron_id back to the transcript ID so we track specific introns 
  mutate(intron_ID = paste0(transcript_ID, "_", intron_ID))


# Visualize the coefficients
ggplot(group_model , aes(coef, fill = coef )) +
  geom_bar()+
  ggtitle("Distribution of differentially spliced introns by age group")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1) )


hist(group_model $Estimate)




# Explore each individual age group
unique(group_model$coef)

twenties <- group_model %>%
  filter(coef == "group21 to 29")
twenties %>%
  ggplot()+
  geom_histogram(aes(x = Estimate)) +
  # ylab("Number of ds introns per Estimate")+
  # xlab("Change in SE by aging and by RT")+
  # #  ggtitle("Introns differentially spliced by both aging and RT")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))+
  theme_cowplot()


thirties <- group_model %>%
  filter(coef == "group30 to 39")

thirties %>%
  ggplot()+
  geom_histogram(aes(x = Estimate)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))+
  theme_cowplot()


fourties <- group_model %>%
  filter(coef == "group40 to 49")

fourties  %>%
  ggplot()+
  geom_histogram(aes(x = Estimate)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))+
  theme_cowplot()


fifties <- group_model %>%
  filter(coef == "group50 to 59")

fifties  %>%
  ggplot()+
  geom_histogram(aes(x = Estimate)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))+
  theme_cowplot()

sixties <- group_model %>%
  filter(coef == "group60 to 69")

sixties %>%
  ggplot()+
  geom_histogram(aes(x = Estimate)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))+
  theme_cowplot()


seventies <- group_model %>%
  filter(coef == "group70 to 79")

seventies %>%
  ggplot()+
  geom_histogram(aes(x = Estimate)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))+
  theme_cowplot()

eighties <- group_model %>%
  filter(coef == "group80 and above")



ego_twenties <- enrichGO(gene = twenties$ensembl_gene_id,
                 keyType = "ENSEMBL",
                universe = gene_counts$gene_id,
                 OrgDb = org.Hs.eg.db, 
                 ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = T)


## Output results from GO analysis to a table
summary_twenties <- data.frame(ego_twenties)

  dotplot(ego_twenties,
             
             font.size = 8, title = "Enriched bioloical processes in ds introns for those in their twenties") +
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))

  
# Eighties 
  ego_eighties <- enrichGO(gene = eighties$ensembl_gene_id,
                           keyType = "ENSEMBL",
                           universe = gene_counts$gene_id,
                           OrgDb = org.Hs.eg.db, 
                           ont = "BP", 
                           pAdjustMethod = "BH", 
                           qvalueCutoff = 0.05, 
                           readable = T)
  
  
  ## Output results from GO analysis to a table
  summary_eighties <- data.frame(ego_eighties)
  
  dotplot(ego_eighties,
          
          font.size = 8, title = "Enriched bioloical processes in ds introns for those in their twenties") +
    theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))
  
    
  
  # Sixties
  
  ego_sixties <- enrichGO(gene = sixties$ensembl_gene_id,
                           keyType = "ENSEMBL",
                           universe = gene_counts$gene_id,
                           OrgDb = org.Hs.eg.db, 
                           ont = "BP", 
                           pAdjustMethod = "BH", 
                           qvalueCutoff = 0.05, 
                           readable = T)
  
  
  ## Output results from GO analysis to a table
  summary_sixties <- data.frame(ego_sixties)
  
  dotplot(ego_sixties,
          
          font.size = 8, title = "Enriched bioloical processes in ds introns for those in their sixties") +
    theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))
  
  
  
  
  # Seventies
  ego_seventies <- enrichGO(gene = seventies$ensembl_gene_id,
                           keyType = "ENSEMBL",
                           universe = gene_counts$gene_id,
                           OrgDb = org.Hs.eg.db, 
                           ont = "BP", 
                           pAdjustMethod = "BH", 
                           qvalueCutoff = 0.05, 
                           readable = T)
  
  
  ## Output results from GO analysis to a table
  summary_seventies <- data.frame(ego_seventies)
  
  dotplot(ego_seventies,
          
          font.size = 8, title = "Enriched bioloical processes in ds introns for those in their seventies") +
    theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))
  
  
  
  # Fifties
  ego_fifties <- enrichGO(gene = fifties$ensembl_gene_id,
                           keyType = "ENSEMBL",
                           universe = gene_counts$gene_id,
                           OrgDb = org.Hs.eg.db, 
                           ont = "BP", 
                           pAdjustMethod = "BH", 
                           qvalueCutoff = 0.05, 
                           readable = T)
  
  
  ## Output results from GO analysis to a table
  summary_fifties <- data.frame(ego_fifties)
  
  dotplot(ego_fifties,
          
          font.size = 8, title = "Enriched bioloical processes in ds introns for those in their fifties") +
    theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))
  