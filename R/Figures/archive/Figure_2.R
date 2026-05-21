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



# This presents the baseline results of the participants when grouped in age groups

# Visualize the distribution of participants in their age groups

# Load the pre-exercise metadata
# Load metadata
all_pre_metadata <- readRDS("data_new/Pre_Exercise/all_prexercise_metadata.RDS") %>%
  mutate(group = case_when(age <=20 ~ "20 and below" ,
   age > 20 & age < 30 ~ "21 to 29",
   age >= 30 & age < 40 ~ "30 to 39", 
   age >= 40 & age < 50 ~ "40 to 49",
  age >= 50 & age < 60 ~ "50 to 59",
  age >= 60 & age < 70 ~ "60 to 69",
  age >= 70 & age < 80 ~ "70 to 79",
 age >= 80 ~ "80 and above")) %>%
  mutate(group = factor(group, levels = c("20 and below", "21 to 29", "30 to 39",
                                          "40 to 49", "50 to 59",  "60 to 69",
                                          "70 to 79", "80 and above" )))


ggplot(all_pre_metadata, aes(group, fill = group )) +
  geom_bar()+
  ggtitle("Distribution of baseline data by age group of participants")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, face = "bold", size = 13),
        axis.text.y = element_text(face = "bold", size = 13))+
  ylab("Number of participants")+
  xlab("Age group of participants")
  




# Load the gene annotation file
gene_annotation <- readRDS("data_new/ensembl_gene_annotation.RDS")


# The second model looks at the change in Se by age groups
group_model <- readRDS("data_new/models/grouped_model.RDS") %>%
  filter(Pr...z..<= 0.05 ) %>%
  separate("target", c("transcript_ID", "intron_ID", "chr"), sep = "_") %>%
  inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version")) %>%
  dplyr::select(coef, Estimate, transcript_ID, intron_ID,  transcript_biotype, ensembl_gene_id, external_gene_name, transcript_length)%>%
  # merge the intron_id back to the transcript ID so we track specific introns 
  mutate(intron_ID = paste0(transcript_ID, "_", intron_ID),
         coef = recode(coef, "group21 to 29" = "21 to 29", 
                       "group30 to 39" = "30 to 39",
                       "group40 to 49" = "40 to 49",
                       "group50 to 59" = "50 to 59",
                       "group60 to 69" =  "60 to 69",
                       "group70 to 79" = "70 to 79",
                       "group80 and above" = "80 and above" )) 

# Visualize the coefficients
all <- ggplot(group_model , aes(coef, fill = coef )) +
  geom_bar()+
  ggtitle("Distribution of differentially spliced introns by age group")+
  xlab(" Age group of participants")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 10))
  
  
  

  
  # Twenties
  twenties <- group_model %>%
    filter(coef == "21 to 29")
 a <-  twenties %>%
    ggplot()+
    geom_histogram(aes(x = Estimate)) +
     ylab("Number of ds introns per Estimate")+
    # xlab("Change in SE by aging and by RT")+
     ggtitle("Distribution of ds introns in participants aged 21 to 29")+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))+
    theme_cowplot()
  
  
  
# Gene ontology analyses
  
# load the batch-corrected gene counts data
gene_counts <- readRDS("data_new/gene_counts/batch_corrected_genecounts.RDS") %>%
 # remove the version number to match gene_id
 mutate(gene_id = gsub("\\..*", "",  gene_id))



# gene ontology (biological process)  analysis of the genes containing introns with improved SE with aging
ego_twenties <- enrichGO(gene = twenties$ensembl_gene_id,
                      keyType = "ENSEMBL",
                      OrgDb = org.Hs.eg.db, 
                      universe = gene_counts$gene_id,
                      ont = "BP", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05, 
                      readable = T)



twenties_plot<-  dotplot(ego_twenties,
             
             font.size = 8, title = "Enriched biological processes for those aged 21 to 29") +
  theme(axis.text = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12))



# Thirties
thirties <- group_model %>%
  filter(coef == "30 to 39")
b <- thirties %>%
  ggplot()+
  geom_histogram(aes(x = Estimate)) +
  ylab("Number of ds introns per Estimate")+
  # xlab("Change in SE by aging and by RT")+
  ggtitle("Distribution of ds introns in participants aged 30 to 39")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))+
  theme_cowplot()






# gene ontology (biological process)  analysis of the genes containing introns with improved SE with aging
ego_thirties <- enrichGO(gene = thirties$ensembl_gene_id,
                         keyType = "ENSEMBL",
                         OrgDb = org.Hs.eg.db, 
                         universe = gene_counts$gene_id,
                         ont = "BP", 
                         pAdjustMethod = "BH", 
                         qvalueCutoff = 0.05, 
                         readable = T)



thirties_plot<-  dotplot(ego_thirties,
                         
                         font.size = 8, title = "Enriched biological processes for those aged 30 to 39") +
  theme(axis.text = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12))





# Fourties
fourties <- group_model %>%
  filter(coef == "40 to 49")
c <- fourties %>%
  ggplot()+
  geom_histogram(aes(x = Estimate)) +
  ylab("Number of ds introns per Estimate")+
  # xlab("Change in SE by aging and by RT")+
  ggtitle("Distribution of ds introns in participants aged 40 to 49")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))+
  theme_cowplot()






# gene ontology (biological process)  analysis of the genes containing introns with improved SE with aging
ego_fourties <- enrichGO(gene = fourties$ensembl_gene_id,
                         keyType = "ENSEMBL",
                         OrgDb = org.Hs.eg.db, 
                         universe = gene_counts$gene_id,
                         ont = "BP", 
                         pAdjustMethod = "BH", 
                         qvalueCutoff = 0.05, 
                         readable = T)



fourties_plot<-  dotplot(ego_fourties,
                         
                         font.size = 8, title = "Enriched biological processes for those aged 40 to 49") +
  theme(axis.text = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12))





# Fifties
fifties <- group_model %>%
  filter(coef == "50 to 59")
d <- fifties %>%
  ggplot()+
  geom_histogram(aes(x = Estimate)) +
  ylab("Number of ds introns per Estimate")+
  # xlab("Change in SE by aging and by RT")+
  ggtitle("Distribution of ds introns in participants aged 50 to 59")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))+
  theme_cowplot()






# gene ontology (biological process)  analysis of the genes containing introns with improved SE with aging
ego_fifties <- enrichGO(gene = fifties$ensembl_gene_id,
                         keyType = "ENSEMBL",
                         OrgDb = org.Hs.eg.db, 
                         universe = gene_counts$gene_id,
                         ont = "BP", 
                         pAdjustMethod = "BH", 
                         qvalueCutoff = 0.05, 
                         readable = T)



fifties_plot<-  dotplot(ego_fifties,
                         
                         font.size = 8, title = "Enriched biological processes for those aged 50 to 59") +
  theme(axis.text = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12))






# Sixties
sixties <- group_model %>%
  filter(coef == "60 to 69")
e <- sixties %>%
  ggplot()+
  geom_histogram(aes(x = Estimate)) +
  ylab("Number of ds introns per Estimate")+
  # xlab("Change in SE by aging and by RT")+
  ggtitle("Distribution of ds introns in participants aged 60 to 69")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))+
  theme_cowplot()






# gene ontology (biological process)  analysis of the genes containing introns with improved SE with aging
ego_sixties <- enrichGO(gene = sixties$ensembl_gene_id,
                         keyType = "ENSEMBL",
                         OrgDb = org.Hs.eg.db, 
                         universe = gene_counts$gene_id,
                         ont = "BP", 
                         pAdjustMethod = "BH", 
                         qvalueCutoff = 0.05, 
                         readable = T)



sixties_plot<-  dotplot(ego_sixties,
                         
                         font.size = 8, title = "Enriched biological processes for those aged 60 to 69") +
  theme(axis.text = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12))



# seventies 
seventies <- group_model %>%
  filter(coef == "70 to 79")
f <- seventies %>%
  ggplot()+
  geom_histogram(aes(x = Estimate)) +
  ylab("Number of ds introns per Estimate")+
  # xlab("Change in SE by aging and by RT")+
  ggtitle("Distribution of ds introns in participants aged 70 to 79")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))+
  theme_cowplot()






# gene ontology (biological process)  analysis of the genes containing introns with improved SE with aging
ego_seventies <- enrichGO(gene = seventies$ensembl_gene_id,
                         keyType = "ENSEMBL",
                         OrgDb = org.Hs.eg.db, 
                         universe = gene_counts$gene_id,
                         ont = "BP", 
                         pAdjustMethod = "BH", 
                         qvalueCutoff = 0.05, 
                         readable = T)



seventies_plot<-  dotplot(ego_seventies,
                         
                         font.size = 8, title = "Enriched biological processes for those aged 70 to 79") +
  theme(axis.text = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12))



# Eighties
eighties <- group_model %>%
  filter(coef == "80 and above")
g <- eighties %>%
  ggplot()+
  geom_histogram(aes(x = Estimate)) +
  ylab("Number of ds introns per Estimate")+
  # xlab("Change in SE by aging and by RT")+
  ggtitle("Distribution of ds introns in participants aged 80 and above")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))+
  theme_cowplot()






# gene ontology (biological process)  analysis of the genes containing introns with improved SE with aging
ego_eighties <- enrichGO(gene = eighties$ensembl_gene_id,
                         keyType = "ENSEMBL",
                         OrgDb = org.Hs.eg.db, 
                         universe = gene_counts$gene_id,
                         ont = "BP", 
                         pAdjustMethod = "BH", 
                         qvalueCutoff = 0.05, 
                         readable = T)



eighties_plot<-  dotplot(ego_eighties,
                         
                         font.size = 8, title = "Enriched biological processes for those aged 80 and above") +
  theme(axis.text = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12))


ggarrange(all, a + rremove("xy.title"), b + rremove("xy.title"), 
          c + rremove("xy.title") , d + rremove("xy.title") , e + rremove("xy.title") , f + rremove("y.title"),
          g + rremove("y.title"), ncol = 2, nrow = 4,
         # align = "v",
          labels = c("A", "B", "C",
                     "D", "E", "F", "G", "H"))
ggsave("Images_tables/Figure2_grouped_baseline.png", bg = "white" ,  scale = 2, dpi = 400)



# plot the gene ontology results
ggarrange(twenties_plot, thirties_plot, fourties_plot,
          fifties_plot, sixties_plot, seventies_plot,
          eighties_plot, align = "hv", ncol = 2, nrow = 4
          )

ggsave("Images_tables/Figure_grouped_baselineGO.png", bg = "white" ,  scale = 2, dpi = 400)
