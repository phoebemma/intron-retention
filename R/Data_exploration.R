# This script explores the splicing data and metadata 

library(dplyr)
library(ggplot2)
library(biomaRt)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

# Load the combined pre exercise metadata and splicing datsets

# metadata
all_pre_metadata <- readRDS("data_new/Pre_Exercise/all_prexercise_metadata.RDS")

# splicing data
all_pre_splice <- readRDS("data_new/Pre_Exercise/all_pre_Exc_splicing_data.RDS")


# Are there sex-specific intron retentions

sex_spec <- all_pre_splice %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID))%>%
  inner_join(all_pre_metadata, by = "seq_sample_id") %>%
  summarise(median = median(SE), 
            min = min(SE), 
            max = max(SE),
            s = sd(SE), 
            .by = c(sex, transcript_ID)) %>%
  filter( s > 0.099) %>%
  separate("transcript_ID", c("transcript_ID", "intron_ID", "chr"), sep = "_") 
hist(sex_spec$s)

# Are there specific introns that charactaristically have low splicing efficiency

low_SE <- all_pre_splice %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID)) %>%
#  inner_join(all_pre_metadata, by = "seq_sample_id") %>%
  summarise(.by = transcript_ID, 
            me = median(SE),
            min = min(SE), 
            max = max(SE), 
            q20 = quantile(SE, 0.2), 
            range = max(SE) - min(SE)) %>%
  filter(max <= 0.6) %>%
  separate("transcript_ID", c("transcript_ID", "intron_ID", "chr"), sep = "_")
length(unique(low_SE$transcript_ID))
hist(low_SE$q20)



High_SE <- all_pre_splice %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID)) %>%
  summarise(.by = transcript_ID, 
            me = median(SE),
            min = min(SE), 
            max = max(SE)) %>%
  filter(min == 1)%>%
  separate("transcript_ID", c("transcript_ID", "intron_ID", "chr"), sep = "_")

#Get the ensemble annotation of genes
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "hsapiens_gene_ensembl",
                   verbose = TRUE)

attributes <- listAttributes(ensembl)

# #extract GENE biotypes and  names
annotation<- getBM(attributes = c("transcript_biotype", "external_transcript_name",
                                  "ensembl_gene_id", "ensembl_transcript_id_version", "external_gene_name"),  mart = ensembl )
annotation_low_SE <- inner_join(low_SE, annotation, by= c("transcript_ID" = "ensembl_transcript_id_version"))

annotation_high_SE <- inner_join(High_SE, annotation, by= c("transcript_ID" = "ensembl_transcript_id_version"))

annotation_high_SE %>%
  ggplot(aes(transcript_biotype, fill = transcript_biotype))+
  geom_bar()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = -0.5)


saveRDS(annotation_low_SE, "data_new/Pre_Exercise/annotated_low_SE_introns.RDS")

saveRDS(annotation_high_SE, "data_new/Pre_Exercise/annotated_perfect_SE_introns.RDS")

saveRDS(annotation, "data_new/ensembl_gene_annotation.RDS")

ego_df <- enrichGO(gene = annotation_low_SE$ensembl_gene_id,
                   keyType = "ENSEMBL",
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = T)

## Output results from GO analysis to a table
cluster_summary <- data.frame(ego_df)

dotplot(ego_df,
        
        font.size = 8, title = "Enriched biological processes in completely spliced out  introns") +
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))

