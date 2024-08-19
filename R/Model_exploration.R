#This is where the models built are explored

library(dplyr)
library(ggplot2)
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(plotrix)

#Load the prexercise model

pre_model <- readRDS("data/model/filtered_full_model.RDS") %>%
  dplyr::select(coef, Estimate, target, Pr...z..)

hist(pre_model$Estimate)
hist(pre_model$Pr...z..)

#The target contains the transcript_id, intron_id and chromosome in thats order.

# see "read_Splice_Q" in Trainome_functions.R
colnames(pre_model)

pre_df <- pre_model %>%
  separate("target", c("transcript_ID", "intron_ID", "chr"), sep = "_") %>%
  subset(coef == "timePostExc")

#unique(pre_df$coef)
ggplot(pre_df, aes(chr)) +
  geom_bar()+
  labs(x = "chromosome")+
  ggtitle("Distribution of chromosomes containing differentially retained introns  by exercise")+
theme(plot.title = element_text(hjust = 0.5))

length(unique(pre_df$transcript_ID))


#Use the ensemble database to get the annotation of the transcripts
listMarts()
ensembl=useMart("ENSEMBL_MART_ENSEMBL")
datasets <- listDatasets(ensembl)
ensembl=useDataset("hsapiens_gene_ensembl",mart=ensembl)
attributes <- listAttributes(ensembl)

# #extract GENE biotypes and  names
annotation<- getBM(attributes = c("transcript_biotype", "external_transcript_name",
                                  "ensembl_gene_id", "ensembl_transcript_id_version"), values =pre_df$transcript_ID, mart = ensembl )
annotation <- inner_join(pre_df, annotation, by= c("transcript_ID" = "ensembl_transcript_id_version"))

#Explore the biotypes of the transcripts
x <- annotation %>%
  group_by(transcript_biotype)%>%
  count()%>% #count the occurence
  ungroup()%>%
  mutate(perc = n/sum(n))%>% #extract the percentage
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))
  

ggplot(x, aes(x = "", y = perc, fill = transcript_biotype)) +
  geom_col()+
  geom_text(aes(label = labels),
            position = position_stack(vjust = 0.8)) +
  coord_polar(theta = "y")
  
pie3D(x$perc, labels = x$transcript_biotype,
      #explode = 0.2,
      labelcex = 0.65,
      radius = 1, theta = pi/1)


#Geneontology
ego_df <- enrichGO(gene = annotation$ensembl_gene_id,
                   keyType = "ENSEMBL",
                   OrgDb = org.Hs.eg.db,
                   ont = "bp",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   readable = T)

dotplot(ego_df, showCategory = 15,
        
        font.size = 5, title = "15 top biological processes of genes containing DE introns by exercise alone") +
  theme(axis.text = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 10))



#Enrichment analyses

#get the entrezid of the unique genes
entrez_ids <- bitr(annotation$ensembl_gene_id,"ENSEMBL", "ENTREZID", org.Hs.eg.db)

kegg_df <- enrichKEGG(gene = entrez_ids$ENTREZID,,
                      organism = "hsa",
                      keyType = "kegg",
                      # OrgDb = org.Hs.eg.db, 
                      #ont = "MF", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05)
barplot(kegg_df, showCategory = 15, title = "10 most enriched pathways of genes containing differeantially retained introns by exercise alone")+
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 13))
