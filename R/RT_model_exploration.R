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


RT_effects <- model %>%
  filter(coef == "timePostExc" | coef == "scaled_age:timePostExc") %>%
  #filter(fcthreshold == "s") %>%
  separate("target", c("transcript_ID", "intron_ID", "chr"), sep = "_") %>%
  inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version")) %>%
  dplyr::select(Estimate, coef,  transcript_ID, transcript_biotype, ensembl_gene_id, external_gene_name, transcript_length)


#hist(sev_above$Estimate)
# Load the gene annotation extracted from Ensemble
gene_annotation <- readRDS("data_new/ensembl_gene_annotation.RDS")



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
hist(filt_exc$Estimate)


ego_df_postEXC <- enrichGO(gene = filt_exc$ensembl_gene_id,
                      keyType = "ENSEMBL",
                      OrgDb = org.Hs.eg.db, 
                      ont = "MF", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05, 
                      readable = T)

## Output results from GO analysis to a table
cluster_summary_postExc <- data.frame(ego_df_postEXC)

dotplot(ego_df_postEXC,
             
             font.size = 8, title = "Enriched molecular functions in genes containing differentially at postexc") +
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))



filt_int <- filt %>%
  filter(coef == "scaled_age:timePostExc")
hist(filt_int$Estimate)

cor.test(filt_int$Estimate, filt_int$transcript_length)
plot(filt_int$Estimate, filt_int$transcript_length)


filt_int_pos <- filt_int %>%
  filter(Estimate > 0)

filt_int_neg <- filt_int %>%
  filter(Estimate < 0)

ego_df_postEXC <- enrichGO(gene = filt_int_neg$ensembl_gene_id,
                           keyType = "ENSEMBL",
                           OrgDb = org.Hs.eg.db, 
                           ont = "MF", 
                           pAdjustMethod = "BH", 
                           qvalueCutoff = 0.05, 
                           readable = T)

## Output results from GO analysis to a table
cluster_summary_postExc <- data.frame(ego_df_postEXC)

a <- dotplot(ego_df_postEXC,
        
        font.size = 8, title = "Enriched molecular functions in genes whose SE reduced postexc") +
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))


ego_df_postEXC2 <- enrichGO(gene = filt_int_pos$ensembl_gene_id,
                           keyType = "ENSEMBL",
                           OrgDb = org.Hs.eg.db, 
                           ont = "MF", 
                           pAdjustMethod = "BH", 
                           qvalueCutoff = 0.05, 
                           readable = T)

## Output results from GO analysis to a table
cluster_summary_postExc2 <- data.frame(ego_df_postEXC2)

b <- dotplot(ego_df_postEXC2,
             
             font.size = 8, title = "Enriched molecular functions in genes whose SE improved postexc") +
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))












# extract those ds by interactionwith age




# Extract the annotation 
sev_annotation <- inner_join(sev_above, En_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version"))


ego_df <- enrichGO(gene = sev_annotation$ensembl_gene_id,
                   keyType = "ENSEMBL",
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = T)

## Output results from GO analysis to a table
# cluster_summary <- data.frame(ego_df)

dotplot(ego_df,
        
        font.size = 8, title = "Enriched biological processes of above 70s at post exercise") +
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))

ggsave("Figures/over_70_not_interaction_BP.png")


# eXtract introns that are differentially spliced by those above 70 at postexercise

ds_70_Rt <- model_3 %>%
  filter(coef == "group>20 & <=30:timePostExc") %>%
  pull(target)

# Load the full splicing data

all_splice_df<- readRDS("data_new/processed_data/all_splice_data.RDS")

ds_70_df <- all_splice_df[all_splice_df$transcript_ID %in% ds_70_Rt,]

# LOad full metadata

all_full_metadata <- readRDS("data_new/processed_data/all_full_metadata.RDS")


# Visualization
ds_70_df_long<- ds_70_df%>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID) )%>%
  inner_join(all_full_metadata, by = "seq_sample_id")
 

hist(ds_70_df_long$SE)


# Plot the relationship between age and splicing efficiency
ds_70_df_long %>%
  group_by(group, time )%>%
  summarise(avg = mean(SE)) %>%
  ggplot(aes(time,avg))+
  geom_point(mapping = aes(colour = group, size = 10))+ 
  geom_smooth()+
  ggtitle("Relationship between age group and splicing efficiency") +
  ylab(" Average splicing efficiency")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_cowplot() 

# Visualize using heatmap

ds_70_df_long %>%
  ggplot(mapping = aes(y = transcript_ID,
                          x = time,
                          fill = SE,
                       colour = group))+
  geom_tile() +
  facet_grid(~group)
