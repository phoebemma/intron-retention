library(dplyr)
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(org.Hs.eg.db)
library(clusterProfiler)

# This model analyses the interaction between  the age_groups and RT, and sex
model <- readRDS("data_new/models/full_data_RT_and_age_model.RDS") %>%
  # filter to only those with adjusted p values at or below 0.05
  filter(adj.p <= 0.05)

length(unique(model$target))

# This model analyses only exercise and group
model_RT <- readRDS("data_new/models/full_data_RT_model.RDS")%>%
  filter(adj.p <= 0.05)
length(unique(model_RT$target))
hist(model_RT$Estimate)
hist(model$Estimate)


model %>%
  ggplot(aes(x = coef)) +
  geom_bar()+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

# This looked at interaction between sex and RT, as well as group and RT

model_3 <- readRDS("data_new/models/full_data_RT_interction_group_sex_model.RDS")%>%
  filter(adj.p <= 0.05)
length(unique(model_3$target))
unique(model_3$coef)

model_3 %>%
  ggplot(aes(x = coef)) +
  geom_bar()+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))


# filter the ds of 70 and above's interaction with time
sev_above <- model_3 %>%
  filter(coef == "group>70:timePostExc")%>%
  # seperate it into the original format 
  separate("target",  c("transcript_ID", "intron_ID", "chr"), sep = "_") %>%
  select(transcript_ID,intron_ID )

#hist(sev_above$Estimate)
# Load the gene annotation extracted from Ensemble
En_annotation <- readRDS("data_new/ensembl_gene_annotation.RDS")


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
