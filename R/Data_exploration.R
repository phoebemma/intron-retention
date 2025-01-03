# This script explores the splicing data and metadata 

library(dplyr)
library(ggplot2)
library(biomaRt)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)


# metadata
#Load the metadata of the the datasets

#Copd 
copd_metadata <- readRDS("data_new/Pre_Exercise/copd_preExc_metadata.RDS") %>%
  select(study, participant, sex, time, seq_sample_id, age)
unique(copd_metadata$time)
colnames(copd_metadata)



#Volume

volume_metadata <- readRDS("data_new/Pre_Exercise/vol_preExc_metadata.RDS")%>%
  select(study, participant, sex, time, seq_sample_id, age)

colnames(volume_metadata)
unique(volume_metadata$time)


# Contratratrain
Contratrain_metadata <- readRDS("data_new/Pre_Exercise/ct_PreExc_metadata.RDS")%>%
  select(study, participant, sex, time, seq_sample_id, age)
colnames(Contratrain_metadata)
unique(Contratrain_metadata$sex)

# SRP102542
SRP102542_metadata <- readRDS("data_new/Pre_Exercise/SRP102542_preExc_metadata.RDS")%>%
  select(study, participant, sex, time, seq_sample_id, age)
colnames(SRP102542_metadata)
unique(SRP102542_metadata$sex)

range(SRP102542_metadata$age)

# SRP280348_metadata <- readRDS("data/preexercise_data/SRP280348_preExc_metadata.RDS")%>%
#   select(study, participant, sex, time, seq_sample_id, age, age_group)
# colnames(SRP280348_metadata)
# unique(SRP280348_metadata$sex)
# #


Alpha_Omega_metadata <- readRDS("data_new/Pre_Exercise/Alpha_Omega_PreExc_metadata.RDS")
# Merge all in one
all_pre_metadata <- rbind(copd_metadata, volume_metadata)%>%
  rbind(Contratrain_metadata) %>%
  rbind(Alpha_Omega_metadata) %>%
  rbind(SRP102542_metadata) %>%
  #Copd and volume and SRP102542 have age data in decimal
  mutate(across(c("age"), round, 0)) %>%
  # Create youps where group 1 = those below 30
  #group 2 is those above 30 but below 51
  # group 3 those above 50 but below 71
  # group 4 is those above 70
  
  mutate(group = case_when(age <=20 ~ "<=20" ,
                           age > 20 & age <= 30 ~ ">20 & <=30",
                           age > 30 & age <= 40 ~ ">30 & <=40", 
                           age > 40 & age <= 50 ~ ">40 & <=50",
                           age > 50 & age <= 60 ~ ">50 & <=60",
                           age > 60 & age <= 70 ~ ">60 & <=70",
                           age > 70 ~ ">70")) %>%
  mutate(sex = factor(sex, levels = c("female", "male")),
         group = factor(group, levels = c("<=20" ,">20 & <=30", ">30 & <=40", ">40 & <=50",
                                          ">50 & <=60",">60 & <=70", ">70"))) 
unique(all_pre_metadata$sex)
length(unique(all_pre_metadata$participant))

# ggplot(all_pre_metadata, aes(age, fill = group )) +
#   geom_bar()+
#   ggtitle("Distribution of baseline data")+
#   theme(plot.title = element_text(hjust = 0.5))


# saveRDS(all_pre_metadata, "data_new/Pre_Exercise/all_prexercise_metadata.RDS")

ggplot(all_pre_metadata, aes(group, fill = group)) +
  geom_bar()+
  ggtitle("Distribution of baseline data")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("Age range of participants")+
  stat_count(geom = "Text", aes(label = ..count..), vjust = 1.5)


# ggplot(all_pre_metadata, aes(age_group, fill = age_group)) +
#   geom_bar()+
#   ggtitle("Distribution of baseline data")+
#   theme(plot.title = element_text(hjust = 0.5))+
#   stat_count(geom = "Text", aes(label = ..count..), vjust = 1.5)+
#   xlab("Age group of participants")



ggplot(all_pre_metadata, aes(study, fill = group)) +
  geom_bar()+
  ggtitle("Distribution of baseline data")+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = 1.5)


ggplot(all_pre_metadata, aes(sex, fill = group)) +
  geom_bar()+
  ggtitle("Distribution of baseline data")+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = 1.5)

length(all_pre_metadata$participant)




copd_data <- readRDS("data_new/Pre_Exercise/copd_preExc_splicing_data.RDS")


volume_data <- readRDS("data_new/Pre_Exercise/vol_preExc_splicing_data.RDS")

contratrain_data <- readRDS("data_new/Pre_Exercise/ct_PreExc_splicing_data.RDS")

SRP102542_data <- readRDS("data_new/Pre_Exercise/SRP102542_preExc_splicing_data.RDS")

Alpha_Omega_data <- readRDS("data_new/Pre_Exercise/Alpha_Omega_PreExc_splicing_data.RDS")


all_pre_splice <- copd_data%>%
  inner_join(volume_data, by = "transcript_ID") %>%
  inner_join(contratrain_data, by = "transcript_ID")%>%
  inner_join(SRP102542_data, by = "transcript_ID") %>%
  inner_join(Alpha_Omega_data, by = "transcript_ID") %>%
  drop_na()




# Visualization
long_df <- all_pre_splice%>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID) )%>%
  inner_join(all_pre_metadata, by = "seq_sample_id")



# Plot the relationship between age and splicing efficiency
long_df %>%
  group_by(group)%>%
  summarise(avg = mean(SE)) %>%
  ggplot(aes(avg, group))+
  geom_point(mapping = aes(colour = group, size = 10))+ 
  geom_smooth()+
  ggtitle("Relationship between age group and splicing efficiency") +
  xlab(" Average splicing efficiency")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_cowplot()


long_df %>%
  group_by(age, sex)%>%
  summarise(avg = mean(SE)) %>%
  ggplot(aes(age, avg))+
  geom_point(mapping = aes(colour = sex ,  size = 5))+ 
  geom_smooth()+
  ggtitle("Relationship between age, gender and splicing efficiency") +
  xlab(" Average splicing efficiency")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_cowplot()




ggarrange(x, z)


saveRDS(all_pre_splice, "data_new/Pre_Exercise/all_pre_Exc_splicing_data.RDS")

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
            .by = c(sex,  transcript_ID)) %>%
  filter( max < 0.5) %>%
  separate("transcript_ID", c("transcript_ID", "intron_ID", "chr"), sep = "_") 
hist(sex_spec$s)
unique(sex_spec$transcript_ID)
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

# Investigate introns with low SE across groups. 
# Appears to be same introns that are low across all samples
low_SE_grouped <- all_pre_splice %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID)) %>%
  inner_join(all_pre_metadata, by = "seq_sample_id") %>%
  summarise(.by = c(transcript_ID, group),
            me = median(SE),
            min = min(SE), 
            max = max(SE), 
            q20 = quantile(SE, 0.2), 
            range = max(SE) - min(SE)) %>%
  filter(max <= 0.6) %>%
  separate("transcript_ID", c("transcript_ID", "intron_ID", "chr"), sep = "_")

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

length(unique(High_SE$transcript_ID))

# Get the ensemble annotation of genes
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

