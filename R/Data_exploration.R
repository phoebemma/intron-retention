# This script explores the splicing data and metadata 

library(dplyr)
library(ggplot2)
library(biomaRt)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)


# metadata
#Load the metadata of the the datasets

#Copd 
copd_metadata <- readRDS("data_new/Pre_Exercise/copd_preExc_metadata.RDS") %>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)
unique(copd_metadata$time)
colnames(copd_metadata)
length(unique(copd_metadata$participant))
length(unique(copd_metadata$seq_sample_id))


#Volume

volume_metadata <- readRDS("data_new/Pre_Exercise/vol_preExc_metadata.RDS")%>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)

colnames(volume_metadata)
unique(volume_metadata$time)
length(unique(volume_metadata$participant))
length(unique(volume_metadata$seq_sample_id))
# Contratratrain
Contratrain_metadata <- readRDS("data_new/Pre_Exercise/ct_PreExc_metadata.RDS")%>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)
colnames(Contratrain_metadata)
unique(Contratrain_metadata$sex)
length(unique(Contratrain_metadata$participant))
length(unique(Contratrain_metadata$seq_sample_id))
range(Contratrain_metadata$age)

# SRP102542
SRP102542_metadata <- readRDS("data_new/Pre_Exercise/SRP102542_preExc_metadata.RDS")%>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)
colnames(SRP102542_metadata)
unique(SRP102542_metadata$sex)
length(unique(SRP102542_metadata$participant))
range(SRP102542_metadata$age)

# SRP280348_metadata <- readRDS("data/preexercise_data/SRP280348_preExc_metadata.RDS")%>%
#   select(study, participant, sex, time, seq_sample_id, age, age_group)
# colnames(SRP280348_metadata)
# unique(SRP280348_metadata$sex)
# #


Alpha_Omega_metadata <- readRDS("data_new/Pre_Exercise/Alpha_Omega_PreExc_metadata.RDS")
length(unique(Alpha_Omega_metadata$seq_sample_id))

length(unique(Alpha_Omega_metadata$participant))
range(Alpha_Omega_metadata$age)
# Merge all in one
all_pre_metadata <- rbind(copd_metadata, volume_metadata)%>%
  rbind(Contratrain_metadata) %>%
  rbind(Alpha_Omega_metadata) %>%
  rbind(SRP102542_metadata) %>%
  #Copd and volume and SRP102542 have age data in decimal
  mutate(across(c("age"), round, 0)) %>%
  # Create age groups where group 1 = those 20 and below
  #group 2 is those above 20 but below 31
  # group 3 those above 30 but below 41
  # etc
  
  mutate(group = case_when(age <=20 ~ "<=20" ,
                           age > 20 & age <= 30 ~ ">20 & <=30",
                           age > 30 & age <= 40 ~ ">30 & <=40", 
                           age > 40 & age <= 50 ~ ">40 & <=50",
                           age > 50 & age <= 60 ~ ">50 & <=60",
                           age > 60 & age <= 70 ~ ">60 & <=70",
                           age > 70 ~ ">70"),
         age_class = case_when(age <= 40 ~ "Young",
                               age > 40 & age <= 60 ~ "Middle aged",
                               age > 60 ~ "Old")) %>%
  mutate(sex = factor(sex, levels = c("female", "male")),
         group = factor(group, levels = c("<=20" ,">20 & <=30", ">30 & <=40", ">40 & <=50",
                                          ">50 & <=60",">60 & <=70", ">70")),
         age_class = factor(age_class, levels = c("Young","Middle aged", "Old" ))) 
unique(all_pre_metadata$sex)
length(unique(all_pre_metadata$participant))
unique(all_pre_metadata$age_class)

all

 saveRDS(all_pre_metadata, "data_new/Pre_Exercise/all_prexercise_metadata.RDS")

a<- ggplot(all_pre_metadata, aes(group, fill = group)) +
  geom_bar()+
#  ggtitle("Distribution of baseline data")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("Age range of participants")+
  stat_count(geom = "Text", aes(label = ..count..), vjust = 1.5)

# ggsave("Figures/Baseline_data.png")

d <- ggplot(all_pre_metadata, aes(age_class, fill = age_class)) +
  geom_bar()+
  #  ggtitle("Distribution of baseline data")+
  theme(plot.title = element_text(hjust = 0.5))+
#  xlab("Age c of participants")+
  stat_count(geom = "Text", aes(label = ..count..), vjust = 1.5)



b <- ggplot(all_pre_metadata, aes(study, fill = group)) +
  geom_bar()+
#  ggtitle("Distribution of baseline data")+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = 1.5)

# ggsave("Figures/Baseline_data_by_study.png")

c <- ggplot(all_pre_metadata, aes(sex, fill = group)) +
  geom_bar()+
#  ggtitle("Distribution of baseline data")+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = 1.5)

# ggsave("Figures/Baseline_data_by_gender.png")
length(all_pre_metadata$participant)

ggarrange(a ,b  , c , d, 
          labels = c("A", "B", "C", "D"), 
          common.legend = T,
          align = "hv",
          hjust = -1,
          legend = "right")

ggsave("Figures/Baseline_data.png", scale = 2, dpi = 400, bg = "white")


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
  ggplot(aes(group,avg))+
  geom_point(mapping = aes(colour = group, size = 10))+ 
  geom_smooth()+
  ggtitle("Relationship between age group and splicing efficiency") +
  ylab(" Average splicing efficiency")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_cowplot()
ggsave("Figures/age_group.png", bg = "white")

long_df %>%
  group_by(age, sex)%>%
  summarise(avg = mean(SE)) %>%
  ggplot(aes(age, avg))+
  geom_point(mapping = aes(colour = sex ,  size = 5))+ 
  geom_smooth()+
  ggtitle("Relationship between age, gender and splicing efficiency") +
  ylab(" Average splicing efficiency")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_cowplot()

ggsave("Figures/age_group_by_gender.png", bg = "white")


# ggarrange(x, z)


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
            mean = mean(SE),
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
  filter(min == 1) %>%
  separate("transcript_ID", c("transcript_ID", "intron_ID", "chr"), sep = "_")

length(unique(High_SE$transcript_ID))

# Get the ensemble annotation of genes
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "hsapiens_gene_ensembl",
                   verbose = TRUE)

attributes <- listAttributes(ensembl)

# #extract GENE biotypes and  names
annotation<- getBM(attributes = c("transcript_biotype", "external_transcript_name",
                                  "ensembl_gene_id", "ensembl_transcript_id_version", "external_gene_name", "transcript_length"),  mart = ensembl )
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



# Gene annotation
ego_df <- enrichGO(gene = annotation_low_SE$ensembl_gene_id,
                   keyType = "ENSEMBL",
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = T)

## Output results from GO analysis to a table
# cluster_summary <- data.frame(ego_df)

dotplot(ego_df,
        
        font.size = 8, title = "Enriched biological processes in poorly spliced out  introns") +
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))
ggsave("Figures/GO_poorly_spliced_introns.png")


ego_df_high <- enrichGO(gene = annotation_high_SE$ensembl_gene_id,
                   keyType = "ENSEMBL",
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = T)

## Output results from GO analysis to a table

dotplot(ego_df_high,
        
        font.size = 8, title = "Enriched biological processes in perfectly spliced out  introns") +
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))

ggsave("Figures/GO_perfectly_spliced_introns.png")



# Explore what happens to those introns at post exercise

#Load the full data of all types


#COPD metadata
copd_metadata <- readRDS("data_new/processed_data/copd_metadata.RDS") %>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)

# Volume_data
Vol_metadata <- readRDS("data_new/processed_data/volume_metadata.RDS")%>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age) 

unique(Vol_metadata$time)

# Contratrain_data
ct_metadata <- readRDS("data_new/processed_data/contratrain_metadata.RDS") %>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)


unique(ct_metadata$time)

# Publicly available data
SRP102542_metadata <- readRDS("data_new/processed_data/SRP102542_metadata.RDS")%>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)
unique(SRP102542_metadata$time)



# Alpha and Omega data

A_Omega_metadata <- readRDS("data_new/processed_data/Alpha_Omega_metadata.RDS")%>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)
unique(A_Omega_metadata$time)


# Merge them all in one

all_full_metadata <- rbind(copd_metadata, Vol_metadata)%>%
  rbind(ct_metadata) %>%
  rbind(A_Omega_metadata) %>%
  rbind(SRP102542_metadata) %>%
  
  mutate(across(c("age"), round, 0)) %>%
  
  mutate(group = case_when(age <=20 ~ "<=20" ,
                           age > 20 & age <= 30 ~ ">20 & <=30",
                           age > 30 & age <= 40 ~ ">30 & <=40", 
                           age > 40 & age <= 50 ~ ">40 & <=50",
                           age > 50 & age <= 60 ~ ">50 & <=60",
                           age > 60 & age <= 70 ~ ">60 & <=70",
                           age > 70 ~ ">70")) %>%
  mutate(sex = factor(sex, levels = c("female", "male")),
         time = factor(time, levels = c("PreExc", "PostExc")),
         group = factor(group, levels = c("<=20" ,">20 & <=30", ">30 & <=40", ">40 & <=50",
                                          ">50 & <=60",">60 & <=70", ">70"))) 

unique(all_full_metadata$sex)
unique(all_full_metadata$study)
unique(all_full_metadata$time)
unique(all_full_metadata$group)

saveRDS(all_full_metadata, "data_new/processed_data/all_full_metadata.RDS")


# Extract only post exercise data
post_metadata <- all_full_metadata %>%
  filter(time == "PostExc")




copd_splice_df <- readRDS("data_new/processed_data/copd_splicing_data.RDS")

vol_splice_df <- readRDS("data_new/processed_data/volume_splicing_data.RDS")
ct_splice_df <- readRDS("data_new/processed_data/contratrain_splicing_data.RDS")
SRP102542_splice_df <- readRDS("data_new/processed_data/SRP102542_splicing_data.RDS")
AOD_splice_df <- readRDS("data_new/processed_data/Alpha_Omega_splicing_data.RDS")

all_splice_df <-copd_splice_df  %>%
  inner_join(vol_splice_df, by = "transcript_ID") %>%
  inner_join(ct_splice_df, by = "transcript_ID") %>%
  inner_join(AOD_splice_df, by = "transcript_ID") %>%
  inner_join(SRP102542_splice_df, by = "transcript_ID") 


# Get the post exercise splicing data
# Getting it before removing the NAs in the full data ensures more introns are captured

post_intersect <- intersect(colnames(all_splice_df),post_metadata$seq_sample_id)

post_splice_df <- all_splice_df %>%
  subset(select = c("transcript_ID", post_intersect))%>%
  drop_na()


# select only the splicing samples captured in the metadata
all_intersect <- intersect(colnames(all_splice_df), all_full_metadata$seq_sample_id)

all_spliceq_df <- all_splice_df %>%
  subset(select = c("transcript_ID", all_intersect))%>%
  drop_na()

saveRDS(all_spliceq_df, "data_new/processed_data/all_splice_data.RDS")





#a <- readr::read_tsv("data_new/Alpha_Omega_SpliceQ_outputs/s100_EKRN240058206.tsv")
# Query the low and perfect spliced introns to see how they fared postexercise

High_at_post <- all_spliceq_df[all_spliceq_df$transcript_ID %in% High_SE$transcript_ID,]


# Visualize the postexercise data to see if same pattern remains in postexercise
# Visualization
long_df_post <- post_splice_df%>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID) )%>%
  inner_join(post_metadata, by = "seq_sample_id")



# Plot the relationship between age and splicing efficiency
 long_df_post %>%
  group_by(group)%>%
  summarise(avg = mean(SE)) %>%
  ggplot(aes(group,avg))+
  geom_point(mapping = aes(colour = group, size = 10))+ 
  geom_smooth()+
  ggtitle("Relationship between age group and splicing efficiency") +
  ylab(" Average splicing efficiency at post exercise")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_cowplot()
ggsave("Figures/SE_post_exc.png", bg= "white")

long_df_post %>%
  group_by(age, sex, group)%>%
  summarise(avg = mean(SE)) %>%
  ggplot(aes(age, avg))+
  geom_point(mapping = aes(colour = sex ,  size = 5))+ 
  geom_smooth()+
  ggtitle("Relationship between age, gender and splicing efficiency") +
  ylab(" Average splicing efficiency at postexercise")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_cowplot()
ggsave("Figures/SE_post_exc_age.png", bg= "white")
# ggarrange(x, y, 
#           legend = "top")

