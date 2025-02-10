# This script explores the splicing data and metadata 

library(dplyr)
library(ggplot2)
library(biomaRt)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)
library(cowplot)
library(gridExtra)
library(grid)
library(gt)


source("R/Trainome_functions.R")
# metadata
#Load the metadata of the the datasets

#Copd 
copd_metadata <- readRDS("data_new/Pre_Exercise/copd_preExc_metadata.RDS") %>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)




#Volume

volume_metadata <- readRDS("data_new/Pre_Exercise/vol_preExc_metadata.RDS")%>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)


# Contratratrain
Contratrain_metadata <- readRDS("data_new/Pre_Exercise/ct_PreExc_metadata.RDS")%>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)


# SRP102542
SRP102542_metadata <- readRDS("data_new/Pre_Exercise/SRP102542_preExc_metadata.RDS")%>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)



# Alpha_and_Omega
Alpha_Omega_metadata <- readRDS("data_new/Pre_Exercise/Alpha_Omega_PreExc_metadata.RDS")
length(unique(Alpha_Omega_metadata$participant))
length(unique(Alpha_Omega_metadata$seq_sample_id))
# Relief

Relief_metadata <- readRDS("data_new/Pre_Exercise/Relief_PreExc_metadata.RDS")


# Merge all in one
all_pre_metadata <- rbind(copd_metadata, volume_metadata)%>%
  rbind(Contratrain_metadata) %>%
  rbind(Alpha_Omega_metadata) %>%
  rbind(SRP102542_metadata) %>%
  rbind(Relief_metadata) %>%
  #Copd and volume and SRP102542 have age data in decimal
  mutate(across(c("age"), round, 0)) %>%
  # Create age groups where we have those above 50, and those below 50 
  
   mutate(group = case_when(age <=50 ~ "Fifty and below" ,
                            age > 50  ~ "Above fifty")) %>%
  mutate(sex = factor(sex, levels = c("female", "male")),
          group = factor(group, levels = c("Fifty and below" ,"Above fifty" ))) 


all_pre_metadata$participant <- paste0(all_pre_metadata$study, "_", all_pre_metadata$participant)



 saveRDS(all_pre_metadata, "data_new/Pre_Exercise/all_prexercise_metadata.RDS")
 
 
 

 ggplot(all_pre_metadata, aes(age, fill = group)) +
  geom_bar()+
  ggtitle("Distribution of baseline data by group")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("Age group of participants")+
  stat_count(geom = "Text", aes(label = ..count..), vjust = 1.5)

 ggsave("Figures/Baseline_data.png")


 
 
 # check age group
 ggplot(all_pre_metadata, aes(sex, fill = group)) +
  geom_bar()+
#  ggtitle("Distribution of baseline data")+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = 1.5)


 

hist(all_pre_metadata$age)


pre_chart <- all_pre_metadata %>%
  ggplot(aes(x= age, fill = sex, colour = sex))+
  geom_histogram(binwidth = 1, colour= "lightblue") +
  ylab("Number of participants")+
  facet_wrap(study ~.)+
  ggtitle(" Distribution of baseline samples across all studies")+
     theme(plot.title = element_text(hjust = 0.5))+
     theme_cowplot() # removes grid

#ggsave("Images_tables/Distribution_baseline_samples.png", bg = "white", scale = 1.5, dpi = 400)





# Load the splicing information
copd_data <- readRDS("data_new/Pre_Exercise/copd_preExc_splicing_data.RDS")


volume_data <- readRDS("data_new/Pre_Exercise/vol_preExc_splicing_data.RDS")

contratrain_data <- readRDS("data_new/Pre_Exercise/ct_PreExc_splicing_data.RDS")

SRP102542_data <- readRDS("data_new/Pre_Exercise/SRP102542_preExc_splicing_data.RDS")

Alpha_Omega_data <- readRDS("data_new/Pre_Exercise/Alpha_Omega_PreExc_splicing_data.RDS")
Relief_data <- readRDS("data_new/Pre_Exercise/Relief_PreExc_splicing_data.RDS")

all_pre_splice <- copd_data%>%
  inner_join(volume_data, by = "transcript_ID") %>%
  inner_join(contratrain_data, by = "transcript_ID")%>%
  inner_join(SRP102542_data, by = "transcript_ID") %>%
  inner_join(Alpha_Omega_data, by = "transcript_ID") %>%
  inner_join(Relief_data, by = "transcript_ID") %>%
  drop_na()



saveRDS(all_pre_splice, "data_new/Pre_Exercise/all_pre_Exc_splicing_data.RDS")








# Are there specific introns that charactaristically have low splicing efficiency

low_SE_df <- all_pre_splice %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID)) %>%
#  inner_join(all_pre_metadata, by = "seq_sample_id") %>%
  summarise(.by = transcript_ID, 
            mode = getmode(SE),
            min = min(SE), 
            max = max(SE), 
            mean = mean(SE)
           # q20 = quantile(SE, 0.2), 
           #  range = max(SE) - min(SE)
           ) %>%
  filter(max <= 0.2)


low_SE <- low_SE_df %>%
  separate("transcript_ID", c("transcript_ID", "intron_ID", "chr"), sep = "_")



table_plot_low <- tableGrob(low_SE)

# Open a PNG device
pdf("Figures/low_se_ints.pdf")

# Draw the table
grid.draw(table_plot_low)

# Close the pdf device
dev.off()







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



# Get the ensemble annotation of genes
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "hsapiens_gene_ensembl",
                   verbose = TRUE)

attributes <- listAttributes(ensembl)

# #extract GENE biotypes and  names
annotation<- getBM(attributes = c("transcript_biotype", "external_transcript_name",
                                  "ensembl_gene_id", "ensembl_transcript_id_version", "external_gene_name", "transcript_length"),  mart = ensembl )
annotation_low_SE <- inner_join(low_SE, annotation, by= c("transcript_ID" = "ensembl_transcript_id_version"))


# save the poorly spliced introns as pdf
table_plot_low <- tableGrob(annotation_low_SE %>%
                              dplyr::select(transcript_ID, intron_ID, transcript_biotype,external_gene_name, mode, min, max, ))

# Open a PNG device
jpeg("Figures/Table1.jpeg", width = 780, height = 480, quality = 100)

# Draw the table
grid.draw(table_plot_low)

# Close the pdf device
dev.off()



# Perfectly spliced introns



annotation_high_SE <- inner_join(High_SE, annotation, by= c("transcript_ID" = "ensembl_transcript_id_version"))

high_distribution <- annotation_high_SE %>%
  ggplot(aes(transcript_biotype, fill = transcript_biotype))+
  geom_bar(width = 1, alpha =.8 )+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = 1)+
  ggtitle("Biotypes of transcripts containing the perfectly spliced introns")+
  ylab("Number of unique transcripts")


# ggsave("Images_tables/Figure_1.jpeg")
# 
# 
# table_plot_high <- tableGrob(annotation_high_SE %>%
#                                dplyr::select(transcript_ID, intron_ID, transcript_biotype,external_gene_name, mode, min, max))
# 
# # Open a PNG device
# png("Images_tables/supplementary_table1.png", width = 700)
# 
# # Draw the table
# grid.draw(table_plot_high)
# 
# # Close the pdf device
# dev.off()



length(unique(annotation_high_SE$transcript_ID))



saveRDS(annotation_low_SE, "data_new/Pre_Exercise/annotated_low_SE_introns.RDS")

saveRDS(annotation_high_SE, "data_new/Pre_Exercise/annotated_perfect_SE_introns.RDS")

saveRDS(annotation, "data_new/ensembl_gene_annotation.RDS")




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
        
        font.size = 8, title = "Enriched biological processes in perfectly spliced out  introns") +
  theme(axis.text = element_text(size = 10), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 10))

ggarrange(high_distribution, high_ontology, align = "hv")
 
 
 
 
 ggsave("Images_tables/Figure_2Ontology_high_se.png", dpi = 400, scale = 2)




# check if there is a relationship between splicing efficiency and transcript length
long_df <- all_pre_splice%>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID) )%>%
  inner_join(all_pre_metadata, by = "seq_sample_id")

df <- long_df %>%
  group_by(transcript_ID) %>%
  summarize(avg = round(mean(SE), digits = 2)) %>%
  separate("transcript_ID", c("transcript_ID", "intron_ID", "chr"), sep = "_") %>%
  inner_join(annotation, by = c("transcript_ID" = "ensembl_transcript_id_version"), copy = T) %>%
  mutate(.by = external_transcript_name, 
           SE_per_gene = mean(avg))



plot(df$SE_per_gene, df$transcript_length)

cor.test(df$SE_per_gene, df$transcript_length)




#COPD metadata
copd_metadata <- readRDS("data_new/processed_data/copd_metadata.RDS") %>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)
length(unique(copd_metadata$participant))

# Volume_data
Vol_metadata <- readRDS("data_new/processed_data/volume_metadata.RDS")%>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age) 
length(unique(Vol_metadata$participant))
unique(Vol_metadata$time)

# Contratrain_data
ct_metadata <- readRDS("data_new/processed_data/contratrain_metadata.RDS") %>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)
length(unique(ct_metadata$participant))

unique(ct_metadata$time)

# Publicly available data
SRP102542_metadata <- readRDS("data_new/processed_data/SRP102542_metadata.RDS")%>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)
unique(SRP102542_metadata$time)



# Alpha and Omega data

A_Omega_metadata <- readRDS("data_new/processed_data/Alpha_Omega_metadata.RDS")%>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)
unique(A_Omega_metadata$time)


Relief_full_meta <- readRDS("data_new/processed_data/Relief_metadata.RDS")
unique(Relief_full_meta$seq_sample_id)
# Merge them all in one
length(unique(Relief_full_meta$participant))

all_full_metadata <- rbind(copd_metadata, Vol_metadata)%>%
  rbind(ct_metadata) %>%
  rbind(A_Omega_metadata) %>%
  rbind(SRP102542_metadata) %>%
  rbind(Relief_full_meta) %>%
  
  mutate(across(c("age"), round, 0)) %>%
  mutate(group = case_when(age <=50 ~ "Fifty and below" ,
                           age > 50  ~ "Above fifty")) %>%
  mutate(sex = factor(sex, levels = c("female", "male")),
         group = factor(group, levels = c("Fifty and below" ,"Above fifty" )),
         time = factor(time, levels = c("PreExc", "PostExc")))
all_full_metadata$participant <- paste0(all_full_metadata$study, "_", all_full_metadata$participant)



saveRDS(all_full_metadata, "data_new/processed_data/all_full_metadata.RDS")


# Extract only post exercise data
post_metadata <- all_full_metadata %>%
  filter(time == "PostExc")




copd_splice_df <- readRDS("data_new/processed_data/copd_splicing_data.RDS")

vol_splice_df <- readRDS("data_new/processed_data/volume_splicing_data.RDS")
ct_splice_df <- readRDS("data_new/processed_data/contratrain_splicing_data.RDS")
SRP102542_splice_df <- readRDS("data_new/processed_data/SRP102542_splicing_data.RDS")
AOD_splice_df <- readRDS("data_new/processed_data/Alpha_Omega_splicing_data.RDS")
Relief_full_splice <- readRDS("data_new/processed_data/Relief_splicing_data")

all_splice_df <-copd_splice_df  %>%
  inner_join(vol_splice_df, by = "transcript_ID") %>%
  inner_join(ct_splice_df, by = "transcript_ID") %>%
  inner_join(AOD_splice_df, by = "transcript_ID") %>%
  inner_join(SRP102542_splice_df, by = "transcript_ID") %>%
  inner_join(Relief_full_splice, by = "transcript_ID")


# Get the post exercise splicing data
# Getting it before removing the NAs in the full data ensures more introns are captured

post_intersect <- intersect(colnames(all_splice_df),post_metadata$seq_sample_id)

post_splice_df <- all_splice_df %>%
  subset(select = c("transcript_ID", post_intersect))%>%
  drop_na()


# select only the splicing samples captured in the metadata
all_intersect <- intersect(colnames(all_splice_df), all_full_metadata$seq_sample_id)

all_splice_df <- all_splice_df %>%
  subset(select = c("transcript_ID", all_intersect))%>%
  drop_na()

saveRDS(all_splice_df, "data_new/processed_data/all_splice_data.RDS")





low_SE_full <- post_splice_df %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID)) %>%
  summarise(.by = transcript_ID, 
            mode = getmode(SE),
            min = min(SE), 
            max = max(SE), 
            mean = mean(SE)) %>%
  filter(max <= 0.2) %>%
  separate("transcript_ID", c("transcript_ID", "intron_ID", "chr"), sep = "_")






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




#a <- readr::read_tsv("data_new/Alpha_Omega_SpliceQ_outputs/s100_EKRN240058206.tsv")
# Query the low and perfect spliced introns to see how they fared postexercise

High_at_post <- post_splice_df[post_splice_df$transcript_ID %in% High_SE_df$transcript_ID,]



intersect_high_SE <- High_SE[High_SE$transcript_ID %in% High_SE_full$transcript_ID, ]





annotation_high_SE_full <- inner_join(High_SE_full, annotation, by= c("transcript_ID" = "ensembl_transcript_id_version"))

annotation_high_SE_full %>%
  ggplot(aes(transcript_biotype, fill = transcript_biotype))+
  geom_bar(width = 1, alpha =.8 )+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = 1)+
  ggtitle("Biotypes of transcripts containing the perfectly spliced introns")+
  ylab("Number of unique transcripts")



# Visualize the postexercise data to see if same pattern remains in postexercise
# Visualization
long_df_post <- all_spliceq_df %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID) )%>%
  inner_join(post_metadata, by = "seq_sample_id")



df_post <- long_df_post %>%
  group_by(transcript_ID) %>%
  summarize(avg = round(mean(SE), digits = 2)) %>%
  separate("transcript_ID", c("transcript_ID", "intron_ID", "chr"), sep = "_") %>%
  inner_join(annotation, by = c("transcript_ID" = "ensembl_transcript_id_version"), copy = T) %>%
  mutate(.by = external_transcript_name, 
         SE_per_gene = mean(avg))



plot(df_post$SE_per_gene, df_post$transcript_length)

 x <- cor.test(df_post$SE_per_gene, df_post$transcript_length) %>%
  as.dataframe() 

 
post_chart <-   all_full_metadata %>%
   ggplot(aes(x= age, fill = sex, colour = sex))+
   geom_histogram(binwidth = 1, colour= "lightblue") +
   ylab("Number of participants")+
   facet_wrap(study ~.)+
   ggtitle(" Distribution of postexercise samples across all studies")+
   theme(plot.title = element_text(hjust = 0.5))+
   theme_cowplot()

ggarrange(pre_chart, post_chart, 
          labels = "A", "B")

ggsave("Images_tables/Distribution_samples.png", bg = "white", scale = 7, dpi = 400, 
       limitsize = F)
