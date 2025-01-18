library(dplyr)
library(ggplot2)



# Load the primary pre-exercise model
# This was built using the formula y ~  scaled_age + (1|study) + (scaled_age+0|study) +(1|participant)
Pre_group <- readRDS("data_new/models/scaled_age_seperate_slope_intercept_model.RDS") # %>%
 # drop_na()
colnames(Pre_group)
unique(Pre_group$coef)

filt_pre_group <- Pre_group %>%
  filter(Pr...z.. <= 0.05 & fcthreshold == "s")



filt_pre_group %>%
  ggplot(aes(x = coef)) +
  geom_bar()+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))



unique(filt_pre_group$coef)
length(unique(filt_pre_group$target))


 # Appears to be the best model
saveRDS(filt_pre_group, "data_new/models/filt_scaled_age_seperate_slope_intercept_model.RDS")



# explore the annotations of the ds introns

# Load the gene annotation file
 gene_annotation <- readRDS("data_new/ensembl_gene_annotation.RDS")



# Extract the transcript_Id of the ds introns
 
 filt <- filt_pre_group %>%
   separate("target", c("transcript_ID", "intron_ID", "chr"), sep = "_") %>%
    inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version")) %>%
   dplyr::select(Estimate, transcript_ID, transcript_biotype, ensembl_gene_id, external_gene_name, transcript_length)

 hist(filt$Estimate)
 # How many unique transcripts
length(unique(filt$transcript_ID))

# How many unique genes
length(unique(filt$ensembl_gene_id))

filt %>%
   ggplot(aes(transcript_biotype))+
   geom_bar()+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
   stat_count(geom = "Text", aes(label = ..count..), vjust = -0.5) +
  ggtitle("Annotation of genes containing differentially spliced introns at baseline") +
  ylab("Number of genes")

cor(filt$Estimate, filt$transcript_length)

cor.test(filt$Estimate, filt$transcript_length)


