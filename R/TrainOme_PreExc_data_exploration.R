library(dplyr)
library(ggplot2)
library(biomaRt)
library(tmod)

# Load the primary pre-exercise model

prim_Pre_group <- readRDS("data/Trainome_data_models/preExc_model_age_group.RDS")
unique(prim_Pre_group$coef)

# Get those differentially spliced by age

ds_age <- prim_Pre_group %>%
  subset(coef == "age_groupOld")%>%
  dplyr::select(target,Estimate, Std..Error, Pr...z.., adj.p,  log2fc ) %>%
  separate("target", c("transcript_ID", "intron_ID", "chr"), sep = "_") 
hist(ds_age$adj.p)
saveRDS(ds_age, "data/Trainome_data_models/ds_age_group.RDS")

hist(ds_age$Estimate, col = "darkgrey", border = "black", xlab = "Estimate of DS introns by age", main = "Distribution of model Estimate in DS introns by age group" )

length(unique(ds_age$transcript_ID))
# extract those with at least logf2c
sig_ds <- prim_Pre_group %>%
  subset(coef == "age_groupOld")%>%
  filter(fcthreshold == "s")%>%
  dplyr::select(target,Estimate, Std..Error, Pr...z.., adj.p,  log2fc ) %>%
  separate("target", c("transcript_ID", "intron_ID", "chr"), sep = "_") 
saveRDS(sig_ds, "data/Trainome_data_models/sig_ds_age_group.RDS")

hist(sig_ds$Estimate, col = "darkgrey", border = "black", xlab = "Estimate of DS introns by age" , 
     main = "Distribution of model Estimate in DS introns(age group) with log2fc" )

length(unique(sig_ds$transcript_ID))

ds_sex <- prim_Pre_group %>%
  subset(coef == "sexmale")%>%
  dplyr::select(target,Estimate, Std..Error, Pr...z.., adj.p,  log2fc ) %>%
  separate("target", c("transcript_ID", "intron_ID", "chr"), sep = "_") 
saveRDS(ds_sex, "data/Trainome_data_models/ds_sex.RDS")


hist(ds_sex$Estimate,col = "darkgrey", border = "black", xlab = "Estimate of DS introns by sex" , 
     main = "Distribution of model Estimate in DS introns  by sex")

ds_sex_age <- prim_Pre_group %>%
  subset(coef == "age_groupOld:sexmale")%>%
  dplyr::select(target,Estimate, Std..Error, Pr...z.., adj.p,  log2fc ) %>%
  separate("target", c("transcript_ID", "intron_ID", "chr"), sep = "_") 
saveRDS(ds_sex_age, "data/Trainome_data_models/ds_age_group_sex.RDS")

# Load the ds introns by age as continous variable
 ds_age_count <- readRDS("data/Trainome_data_models/preExc_model_age.RDS")%>%
   dplyr::select(target,Estimate, Std..Error, Pr...z.., adj.p,  log2fc ) %>%
   separate("target", c("transcript_ID", "intron_ID", "chr"), sep = "_") 
 
 hist(ds_sex_age$Estimate,col = "darkgrey", border = "black", xlab = "Estimate of DS introns by sex and age" , 
      main = "Distribution of model Estimate in DS introns  by sex and age")
 
 # Any intersect between the ds in age_group and age as a continous variable ?
 
 intersect(ds_age$target, ds_age_count$target)
 
 
 
 #Use the ensemble database to get the annotation of the transcripts
 # The vversion should correspond to the annotation version of the reference genome used
 ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                    dataset = "hsapiens_gene_ensembl",
                    host = "https://apr2022.archive.ensembl.org", 
                    verbose = TRUE)
 
 attributes <- listAttributes(ensembl)
 
 
 annotation<- getBM(attributes = c("transcript_biotype", "external_transcript_name",
                                   "ensembl_gene_id", "ensembl_transcript_id_version", "external_gene_name"),  mart = ensembl )

 # Save annotationso we would not need to run annotation again
 
 saveRDS(annotation, "data/ENSEMBL_gene_annotation_april_2022_version.RDS")
 
 
 
 annotation_age_group <- inner_join(ds_age, annotation, by= c("transcript_ID" = "ensembl_transcript_id_version"))
annotation_age_group %>%
  ggplot(aes(transcript_biotype, fill = transcript_biotype))+
   geom_bar()+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 18, , face = "bold"),
         axis.text.y = element_text( size = 15),
         plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))+
   stat_count(geom = "Text", aes(label = ..count..), vjust = -0.5)+
  ggtitle("Annotation of genes containing DS introns by age group")

saveRDS(annotation_age_group, "data/Trainome_data_models/annotation_ds_age_group.RDS")

annotation_sex <- inner_join(ds_sex, annotation, by= c("transcript_ID" = "ensembl_transcript_id_version"))

annotation_sex %>%
  ggplot(aes(transcript_biotype, fill = transcript_biotype))+
  geom_bar(width = 0.2)+
  theme(axis.text.x = element_text(angle = 90,  size = 18, , face = "bold"),
        axis.text.y = element_text( size = 15),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = -0.5)+
  ggtitle("Annotation of genes containing DS introns by sex")

saveRDS(annotation_sex, "data/Trainome_data_models/annotation_ds_sex.RDS")
