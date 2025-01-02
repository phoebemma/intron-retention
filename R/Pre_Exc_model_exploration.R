library(dplyr)
library(ggplot2)
library(biomaRt)
library(tmod)


# Load the primary pre-exercise model

Pre_group <- readRDS("data_new/models/preExc_group_only_model.RDS") # %>%
 # drop_na()
colnames(Pre_group)
unique(Pre_group$coef)

filt_pre_group <- Pre_group %>%
  filter(coef %in% c("group>20 & <=30",
                     "group>30 & <=40",
                     "group>40 & <=50",
                     "group>50 & <=60",
                     "group>60 & <=70",
                     "group>70" )) %>%
  # extract the adjusted p value, log2fc and fcthreshold of each coefficient
  mutate(.by = coef,
         adj.p = p.adjust(Pr...z.., method = "fdr"),
                  log2fc = Estimate/log(2),
                  fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns"))%>%
  filter(adj.p <= 0.05)


unique(filt_pre_group$coef)





# Extract ds introns in the model that looked at age group and sex


pre_group_sex <- readRDS("data_new/models/preExc_group_and_sex_model.RDS")

unique(pre_group_sex$coef)

filt_pre_grop_sex <- pre_group_sex %>%
  filter(coef %in% c("group>20 & <=30",
                     "group>30 & <=40",  "group>40 & <=50",
                     "group>50 & <=60",
                     "group>60 & <=70",
                     "group>70",
                     "sexmale",
                     "group>20 & <=30:sexmale",
                     "group>30 & <=40:sexmale",
                     "group>40 & <=50:sexmale",
                     "group>50 & <=60:sexmale",
                     "group>60 & <=70:sexmale",
                     "group>70:sexmale" )) %>%
mutate(.by = coef,
       adj.p = p.adjust(Pr...z.., method = "fdr"),
       log2fc = Estimate/log(2),
       fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns"))%>%
  filter(adj.p <= 0.05)


unique(filt_pre_grop_sex$coef)







pre_group_sex_no_int <- readRDS("data_new/models/preExc_group_and_sex_without_interaction_model.RDS")

unique(pre_group_sex_no_int$coef)

filt_pre_grop_sex_no_int <- pre_group_sex_no_int %>%
  filter(coef %in% c("group>20 & <=30",
                     "group>30 & <=40", 
                     "group>40 & <=50",
                     "group>50 & <=60",
                     "group>60 & <=70",
                     "group>70",
                     "sexmale"
                     )) %>%
  mutate(.by = coef,
         adj.p = p.adjust(Pr...z.., method = "fdr"),
         log2fc = Estimate/log(2),
         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns"))%>%
  filter(adj.p <= 0.05)

saveRDS(filt_pre_grop_sex_no_int, "data_new/models/filt_preExc_group_sex_no_int.RDS")

length(unique(filt_pre_grop_sex$target))

length(unique(filt_pre_group$target))


length(unique(filt_pre_grop_sex_no_int$target))





# ds_20_30 <- Pre_group %>%
#   subset(coef == "group>30 & <=40:sexmale")%>%
#   mutate(adj.p = p.adjust(Pr...z.., method = "fdr") ,
#          log2fc = Estimate/log(2),
#          fcthreshold = if_else(abs(log2fc) > 1, "s", "ns")) %>%
#   filter(adj.p <= 0.05 )
# # 1801 introns are differentially spliced by age
# saveRDS(ds_age, "data/re_models/primary_model_extracts/diff_spliced_ints_by_age_group.RDS")
# 
# # Those with an absolute log2fc of 1
# ds_age_sig <- ds_age %>%
#   filter(fcthreshold == "s")
# 
# saveRDS(ds_age_sig, "data/re_models/primary_model_extracts/sig_DS_introns_ageGroup.RDS")
# # 
# hist(ds_age$Estimate)
# hist(ds_age_sig$Estimate)
# hist(ds_age$Pr...z..)
# 
# hist(ds_age$adj.p)
# 
# # Get those diferentially expressed by sex
# ds_sex <- prim_Pre_group %>%
#   subset(coef == "sexmale")%>%
#   # filter(Pr...z.. <= 0.05) %>%
#   mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
#          log2fc = Estimate/log(2),
#          fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) %>%
#   filter( adj.p <= 0.05 )
# # 3 are differentially spliced by sex
# 
# saveRDS(ds_sex, "data/re_models/primary_model_extracts/DS_introns_by_sex.RDS")
# # Get the interaction between age and sex
# 
# ds_sex_age <- prim_Pre_group %>%
#   subset(coef == "age_groupOld:sexmale")%>%
#  #  filter(Pr...z.. <= 0.05) %>%
#   mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
#          log2fc = Estimate/log(2),
#          fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) %>%
#   filter( adj.p <= 0.05 )
# 
# 
# 
# 
# # Checkthe counts model
# prim_pre_count <- readRDS("data/re_models/primary_preExc_count_interaction_model.RDS")
# 
# # This suggests those that are increased with each increase in age
# age_ds_count <- prim_pre_count %>%
#   subset(coef == "age")%>%
#   mutate(adj.p = p.adjust(Pr...z.., method = "fdr") ,
#          log2fc = Estimate/log(2),
#          fcthreshold = if_else(abs(log2fc) > 1, "s", "ns")) %>%
#   filter(adj.p <= 0.05 )
# 
# 
# 
# # Any intersect between the ds in age_group and age as a continous variable ?
# 
# intersect(ds_age$target, age_ds_count$target)
# # Using the filter characteristics, only age showed differentialy splice introns
# 
# # Exploring this further 
# age_df <- ds_age %>%
#   dplyr::select(target,Estimate, Std..Error, Pr...z.., adj.p,  log2fc ) %>%
#   separate("target", c("transcript_ID", "intron_ID", "chr"), sep = "_") 
# 
# 
# age_df_sig <- ds_age_sig%>%
#   dplyr::select(target,Estimate, Std..Error, Pr...z.., adj.p,  log2fc ) %>%
#   separate("target", c("transcript_ID", "intron_ID", "chr"), sep = "_") 
# 
# age_count <- age_ds_count %>%
#   dplyr::select(target,Estimate, Std..Error, Pr...z.., adj.p,  log2fc ) %>%
#   separate("target", c("transcript_ID", "intron_ID", "chr"), sep = "_") 
# 
# 
# 
# sex_ds <- ds_sex %>%
#   dplyr::select(target,Estimate, Std..Error, Pr...z.., adj.p,  log2fc) %>%
#   separate("target", c("transcript_ID", "intron_ID", "chr"), sep = "_")
# #unique(pre_df$coef)
# ggplot(age_df, aes(chr)) +
#   geom_bar()+
#   labs(x = "Chromosomes")+
#   ggtitle("chromosomes containing DR introns in baseline model")+
#   theme(plot.title = element_text(hjust = 0.5))
# 
# ggplot(age_df_sig, aes(chr)) +
#   geom_bar()+
#   labs(x = "Chromosomes")+
#   ggtitle("chromosomes containing DR introns in baseline model")+
#   theme(plot.title = element_text(hjust = 0.5))
# 
# 
# #Use the ensemble database to get the annotation of the transcripts
# # The vversion should correspond to the annotation version of the reference genome used
# ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
#                    dataset = "hsapiens_gene_ensembl",
#                    host = "https://apr2022.archive.ensembl.org", 
#                    verbose = TRUE)
# 
# attributes <- listAttributes(ensembl)
# 
# # #extract GENE biotypes and  names
# annotation<- getBM(attributes = c("transcript_biotype", "external_transcript_name",
#                                   "ensembl_gene_id", "ensembl_transcript_id_version", "external_gene_name"),  mart = ensembl )
# annotation_all <- inner_join(age_df, annotation, by= c("transcript_ID" = "ensembl_transcript_id_version"))
# 
# annotation_sig <- inner_join(age_df_sig, annotation, by= c("transcript_ID" = "ensembl_transcript_id_version"))
# annotation_count <-  inner_join(age_count, annotation, by= c("transcript_ID" = "ensembl_transcript_id_version"))
# annotation_sex <- inner_join(sex_ds, annotation, by= c("transcript_ID" = "ensembl_transcript_id_version"))
# annotation_all %>%
#   ggplot(aes(transcript_biotype, fill = transcript_biotype))+
#   geom_bar()+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
#   stat_count(geom = "Text", aes(label = ..count..), vjust = -0.5)
# saveRDS(annotation_all, "data/re_models/primary_model_extracts/annotated_DS_introns_age.RDS")
# 
# annotation_sig %>%
#   ggplot(aes(transcript_biotype, fill = transcript_biotype))+
#   geom_bar()+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
#   stat_count(geom = "Text", aes(label = ..count..), vjust = -0.5)
# 
# saveRDS(annotation_sig, "data/re_models/primary_model_extracts/annotated_sig_DS_introns_age.RDS")
# 
# annotation_count %>%
#   ggplot(aes(transcript_biotype, fill = transcript_biotype))+
#   geom_bar()+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
#   stat_count(geom = "Text", aes(label = ..count..), vjust = -0.5)
# 
# saveRDS(annotation_count, "data/re_models/primary_model_extracts/annotated_DS_introns_age_continous.RDS")
# # How many unique transcripts
# annotation_sex %>%
#   ggplot(aes(transcript_biotype, fill = transcript_biotype))+
#   geom_bar()+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
#   stat_count(geom = "Text", aes(label = ..count..), vjust = -0.5)
# 
# length(unique(annotation$transcript_ID))
# 
