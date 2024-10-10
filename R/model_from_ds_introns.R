# This script contains codes that builds a model from differentially spliced introns at baseline

library(dplyr)
library(tidyverse)
library(seqwrap)
library(gridExtra)
library(ggpubr)
library(biomaRt)

#COPD metadata
copd_metadata <- readRDS("data/processed_data/copd_metadata.RDS")
# hist(copd_metadata$age)
# max(copd_metadata$age)
# unique(copd_metadata$time)
# 
# 
# Volume_data
Vol_metadata <- readRDS("data/processed_data/volume_metadata.RDS")

# unique(Vol_metadata$time)


ct_metadata <- readRDS("data/processed_data/contratrain_metadata.RDS")

all_full_metadata <- rbind(copd_metadata, Vol_metadata)%>%
  rbind(ct_metadata) %>%
  #Copd and volume has data in decimal
  mutate(across(c("age"), round, 0)) %>%
  mutate(age_group = case_when(study == "copd" ~ "Old", 
                               study == "vol" ~ "Young",
                               study == "ct" ~ "Young")) %>%
  mutate(sex = factor(sex, levels = c("female", "male")),
         age_group = factor(age_group, levels= c("Young", "Old")),
         time = factor(time, levels = c("PreExc", "PostExc"))) 

# Load splicing data

copd_splice_df <- readRDS("data/processed_data/copd_splicing_data.RDS")

vol_splice_df <- readRDS("data/processed_data/volume_splicing_data.RDS")
ct_splice_df <- readRDS("data/processed_data/contratrain_splicing_data.RDS")


all_splice_df <- copd_splice_df%>%
  inner_join(vol_splice_df, by = "transcript_ID") %>%
  inner_join(ct_splice_df, by = "transcript_ID") %>%
  drop_na()


#select only the splicing samples captured in the metadata
splice_intersect <- intersect(colnames(all_splice_df), all_full_metadata$seq_sample_id)

all_splice_df <- all_splice_df %>%
  subset(select = c("transcript_ID", splice_intersect))
# Relace all the 1s in the dataframe to 0.99

all_splice_df[all_splice_df == 1 ] <- 0.999


# For the model, we would use only the ds intorns
ds_introns <- readRDS("data/re_models/primary_model_extracts/diff_spliced_ints_by_age_group.RDS")
introns_of_int<- all_splice_df[all_splice_df$transcript_ID %in% ds_introns$target,]

# 108 introns are not captured
# Probably due to missing values


args<- list(formula = y ~  time + (1|study) +(1|participant), 
             #ziformula = ~1,
             family = glmmTMB::beta_family())


full_group_model <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                            arguments = args,
                            data = introns_of_int,
                            metadata = all_full_metadata,
                            samplename = "seq_sample_id",
                            summary_fun = sum_fun,
                            eval_fun = eval_mod,
                            exported = list(),
                            save_models = FALSE,
                            return_models = FALSE,
                            cores = ncores-2)

excl_group <- names(which(full_group_model$summaries == "NULL"))

geneids_group <- names(which(full_group_model$summaries != "NULL"))


mod_sum_group <- bind_rows(within(full_group_model$summaries, rm(excl_group))) %>%
  mutate(target = rep(geneids_group, each = 2))  %>%
subset(coef != "(Intercept)")%>%
mutate(adj.p = p.adjust(Pr...z.., method = "fdr")) %>%
  filter(adj.p <= 0.05 )



mod_eval_group <- bind_rows(within(full_group_model$evaluations, rm(excl_group)))%>%
  mutate(target = geneids_group)

#merge the model evaluation and summary dataframes

model_full_group <- mod_sum_group %>%
  inner_join(mod_eval_group, by = "target")


saveRDS(model_full_group, "data/re_models/primary_model_extracts/filtered_model_from_ds_intronsAge_at_baseline.RDS")



hist(model_full_group$Estimate)

# model_full_group <- readRDS("data/re_models/primary_model_extracts/filtered_model_from_ds_intronsAge_at_baseline.RDS")

ds_time <- model_full_group%>%
  dplyr::select(target,Estimate, Std..Error, Pr...z.., adj.p) %>%
  separate("target", c("transcript_ID", "intron_ID", "chr"), sep = "_") 
#Use the ensemble database to get the annotation of the transcripts
listMarts()
ensembl=useMart("ENSEMBL_MART_ENSEMBL")
datasets <- listDatasets(ensembl)
ensembl=useDataset("hsapiens_gene_ensembl",mart=ensembl)
attributes <- listAttributes(ensembl)

# #extract GENE biotypes and  names
annotation<- getBM(attributes = c("transcript_biotype", "external_transcript_name",
                                  "ensembl_gene_id", "ensembl_transcript_id_version"),  mart = ensembl )

improved <- inner_join(ds_time, annotation, by= c("transcript_ID" = "ensembl_transcript_id_version"))

improved %>%
  ggplot(aes(transcript_biotype, fill = transcript_biotype))+
  geom_bar()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = -0.5)
saveRDS(improved, "data/re_models/primary_model_extracts/annotated_improved_ds_introns_by_RT.RDS")
