
library(dplyr)
library(tidyverse)
library(scales)
library(seqwrap)
library(glmmTMB)
library(marginaleffects)
library(purrr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(cowplot)
library(ggVennDiagram)
library(clusterProfiler)
library(org.Hs.eg.db)

# define color scale for images
## Color scale ##
colors <-  c("#d7191c", "#fdae61", "#ffffbf", "#abd9e9", "#2c7bb6", "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "grey70")




effect_colors <- c(
  "Improved SE" = colors[6],
  "Reduced SE"  = colors[1],
  "No effect"   = "grey70"
)



# ANNOTATION AND HELPER FUNCTIONS
#

gene_annotation <- readRDS("data/ensembl_gene_annotation.RDS")

gene_exp_df <- readRDS("data_new/gene_counts/batch_corrected_genecounts.RDS")

intron_length <- readr::read_tsv("data_new/Alpha_Omega_SpliceQ_outputs/A_102.tsv") %>%
  distinct(across(6:ncol(.)), .keep_all = TRUE) %>%
  mutate(
    transcript_ID  = paste0(transcript_ID, "_", intron_ID, "_", chr),
    intron_length  = abs((sj3start - sj5end) + 1),
    intron_ID      = as.character(intron_ID),
    number_introns = n()
  ) %>%
  dplyr::select(transcript_ID, intron_ID, gene_ID, intron_length, number_introns)

# Annotates outputs with gene names and intron metadata
# flip = TRUE reverses estimate from retention to SE scale
annotate_introns <- function(df, gene_annotation, intron_length, flip = TRUE) {
  df %>%
    mutate(
      estimate       = if (flip) -estimate else estimate,
      conf.low_se    = if (flip) -conf.high else conf.low,
      conf.high_se   = if (flip) -conf.low else conf.high,
      transcript_ID  = str_split(target, "_", simplify = TRUE)[, 1]
    ) %>%
    inner_join(gene_annotation,
               by = c("transcript_ID" = "ensembl_transcript_id_version")) %>%
    separate(target, into = c(NA, "intron_ID", NA),
             sep = "_", remove = FALSE) %>%
    inner_join(intron_length %>%
                 mutate(intron_ID = as.character(intron_ID)),
               by = c("intron_ID",
                      "ensembl_gene_id_version" = "gene_ID",
                      "target" = "transcript_ID")) %>%
    mutate(
      gene_label  = ifelse(
        is.na(external_gene_name) | external_gene_name == "",
        ensembl_gene_id, external_gene_name
      ),
      gene_intron = paste(gene_label, intron_ID, sep = " : ")
    ) %>%
    {if ("rank_score" %in% colnames(.)) arrange(., desc(rank_score)) else .}
}




all_splice_df <- readRDS("data/Trainome_all_splice_df.RDS")



metadata <- readRDS("data/Trainome_metadata.RDS")

# Store age range for back-transformation
age_min <- min(metadata$age, na.rm = TRUE)
age_max <- max(metadata$age, na.rm = TRUE)

# Load the gene annotation file
gene_annotation <- readRDS("data/ensembl_gene_annotation.RDS")


#  load the gene expression dataset
gene_exp_df <- readRDS("data_new/gene_counts/batch_corrected_genecounts.RDS") 




# Relief Model
Relief_model <- readRDS("data/Relief_zi_model.RDS")


# Predictions on the full data
#  Predicted SE trajectories (response scale) 
 zi_predictions <- readRDS("data/zi_predictions.RDS")
 
 
 #  Probability of perfect splicing (zi component)  
 zi_zprob <- readRDS("data/zi_zprob.RDS")
 
 
 #  Degree of retention among partially retained introns (beta component)
 zi_conditional <-   readRDS("data/zi_conditional.RDS")
 
 
 
 # Age slopes at 0.10 increments, separately per timepoint 
 
 zi_age_slopes <-   readRDS("data/zi_age_slopes.RDS")
 
 # Age slopes FDR within each timepoint ---
 zi_age_slopes_fdr <- zi_age_slopes %>%
   group_by(time) %>%
   mutate(adj.p = p.adjust(p.value, method = "fdr")) %>%
   ungroup() %>%
   mutate(
     sig           = adj.p <= 0.05,
     sig_ci        = conf.low > 0 | conf.high < 0,
     neg_log10_fdr = -log10(adj.p),
     rank_score    = -log10(adj.p) * abs(estimate),
     effect        = case_when(
       conf.high < 0 & sig ~ "Improved SE",
       conf.low  > 0 & sig ~ "Reduced SE",
       TRUE                 ~ "No effect"
     )
   ) %>%
   annotate_introns(gene_annotation, intron_length, flip = TRUE)
 
 #  Exercise effect (PostExc - PreExc) at each age anchor 
 zi_time_effects <- readRDS("data/zi_time_effects.RDS")
 
 
 # Exercise effects  FDR within each age anchor 
 zi_time_effects_fdr <- zi_time_effects %>%
   group_by(scaled_age) %>%
   mutate(adj.p = p.adjust(p.value, method = "fdr")) %>%
   ungroup() %>%
   mutate(
     real_age      = scaled_age * (age_max - age_min) + age_min,
     sig           = adj.p <= 0.05,
     sig_ci        = conf.low > 0 | conf.high < 0,
     neg_log10_fdr = -log10(adj.p),
     rank_score    = -log10(adj.p) * abs(estimate),
     effect        = case_when(
       conf.high < 0 & sig ~ "Improved SE",
       conf.low  > 0 & sig ~ "Reduced SE",
       TRUE                 ~ "No effect"
     )
   ) %>%
   annotate_introns(gene_annotation, intron_length, flip = TRUE)
 
 
 
 # Load the RELIEF model
 
 # The marginal effects predictions were built into thesummary function and thus
 # Extracting the data will follow the traditional seqwrap method
 Relief_results <- seqwrap_summarise(Relief_model)
 
 Relief_results <- Relief_results$summaries
 
 
 
 
 contrasts_of_interest <- c("train_young", "train_old", 
                            "age_effect_pre", "age_effect_post", 
                            "train_old_young")
 
 relief_contrasts <- Relief_results %>%
   filter(hypothesis %in% contrasts_of_interest) %>%
   group_by( hypothesis, component) %>%
   mutate(adj.p = p.adjust(p.value, method = "fdr"),
          component = recode(component, 
                             "response"    = "Overall retention",
                             "zprob"       = "Probability of Perfect splicing",
                             "conditional" = "Degree of retention"),
          hypothesis = recode(hypothesis, "train_young"     = "Exercise effect among young participants",
                              "train_old"       = "Exercise effect among older participants",
                              "age_effect_pre"  = "Age effect at baseline",
                              "age_effect_post" = "Age effect postexercise",
                              "train_old_young" = "Interaction effect of training and aging")) %>%
   ungroup() %>%
   mutate(
     sig           = adj.p <= 0.05,
     sig_ci        = conf.low > 0 | conf.high < 0,
     neg_log10_fdr = -log10(adj.p),
     rank_score    = -log10(adj.p) * abs(estimate),
     effect        = case_when(
       conf.high < 0 & sig ~ "Improved SE",
       conf.low  > 0 & sig ~ "Reduced SE",
       TRUE             ~ "No effect"
     )
   ) %>%
   annotate_introns(gene_annotation, intron_length, flip = TRUE) %>%
   filter(component != "Overall retention")
 
 
 #
 

 
