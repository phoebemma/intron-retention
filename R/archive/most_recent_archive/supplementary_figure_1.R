# load the file with the data
source("R/Load_data_for_visualisation.R")

library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggrepel)




# extract summaries from beta-binomial model 
full_model_outputs <- beta_binom_model$summaries %>% 
  dplyr::select(-group) %>%
  inner_join(intron_length, by = c("target" = "transcript_ID")) %>%
  filter(term != "(Intercept)", term != "sexmale") %>%
  drop_na() %>%
  group_by(term) %>%
  mutate(
    adj.p = p.adjust(p.value, method = "fdr"),
    term = recode(term,
                  "scaled_age" = "Aging",
                  "timePostExc" = "Resistance Training"),
    effect = case_when(estimate > 0 & adj.p <= 0.05 ~ "Improved SE", 
                       estimate < 0 & adj.p <= 0.05 ~ "Reduced SE" ,
                       estimate < 0 & adj.p > 0.05 ~ "No effect",
                       estimate > 0 & adj.p > 0.05 ~ "No effect"),
    transcript_ID = str_split(target, "_",simplify= T) [,1]) %>%
  ungroup() %>%
  mutate(
    sig = adj.p <= 0.05,
    neg_log10_fdr = -log10(adj.p)
  ) %>%
  inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version")) %>%
  separate(target, into = c(NA, "intron_ID", NA), sep = "_", remove = F) %>%
  mutate(
    gene_label = ifelse(
      is.na(external_gene_name) | external_gene_name == "",
      ensembl_gene_id,
      external_gene_name
    ), # If gene-name isnt available, use ensembl_gene_id
    # simplify visualisation by generating gene_intro
    # this uses the gene name followed by semicolon and intron_id number
    # it makes identifying it easier in charts
    gene_intron = paste(gene_label, intron_ID, sep = " : "))


# extract model summary from binomial model
binom_model_outputs <- binom_model$summaries %>% 
  dplyr::select(-group) %>%
  inner_join(intron_length, by = c("target" = "transcript_ID")) %>%
  filter(term != "(Intercept)", term != "sexmale") %>%
  drop_na() %>%
  group_by(term) %>%
  mutate(
    adj.p = p.adjust(p.value, method = "fdr"),
    term = recode(term,
                  "scaled_age" = "Aging",
                  "timePostExc" = "Resistance Training"),
    effect = case_when(estimate > 0 & adj.p <= 0.05 ~ "Improved SE", 
                       estimate < 0 & adj.p <= 0.05 ~ "Reduced SE" ,
                       estimate < 0 & adj.p > 0.05 ~ "No effect",
                       estimate > 0 & adj.p > 0.05 ~ "No effect"),
    transcript_ID = str_split(target, "_",simplify= T) [,1]) %>%
  ungroup() %>%
  mutate(
    sig = adj.p <= 0.05,
    neg_log10_fdr = -log10(adj.p)
  ) %>%
  inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version")) %>%
  separate(target, into = c(NA, "intron_ID", NA), sep = "_", remove = F) %>%
  mutate(
    gene_label = ifelse(
      is.na(external_gene_name) | external_gene_name == "",
      ensembl_gene_id,
      external_gene_name
    ), # If gene-name isnt available, use ensembl_gene_id
    gene_intron = paste(gene_label, intron_ID, sep = " : "))




# plot distribution of biotypes
biotype_binom <- binom_model_outputs %>%
  filter(adj.p <= 0.05) %>%
  distinct(target, transcript_biotype, .keep_all = T) %>%
  count(transcript_biotype) %>%
  ggplot(aes(transcript_biotype, n, fill = transcript_biotype))+
  geom_col(width = 0.4)+
  geom_text(aes(label= n), vjust = -0.1, size = 5)+
  # scale_x_discrete(expand = expansion(mult = c(0.15, 0.15)))+
  # scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Biotypes of genes containing ds introns in the binomial model",
    x = "Biotype",
    y = NULL
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face= "bold", size = 14),
    axis.text.y = element_text(size = 14, face= "bold"),
    axis.text.x = element_text(size = 14, face= "bold", angle = 90, vjust = 0.5, hjust = 1),
    axis.title = element_text(size = 14, face= "bold"),
    legend.position = "none"
  )




# plot distribution of intron length
intron_length_binom <- binom_model_outputs %>%
  filter(adj.p <= 0.05) %>%
  dplyr::select(target, intron_length, effect) %>%
  ggplot(aes(x = intron_length, fill = effect)) +
  geom_histogram() +
  labs(
    title = "Distribution of intron length of ds introns in the binomial model",
    x = "Intron Length",
    y = NULL
  ) +
  # scale_x_continuous(limits = c(70, 30000)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.text = element_text(face= "bold", size = 14),
        axis.text.y = element_text(size = 14, face= "bold"),
        axis.text.x = element_text(size = 14, face= "bold", angle = 90, vjust = 0.5, hjust = 1),
        axis.title = element_text(size = 14, face= "bold"))



# Evaluating the betabinomial model

full_model <- full_model_outputs %>%
  filter(adj.p <= 0.05)



# plot distribution of phenotypes
biotype_full <- full_model_outputs %>%
  filter(adj.p <= 0.05) %>%
  distinct(target, transcript_biotype, .keep_all = T) %>%
  count(transcript_biotype) %>%
  ggplot(aes(transcript_biotype, n, fill = transcript_biotype))+
  geom_col(width = 0.4)+
  geom_text(aes(label= n), vjust = -0.1, size = 5)+
  labs(
    title = "Biotypes of genes containing ds introns in the betabinomial model",
    x = "Biotype",
    y = "Number of introns per biotype"
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.text = element_text(face= "bold", size = 14),
        axis.text.y = element_text(size = 14, face= "bold"),
        axis.text.x = element_text(size = 14, face= "bold", angle = 90, vjust = 0.5, hjust = 1),
        axis.title = element_text(size = 14, face= "bold"), 
        legend.position = "none")



# Plot the intron lengths
intron_length_full <- full_model %>%
  dplyr::select(target, intron_length, effect) %>%
  ggplot(aes(x = intron_length, fill = effect)) +
  geom_histogram() +
  labs(
    title = "Distribution of Intron Lengths of ds introns in the betabinomial model",
    x = "Intron Length",
    y = "Number of introns"
  ) +
  # scale_x_continuous(limits = c(70, 30000)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.text = element_text(face= "bold", size = 14),
        axis.text.y = element_text(size = 14, face= "bold"),
        axis.text.x = element_text(size = 14, face= "bold", angle = 90, vjust = 0.5, hjust = 1),
        axis.title = element_text(size = 14, face= "bold"),
        legend.position = "none")





sup_Fig1 <-  biotype_full + biotype_binom + intron_length_full +intron_length_binom 
sup_Fig1 +
  plot_annotation(tag_levels = "A") &
  #plot_layout(widths = c(1.1, 1))
  theme(
    plot.tag = element_text(size = 14, face = "bold")
    ,
    plot.tag.position = c(0.08, 0.98)
  ) 

# ggsave("Figures/Trainome_supp_Figure_1.png", bg = colors[4], width = 16, height = 14, dpi = 400)




