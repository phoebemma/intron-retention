library(dplyr)
library(tidyverse)
library(seqwrap)



Relief_full_meta <- readRDS("data/Relief_metadata.RDS")%>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age) %>%
  mutate(sex = factor(sex, levels = c("female", "male")),
         time = factor(time, levels = c("PreExc", "PostExc")),
         age_group = case_when(age < 40 ~ "Young",
                               age > 40 ~ "Old"),
         age_group = factor(age_group, levels = c("Young", "Old")))


Relief_full_splice <- readRDS("data/Relief_splicing_data.RDS")  %>%
  drop_na()



# REORDER THE SEQUENCE ID TO MATCH BOTH DATAFRAMMES
all_splice_reordered <- Relief_full_splice[, c("transcript_ID",Relief_full_meta$seq_sample_id)] 

# Check if everything matches except the transcript_id
match(colnames(all_splice_reordered), Relief_full_meta$seq_sample_id)


# visualise the data

# Build binomial model
# This model investigates the question, "given an intron,
# what is the probability of perfect splicing as a function of age and resistance exercise training"

# derive a matrix that indicates 0 if SE is not 1
one_inflated_mat <- all_splice_reordered

one_inflated_mat[-1] <- lapply(
  one_inflated_mat[-1],
  function(x) as.integer(x == 1)
)




# Intialise argument
args_binom <- list( formula = y ~ age_group + time + sex +
                      (1 | participant), family  = binomial)

# containerise using seqwrap_compose
binom <- seqwrap_compose(data       = one_inflated_mat,
                         metadata   = Relief_full_meta,
                         samplename = "seq_sample_id",
                         modelfun   = glmmTMB::glmmTMB,
                         arguments  = args_binom)

# build model
binom_results <- seqwrap(binom,
                         return_models = FALSE,
                         cores = 10)

saveRDS(binom_results, "data/Relief_binom_model.RDS")

#binom_results <- readRDS("data/Relief_binom_model.RDS") 


# The second model
# This model accepts as input the full spectrum of SE values. 
# It investigates the impact of resistance training and aging 
# on the slightest SE variations of introns.

# convert the 1.0 to 0.999. This is becasue beta-model accepts only values between 0 and one
all_splice_reordered[all_splice_reordered == 1 ] <- 0.999



# initialise the argument. This time we check the interaction of age and time
args_full <-list(formula = y ~ age_group + time + sex +
                   (1 | participant), 
                 family = glmmTMB::beta_family(link = "logit"))




# check the functions and datasets
container <- seqwrap_compose(data = all_splice_reordered,
                             metadata = Relief_full_meta,
                             samplename = "seq_sample_id",
                             modelfun = glmmTMB::glmmTMB,
                             arguments = args_full)


# build model
full_model <- seqwrap(container,
                      # summary_fun = sum_with_pred,
                      #eval_fun = eval_mod,
                      return_models = F,
                      # subset = 1:150,
                      cores = 2)

saveRDS(full_model, "data/Relief_full_model.RDS")

# full_model<- readRDS("data/Relief_full_model.RDS")



# extract the model summaries in the beta binomial model
Relief_binom <- seqwrap_summarise(binom_results)


# filter significantly differantially spliced introns
Relief_binom_outputs <- Relief_binom$summaries %>% 
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
  separate(target, into = c(NA, "intron_ID", NA), sep = "_") %>%
  # create a gene and intron label using the gene name and intron number
  mutate(
    gene_label = ifelse(
      is.na(external_gene_name) | external_gene_name == "",
      ensembl_gene_id,
      external_gene_name
    ),
    gene_intron = paste(gene_label, intron_ID, sep = " : ")
  ) %>%
  arrange(gene_label, estimate) %>%
  mutate(gene_intron = factor(gene_intron, levels = unique(gene_intron))) %>%
  dplyr::select(intron_ID, term, effect, adj.p, estimate, p.value, intron_length,
                sig, neg_log10_fdr, gene_intron, gene_label, transcript_biotype)



top10_labels <- Relief_binom_outputs %>%
  filter(sig) %>%
  group_by(term) %>%
  slice_max(abs(estimate), n = 6, with_ties = FALSE) %>%
  ungroup()


term_summary <- Relief_binom_outputs %>%
  group_by(term) %>%
  summarise(
    n_total = n(),
    n_sig = sum(sig),
    perc_sig = 100 * n_sig / n_total,
    x = max(estimate, na.rm = TRUE),
    y = max(neg_log10_fdr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      "ds introns: ",
      n_sig, "/", n_total,
      " (", round(perc_sig, 1), "%)"
    )
  )

 ggplot(Relief_binom_outputs, aes(estimate, neg_log10_fdr, colour = effect)) +
  geom_point(aes(colour = effect), alpha = 0.7, size = 2) +
  
  geom_text_repel(
    data = top10_labels,
    aes(label = gene_intron),
    size = 4,
    max.overlaps = Inf
  ) +
  
  # geom_text(
  #   data = term_summary,
  #   aes(x = -Inf, y = Inf, label = label),
  #   inherit.aes = FALSE,
  #   hjust = -0.1,
  #   vjust = 1.2,
  #   size = 3.5,
  #   fontface = "bold"
  # ) +
  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  
  scale_color_manual(values = c("Improved SE" =  colors[6], "Reduced SE" = colors[1]),
                     name = "Effect") +
  
  facet_wrap(~term, scales = "fixed") +
  coord_cartesian(xlim = c(-3, 3))+
  
  labs(
    title = "Differentially spliced introns due to Aging and Resistance Training",
    subtitle = "Binomial model (splicing efficiency coded as 0/1)",
    x = "Effect size",
    y = expression(-log[10]("FDR value"), clip = "off")
  ) +
  
  # theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 13),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, face = "bold"), 
    
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"), 
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
    
    
  )



# extract the results of the non-binarised model

full_model_sum <- seqwrap_summarise(full_model)

# The non-binarised model 
Relief_beta_binom_outputs <- full_model_sum$summaries %>% 
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
  separate(target, into = c(NA, "intron_ID", NA), sep = "_") %>%
  # create a gene and intron label using the gene name and intron number
  mutate(
    gene_label = ifelse(
      is.na(external_gene_name) | external_gene_name == "",
      ensembl_gene_id,
      external_gene_name
    ),
    gene_intron = paste(gene_label, intron_ID, sep = " : ")
  ) %>%
  arrange(gene_label, estimate) %>%
  mutate(gene_intron = factor(gene_intron, levels = unique(gene_intron))) %>%
  dplyr::select(intron_ID, term, effect, adj.p, estimate, p.value, intron_length,
                sig, neg_log10_fdr, gene_intron, gene_label, transcript_biotype)

relief_x <- Relief_beta_binom_outputs %>%
  filter(adj.p <= 0.05)
ggplot(relief_x, aes(x = estimate, y = gene_intron, color = effect)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ term, scales = "free_y") +
  scale_color_manual(values = c("Improved SE" = colors[6],
                                "Reduced SE" = colors[1]),
                     name = "Effect") +
  labs(
    x = "Effect size",
    y = NULL,
    title = "DS introns due to Aging and Resistance Training (RT)",
    subtitle = "Beta-binomial model (0,1)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 8),
    strip.text = element_text(face= "bold")
  )



relief_top10_labels <- Relief_beta_binom_outputs %>%
  filter(sig) %>%
  group_by(term) %>%
  slice_max(abs(estimate), n = 15, with_ties = FALSE) %>%
  ungroup()


relief_term_summary <- Relief_beta_binom_outputs %>%
  group_by(term) %>%
  summarise(
    n_total = n(),
    n_sig = sum(sig),
    perc_sig = 100 * n_sig / n_total,
    x = max(estimate, na.rm = TRUE),
    y = max(neg_log10_fdr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      "ds introns: ",
      n_sig, "/", n_total,
      " (", round(perc_sig, 1), "%)"
    )
  )

ggplot(Relief_beta_binom_outputs, aes(estimate, neg_log10_fdr, colour = effect)) +
  geom_point(aes(colour = effect), alpha = 0.7, size = 2) +
  
  geom_text_repel(
    data = relief_top10_labels,
    aes(label = gene_intron),
    size = 4,
    max.overlaps = Inf
  ) +
  
  # geom_text(
  #   data = term_summary,
  #   aes(x = -Inf, y = Inf, label = label),
  #   inherit.aes = FALSE,
  #   hjust = -0.1,
  #   vjust = 1.2,
  #   size = 3.5,
  #   fontface = "bold"
  # ) +
  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  
  scale_color_manual(values = c("Improved SE" =  colors[6], "Reduced SE" = colors[1]),
                     name = "Effect") +
  
  facet_wrap(~term, scales = "fixed") +
  coord_cartesian(xlim = c(-3, 3))+
  
  labs(
    title = "Differentially spliced introns due to Aging and Resistance Training",
    subtitle = "Binomial model (splicing efficiency coded as 0/1)",
    x = "Effect size",
    y = expression(-log[10]("FDR value"), clip = "off")
  ) +
  
  # theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 13),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, face = "bold"), 
    
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"), 
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
    
    
  )


