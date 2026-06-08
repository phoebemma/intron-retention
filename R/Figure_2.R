# Load the file with the model and data
source("R/Load_data_for_visualisation.R")

library(patchwork)
library(ggrepel)






binom_df <- binom_model$summaries %>%
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



top10_labels <- binom_df %>%
  filter(sig) %>%
  group_by(term) %>%
  slice_max(abs(estimate), n = 15, with_ties = FALSE) %>%
  ungroup()


term_summary <- binom_df %>%
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


binom_fig <- ggplot(binom_df, aes(estimate, neg_log10_fdr, colour = effect)) +
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


## Repeat same for the beta-binomial model


beta_binom_df <- beta_binom_model$summaries %>%
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



beta_top10_labels <- beta_binom_df %>%
  filter(sig) %>%
  group_by(term) %>%
  slice_max(abs(estimate), n = 20, with_ties = FALSE) %>%
  ungroup()


beta_term_summary <- beta_binom_df %>%
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


beta_binom_fig <- ggplot(beta_binom_df, aes(estimate, neg_log10_fdr, colour = effect)) +
  geom_point(aes(colour = effect), alpha = 0.7, size = 2) +
  
  # geom_text_repel(
  #   data = beta_top10_labels,
  #   aes(label = gene_intron),
  #   size = 4,
  #   max.overlaps = Inf
  # ) +
  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  
  scale_color_manual(values = c("Improved SE" =  colors[6], "Reduced SE" = colors[1]),
                     name = "Effect") +
  
  facet_wrap(~term, scales = "fixed") +
  coord_cartesian(xlim = c(-1, 1))+
  
  labs(
    title = "Differentially spliced introns due to Aging and Resistance Training",
    subtitle = "Beta-binomial model (splicing efficiency coded as 0,1)",
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



# Plot the introns differentially spliced by both aging and resistance training

shared <- beta_binom_df %>%
  group_by(gene_intron) %>%
  filter(n_distinct(term) == 2) %>%   # keeps only those seen in both terms
  ungroup()

shared_sig_introns <- shared %>%
  dplyr::group_by(gene_intron) %>%
  dplyr::filter(all(adj.p <= 0.05)) %>%   # significant in BOTH terms
  dplyr::ungroup()
  



both <- ggplot(shared_sig_introns, aes(x = estimate, y = gene_intron, color = effect)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ term, scales = "fixed") +
  scale_color_manual(values = c("Improved SE" = colors[6],
                                "Reduced SE" = colors[1]),
                     name = "Effect") +
  labs(
    x = "Effect size",
    y = NULL,
    title = "Introns Affected by Both Age and Training"
  ) +
  theme_minimal() +
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



Figure_1 <- binom_fig + beta_binom_fig 

Figure_1 +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(
    ,
    plot.tag = element_text(size = 14, face = "bold"),
    plot.tag.position = c(0.08, 0.98)
  ) &
  guides(color = guide_legend(title = "Effect"))

