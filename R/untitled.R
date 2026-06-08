

df <- beta_binom_model$summaries %>%
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
  filter(adj.p <= 0.05) 
  


 # extract the splicing data
# all_splice_data <- all_splice_df %>%
#   dplyr::filter(transcript_ID %in% df$target) %>%
#   pivot_longer(names_to = "seq_sample_id",
#                values_to = "SE",
#                cols = -(transcript_ID) ) %>%
#   inner_join(df, by = c("transcript_ID" = "target"))


exp_df <- gene_exp_df %>%
  dplyr::filter(gene_name %in% df$external_gene_name) %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "gene_count",
               cols = -(gene_name) ) %>%
  inner_join(all_full_metadata, by = "seq_sample_id") %>%
  inner_join(df, by = c("gene_name" = "external_gene_name")) %>%
  group_by(gene_name, seq_sample_id) %>%
  summarise(mean_estimate = mean(estimate, na.rm = TRUE), .groups = "drop")




exp_summary <- gene_exp_df %>%
  dplyr::filter(gene_name %in% df$external_gene_name) %>%
  
  pivot_longer(
    cols = -gene_name,
    names_to = "seq_sample_id",
    values_to = "gene_count"
  ) %>%
  
  inner_join(all_full_metadata, by = "seq_sample_id") %>%
  inner_join(df, by = c("gene_name" = "external_gene_name")) %>%
  
  group_by(gene_name) %>%
  summarise(
    mean_estimate = mean(estimate, na.rm = TRUE),
    mean_expression = mean(gene_count, na.rm = TRUE),
    .groups = "drop"
  )


ggplot(exp_summary, aes(x = log2(mean_expression + 1), y = mean_estimate)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", color = "black") +
  
  labs(
    x = "Mean gene expression (log2 counts)",
    y = "Mean intron effect size",
    title = "Relationship between gene expression and splicing efficiency"
  ) +
  
  theme_minimal()

cor(exp_summary$mean_expression,
    exp_summary$mean_estimate,
    method = "spearman")



