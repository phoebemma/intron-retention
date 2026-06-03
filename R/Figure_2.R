

full_model_sum <- seqwrap_summarise(full_model)

volcano_df <- full_model_sum$summaries %>%
  dplyr::select(-group) %>%
  filter(term != "(Intercept)", term != "sexmale") %>%
  drop_na() %>%
  group_by(term) %>%
  mutate(
    adj.p = p.adjust(p.value, method = "fdr"),
    term = recode(term,
                  "scaled_age" = "Aging",
                  "timePostExc" = "Resistance Training")
  ) %>%
  ungroup() %>%
  mutate(
    sig = adj.p <= 0.05,
    neg_log10_fdr = -log10(adj.p)
  )


top10_labels <- volcano_df %>%
  filter(sig) %>%
  group_by(term) %>%
  slice_max(abs(estimate), n = 10, with_ties = FALSE) %>%
  ungroup()


term_summary <- volcano_df %>%
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


ggplot(volcano_df, aes(estimate, neg_log10_fdr)) +
  geom_point(aes(colour = sig), alpha = 0.7, size = 2) +
  
  geom_text_repel(
    data = top10_labels,
    aes(label = target),
    size = 3,
    max.overlaps = Inf
  ) +
  
  geom_text(
    data = term_summary,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1.1,
    vjust = 1.2,
    size = 3.5,
    fontface = "bold"
  ) +
  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  
  scale_colour_manual(values = c("grey70", "red")) +
  
  facet_wrap(~term, scales = "free") +
  
  labs(
    title = "Model estimates versus statistical significance of genes",
    x = "Model estimate",
    y = expression(-log[10]("FDR value"))
  ) +
  
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
