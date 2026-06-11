
# Load the file with the model and data
source("R/Load_data_for_visualisation.R")

library(patchwork)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# Load the beta-binomial model output 

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
  separate(target, into = c(NA, "intron_ID", NA), sep = "_", remove = F) %>%
  # create a gene and intron label using the gene name and intron number
  mutate(
    gene_label = ifelse(
      is.na(external_gene_name) | external_gene_name == "",
      ensembl_gene_id,
      external_gene_name
    ),
    gene_intron = paste(gene_label, intron_ID, sep = " : ")
  ) 



RT_table <- beta_binom_df %>%
  filter(term == "Resistance Training") %>%
  mutate(rank_score = -log10(adj.p) * abs(estimate)) %>%
  arrange(desc(rank_score))



RT_summary <- RT_table %>%
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
      "Ds introns: ",
      n_sig, "/", n_total,
      " (", round(perc_sig, 1), "%)"
    )
  )


# select the top 6 most age-affected introns
top6_RT_introns <- RT_table %>% slice_head(n = 10)

top6_RT_genes <- unique(top6_RT_introns$gene_label)



RT_introns_df <- all_splice_df %>%
  dplyr::filter(transcript_ID %in% top6_RT_introns$target) %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID) ) %>%
  inner_join(top6_RT_introns, by = c("transcript_ID" = "target")) %>%
  inner_join(metadata, by = "seq_sample_id")%>% 
  group_by( time, gene_intron, transcript_ID) %>%
  summarise(mean_SE = mean(SE, na.rm = TRUE), .groups = "drop")

# Create SE plots, one per intron




RT_intron_plots <- lapply(unique(RT_introns_df$gene_intron), function(intron) {
  ggplot(
    RT_introns_df %>% filter(gene_intron == intron),
    aes(x = time, y = mean_SE, colour = time, e)
  ) +
    geom_point(size = 3, alpha = 0.9) +
    geom_smooth(method = "lm", se = FALSE, size = 0.6) +
    labs(title = intron, x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      axis.title = element_blank()
    )
})


# stack them into a panel


RT_mini_panel <- plot_grid(
  plotlist = RT_intron_plots,
  ncol = 3
)


RT_volcano_plot <- ggplot((beta_binom_df %>%
                             filter(term == "Resistance Training")), aes(estimate, neg_log10_fdr, colour = effect)) +
  geom_point(aes(colour = effect), alpha = 0.7, size = 2) +
  
  geom_text_repel(
    data = top6_RT_introns,
    aes(label = gene_intron),
    size = 4,
    max.overlaps = Inf
  ) +
  
  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text(
    data = RT_summary,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.1,
    vjust = 1.0,
    size = 5,
    fontface = "italic"
  )+
  
  scale_color_manual(values = c("Improved SE" =  colors[6], "Reduced SE" = colors[1]),
                     name = "Effect") +
  
  # facet_wrap(~term, scales = "free") +
  # coord_cartesian(xlim = c(-1, 1))+
  
  labs(
    title = "Differentially spliced introns due to Resistance Training",
    subtitle = "Beta-binomial model (splicing efficiency coded as 0,1)",
    x = "Effect size",
    y = expression(-log[10]("FDR value"), clip = "off")
  ) +
  
  # theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 13),
    legend.title = element_blank(),
    legend.text = element_text(size = 9, face = "bold"), 
    
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"), 
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )


# combine volcano + insets

RT_final_volcano <- ggdraw() +
  draw_plot(RT_volcano_plot) +
  draw_plot(
    RT_mini_panel,
    x = 0.55,   # adjust horizontally
    y = 0.55,  # adjust vertically
    width = 0.38,
    height = 0.4
  )

RT_final_volcano


RT_df <- RT_table %>%
  filter(adj.p <= 0.05)

# Functional annotation of the genes affected
ego_RT <- enrichGO(gene =  RT_df$external_gene_name,
                   keyType = "SYMBOL",
                   universe = gene_exp_df$gene_name,
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = T)


## Output results from GO analysis to a table
cluster_RT <- data.frame(ego_RT)

go_RT <- dotplot(ego_RT,
                 showCategory= 5,
                 font.size = 12, title = "Enriched biological processes in genes containing introns with RT-associated SE") +
  theme(axis.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold") )

print(go_RT)






Figure_3 <-  RT_final_volcano / go_RT 
Figure_3 +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect" , heights = c(1.2, 1)) &
  theme(
    plot.tag = element_text(size = 14, face = "bold"),
    plot.tag.position = c(0.02, 0.98)
  )

ggsave( "Figures/Trainome_fig_4.png",Figure_3,  width = 17, height = 12, dpi = 300)



