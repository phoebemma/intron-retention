# Load the file with the model and data
source("R/Load_data_for_visualisation.R")

library(patchwork)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(cowplot)
library(ggpubr)





# Load the binomial model output
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
  separate(target, into = c(NA, "intron_ID", NA), sep = "_", remove = F) %>%
  # create a gene and intron label using the gene name and intron number
  mutate(
    gene_label = ifelse(
      is.na(external_gene_name) | external_gene_name == "",
      ensembl_gene_id,
      external_gene_name
    ),
    gene_intron = paste(gene_label, intron_ID, sep = " : ")) # %>%
  # arrange(gene_label, estimate) %>%
  # mutate(gene_intron = factor(gene_intron, levels = unique(gene_intron))) %>%
  # dplyr::select(intron_ID, term, effect, adj.p, estimate, p.value, intron_length,
  #               sig, neg_log10_fdr, gene_intron, gene_label, transcript_biotype)





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



aging_table <- beta_binom_df %>%
  filter(term == "Aging") %>%
  mutate(rank_score = -log10(adj.p) * abs(estimate)) %>% # create a ranking variable
  arrange(desc(rank_score)) #%>%
#dplyr::select(gene_intron, effect, estimate, adj.p, gene_label, rank_score)

#  saveRDS(aging_table, "tables/Trainome_aging_table.rds")

summary <- aging_table %>%
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



 top6_introns <- aging_table %>% slice_head(n = 6)


top6_genes <- unique(top6_introns$gene_label)

age_introns_df <- all_splice_df %>%
  dplyr::filter(transcript_ID %in% top6_introns$target) %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID) ) %>%
  inner_join(top6_introns, by = c("transcript_ID" = "target")) %>%
  inner_join(metadata, by = "seq_sample_id")%>% 
  group_by(scaled_age, gene_intron, transcript_ID) %>%
  summarise(mean_SE = mean(SE, na.rm = TRUE), .groups = "drop")

# Create SE plots, one per intron




intron_plots <- lapply(unique(age_introns_df$gene_intron), function(intron) {
  ggplot(
    age_introns_df %>% filter(gene_intron == intron),
    aes(x = scaled_age, y = mean_SE)
  ) +
    geom_point(size = 1.5, alpha = 0.9, colour = "black") +
    geom_smooth(method = "lm", se = FALSE, colour = "red", linewidth = 0.6) +
    labs(title = intron, x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      axis.text = element_text(size = 10),
      axis.title = element_blank()
    )
})


# stack them into a panel


mini_panel <- plot_grid(
  plotlist = intron_plots,
  ncol = 3
)


volcano_plot <- ggplot((beta_binom_df %>%
                          filter(term == "Aging")), aes(estimate, neg_log10_fdr, colour = effect)) +
  geom_point(aes(colour = effect), alpha = 0.7, size = 2) +
  
  geom_text_repel(
    data = top6_introns,
    aes(label = gene_intron),
    size = 6,
    max.overlaps = Inf
  ) +
  
  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  
  scale_color_manual(values = c("Improved SE" =  colors[6], "Reduced SE" = colors[1]),
                     name = "Effect") +
  
  geom_text(
    data = summary,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.1,
    vjust = 1.0,
    size = 5,
    fontface = "italic"
  ) +
  
  labs(
    title = "Differentially spliced(ds) introns due to Aging",
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

final_volcano <- ggdraw() +
  draw_plot(volcano_plot) +
  draw_plot(
    mini_panel,
    x = 0.55,   # adjust horizontally
    y = 0.45,  # adjust vertically
    width = 0.38,
    height = 0.48
  )

final_volcano



# functional annotation of genes with introns whose SE are aging-associated

aging_df <- beta_binom_df %>%
  filter(term == "Aging") %>%
  dplyr::filter(adj.p <= 0.05)
# Functional annotation of the genes affected
ego_aging <- enrichGO(gene =  aging_df$external_gene_name,
                      keyType = "SYMBOL",
                      universe = gene_exp_df$gene_name,
                      OrgDb = org.Hs.eg.db, 
                      ont = "BP", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05, 
                      readable = T)


## Output results from GO analysis to a table
cluster_aging <- data.frame(ego_aging)

go_aging <- dotplot(ego_aging,
                    showCategory = 5,
                    font.size = 12, title = "Enriched biological processes in genes containing introns with aging-associated SE") +
  theme(axis.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold") )

go_aging 



# plot splicing efficiency by gene expression


# First load the gene expression dataset


# select the top 6 most age-affected introns
top6_introns <- aging_table %>% slice_head(n = 6)

top6_genes <- unique(top6_introns$gene_label)

splice_df_long <- all_splice_df %>%
  pivot_longer(
    cols = -transcript_ID,
    names_to = "seq_sample_id",
    values_to = "SE"
  )




age_exp_df <- gene_exp_df %>%
  filter(gene_name %in% beta_binom_df$external_gene_name) %>%
  pivot_longer(
    cols = -gene_name,
    names_to = "seq_sample_id",
    values_to = "gene_count"
  ) %>%
  inner_join(metadata, by = "seq_sample_id") %>%
  inner_join(beta_binom_df, by = c("gene_name" = "external_gene_name")) %>%
  inner_join(splice_df_long, by = c("target" = "transcript_ID", "seq_sample_id"))%>%
  group_by(gene_name, seq_sample_id) %>%
  summarise(
    gene_count = dplyr::first(gene_count),
    mean_SE = mean(SE, na.rm = TRUE)
  ) %>%
  ungroup()



plot_df <- age_exp_df %>%
  filter(gene_name %in% top6_genes)
plot_df$log_expr <- log2(plot_df$gene_count + 1)

cor_plot <- ggplot(plot_df, aes(x = mean_SE, y = log_expr, colour = gene_name) ) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "red") +
  facet_wrap(~ gene_name, scales = "free") +
  stat_cor(method = "spearman", size = 5) +
  theme_minimal() +
  labs(
    x = "Splicing efficiency (SE)",
    y = NULL,
    title = "SE vs expression in genes with most aging-associated introns"
  )+
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
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

# 
# age_exp_df <- age_exp_df %>%
#   mutate(highlight = ifelse(gene_name %in% top6_genes, gene_name, "Other"))
# 
# ggplot() +
#   geom_point(data = age_exp_df,
#              aes(x = mean_SE, y = gene_count),
#              color = "grey80", alpha = 0.4) +
#   
#   geom_point(data = subset(age_exp_df, gene_name %in% top6_genes),
#              aes(x = mean_SE, y = gene_count, color = gene_name),
#              size = 2) +
#   
#   theme_minimal() +
#   labs(
#     x = "Splicing efficiency (SE)",
#     y = "Normalized gene expression",
#     color = "Top genes"
#   )
# 
# 
# 
# 
# ggplot(age_exp_df, aes(x = SE, y = gene_count)) +
#   geom_point(aes(color = highlight), alpha = 0.6) +
#   scale_color_manual(values = c(rep("red", 6), "grey70")) +
#   theme_minimal() +
#   labs(
#     x = "Splicing efficiency (SE)",
#     y = "Normalized gene expression",
#     color = "Gene"
#   )




# 
# layout <- "
# AB
# AC
# "
# 
# 
# 
# 
# Figure_1 <- final_volcano + go_aging  + cor_plot+
#   plot_layout(design = layout, widths = c(1.2,1))
# 
# Figure_1 +
#   plot_annotation(tag_levels = "A", 
#                   
#                   tag_prefix = "",
#                   tag_suffix = ""
#   ) +
#   plot_layout(guides = "collect") &
#   theme(
#     plot.tag = element_text(size = 16, face = "bold"),
#     plot.tag.position = c(0.08, 0.98)
#   ) 



Figure_1 <-  final_volcano / (go_aging | cor_plot)
Figure_1 +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect" , widths = c(0.25, 1, 0.5)) &
  theme(
    plot.tag = element_text(size = 14, face = "bold"),
    plot.tag.position = c(0.02, 0.98)
  )

ggsave( "Figures/1Trainome_fig_3.png",Figure_1,  width = 17, height = 12, dpi = 300)

