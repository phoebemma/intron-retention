# Load the file with the model and data
source("R/Load_data_for_visualisation.R")

library(patchwork)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(cowplot)
library(ggpubr)



zi_model <- zi_model$summaries
zi_df <- zi_pred_resp %>%
  dplyr::select(target, estimate, p.value)

length(unique(zi_model$target))

unique(zi_age_slopes$time)
# From the aging slopes extract a ranking variable

# Loage the average slope estimate per intron
zi_age_df <- zi_age_slopes %>%
  mutate(adj.p = p.adjust(p.value, method = "fdr"),
         neg_log10_fdr = -log10(adj.p),
         effect = case_when(
           estimate > 0 & adj.p <= 0.05 ~ "Decreased Splicing Efficiency",
           estimate < 0 & adj.p <= 0.05 ~ "Improved Splicing Efficiency",
           TRUE ~ "No effect"),
         rank_score = -log10(adj.p) * abs(estimate),
         sig = adj.p <= 0.05) %>% 
  annotate_introns(gene_annotation , intron_length) 


x <- zi_age_df %>%
  filter(sig) %>%
  count(term, effect) %>%
  filter(n > 1)

table(zi_age_slopes$time)
# extract summary statistics for annotating plot

summary <- zi_age_df %>%
  summarise(
    n_total = n(),
    n_sig = sum(sig),
    
    n_dec = sum(effect == "Improved Splicing Efficiency" & sig),
    n_inc = sum(effect == "Decreased Splicing Efficiency" & sig),
    
    perc_sig = 100 * n_sig / n_total,
    perc_dec = 100 * n_dec / n_total,
    perc_inc = 100 * n_inc / n_total,
    
    x = max(estimate, na.rm = TRUE),
    y = max(neg_log10_fdr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      " Improved SE: ", round(perc_dec, 1), "%\n",
      "Decreased SE: ", round(perc_inc, 1), "%)"
    )
  )

# extract the top 9 introns based on the ranking parameter (in the annotate_introns" function)
top9_introns <- zi_age_df %>%
  group_by(time) %>%
slice_head(n = 9) %>%
  ungroup()



volc_aging <- zi_age_df %>%
  ggplot(aes(x = estimate, y = neg_log10_fdr, colour = effect)) +
  geom_point(alpha = 0.9, size = 2) +
  geom_text_repel(
    data = top9_introns,
    aes(label = gene_intron),
    size = 6,
    max.overlaps = Inf
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+ 
  scale_color_manual(values = effect_colors) +
  geom_text(
    data = summary,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.1,
    vjust = 1.0,
    size = 5,
    fontface = "italic"
  ) +
  coord_cartesian(xlim = c(-0.02, 0.02))  + # adjust limits to where bulk of points are
  facet_wrap( ~ time) +
  
  labs(
    title = "Differentially Spliced(ds) Introns Due To Aging",
    x = "Effect size",
    y = expression(-log[10]("FDR value"), clip = "off")
  ) +
  
  # theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 13),
    legend.title = element_blank(),
    legend.text = element_text(size = 12, face = "bold"), 
    
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"), 
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )


# ggsave("SVG_files/Volcano_plot_aging.svg", plot = volc_aging, width = 14, height = 12)











top9_genes <- unique(top9_introns$gene_label)

age_introns_df <- all_splice_df %>%
  dplyr::filter(transcript_ID %in% top9_introns$target) %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID) ) %>%
  inner_join(top9_introns, by = c("transcript_ID" = "target")) %>%
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
      plot.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 14),
      axis.title = element_blank(),
      panel.grid = element_blank(), # removes grid lines 
      
      panel.border = element_rect( # add a subtle border
        colour = "grey80",   
        fill = NA,
        linewidth = 0.5
      )
    )
})


# stack them into a panel


mini_panel <- plot_grid(
  plotlist = intron_plots,
  ncol = 3
)



# ggsave("SVG_files/Top_6_aging_introns.svg", plot = mini_panel, width = 14, height = 12)

# Get the original age range
age_min <- min(metadata$age, na.rm = TRUE)
age_max <- max(metadata$age, na.rm = TRUE)

# Convert scaled_age back to real age
zi_pred_resp <- zi_pred_resp %>%
  mutate(real_age = scaled_age * (age_max - age_min) + age_min)





 zi_pred_resp <- zi_pred_resp %>%
  mutate(SE = 1 - estimate) %>%
group_by(time, scaled_age) %>%
  mutate( adj.p = p.adjust(p.value, method = "fdr")) %>%
  filter(adj.p <= 0.05) %>%
  ungroup() %>%
  group_by(real_age, time) %>%
  summarise(mean_SE = mean(SE, na.rm = TRUE),
            se      = sd(SE, na.rm = TRUE) / sqrt(n()),
            .groups = "drop") %>%
   ungroup()
 
 zi_pred_resp%>%
  ggplot(aes(x = real_age, y = mean_SE, colour = time, fill = time)) +
  geom_ribbon(aes(ymin = mean_SE - se, ymax = mean_SE + se), 
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_smooth(method = "lm", se = FALSE, 
              linetype = "dashed", linewidth = 0.6) +
  scale_colour_manual(values = c("PreExc" = colors[5], "PostExc" = colors[1])) +
  scale_fill_manual(values   = c("PreExc" = colors[5], "PostExc" = colors[1])) +
  labs(
    title    = "Global splicing efficiency trajectory across age",
    subtitle = "Solid = spline fit, dashed = linear fit",
    x        = "Scaled age",
    y        = "Mean splicing efficiency",
    colour   = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


# To classify introns by trajectory shape across four age windows
# Use the zpredictions to extract the direction of change at different age windows:

trajectory_class <- zi_pred_resp %>%
 # mutate(SE = 1 - estimate) %>%
  filter(time == "PreExc") %>%
 # group_by(target) %>%
  summarise(
    young       = mean(SE[scaled_age >= 0    & scaled_age < 0.25]),
    young_mid   = mean(SE[scaled_age >= 0.25 & scaled_age < 0.50]),
    mid_old     = mean(SE[scaled_age >= 0.50 & scaled_age < 0.75]),
    old         = mean(SE[scaled_age >= 0.75 & scaled_age <= 1.0]),
    .groups = "drop"
  ) %>%
  mutate(
    # slopes between consecutive groups
    slope_1 = young_mid - young,      # young to young-middle
    slope_2 = mid_old   - young_mid,  # young-middle to middle-old
    slope_3 = old       - mid_old,    # middle-old to old
    trajectory = case_when(
      # rises then falls
      slope_1 > 0 & slope_2 > 0 & slope_3 < 0 ~ "Rise-Rise-Decline",
      slope_1 > 0 & slope_2 < 0 & slope_3 < 0 ~ "Rise-Decline-Decline",
      slope_1 < 0 & slope_2 < 0 & slope_3 > 0 ~ "Decline-Decline-Rise",
      slope_1 < 0 & slope_2 > 0 & slope_3 > 0 ~ "Decline-Rise-Rise",
      slope_1 > 0 & slope_2 < 0 & slope_3 > 0 ~ "Rise-Decline-Rise",
      slope_1 < 0 & slope_2 > 0 & slope_3 < 0 ~ "Decline-Rise-Decline",
      slope_1 > 0 & slope_2 > 0 & slope_3 > 0 ~ "Consistent Improvement",
      slope_1 < 0 & slope_2 < 0 & slope_3 < 0 ~ "Consistent Decline",
      TRUE                                      ~ "Flat"
    )
  )

# Count introns by trajectory class
trajectory_class %>%
  count(trajectory) %>%
  arrange(desc(n)) %>%
  ggplot(aes(x = reorder(trajectory, n), y = n, fill = trajectory)) +
  geom_col() +
  geom_text(aes(label = n), hjust = -0.2) +
  coord_flip() +
  labs(
    title = "Intron trajectory patterns across four age groups",
    x     = NULL,
    y     = "Number of introns"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )



# Visualise trajectory per class
age_labels <- data.frame(
  scaled_age = c(0.125, 0.375, 0.625, 0.875),
  label      = c("Young", "Young-Mid", "Mid-Old", "Old"),
  mean_SE    = Inf
)

zi_pred_resp %>%
  mutate(SE = 1 - estimate) %>%
  filter(time == "PreExc") %>%
  inner_join(trajectory_class, by = "target") %>%
  group_by(trajectory, scaled_age) %>%
  summarise(
    mean_SE = mean(SE, na.rm = TRUE),
    se      = sd(SE, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = scaled_age, y = mean_SE, colour = trajectory, fill = trajectory)) +
  geom_ribbon(aes(ymin = mean_SE - se, ymax = mean_SE + se),
              alpha = 0.1, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = c(0.25, 0.50, 0.75),
             linetype = "dashed", colour = "grey70") +
  geom_text(data = age_labels, 
            aes(x = scaled_age, y = mean_SE, label = label),
            vjust = 1.5, size = 3.5, fontface = "italic",
            colour = "grey40", inherit.aes = FALSE) +
  facet_wrap(~ trajectory, scales = "free") +
  labs(
    title  = "Mean SE trajectory by intron class across four age groups",
    x      = "Scaled age",
    y      = "Mean splicing efficiency",
    colour = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    strip.text      = element_text(face = "bold")
  )



# Enrichment per trajectory class
trajectory_go <- trajectory_class %>%
  inner_join(
    zi_age_df %>% 
      dplyr::select(target, external_gene_name) %>% distinct(),
    by = "target"
  ) %>%
  group_by(trajectory) %>%
  summarise(genes = list(unique(external_gene_name)), .groups = "drop")

# Run GO for each trajectory class
go_results <- trajectory_go %>%
  mutate(
    ego = map(genes, ~ tryCatch(
      enrichGO(
        gene          = .x,
        keyType       = "SYMBOL",
        universe      = gene_exp_df$gene_name,
        OrgDb         = org.Hs.eg.db,
        ont           = "BP",
        pAdjustMethod = "BH",
        qvalueCutoff  = 0.05,
        readable      = TRUE
      ),
      error = function(e) NULL
    ))
  )

# Plot GO results per class
go_results %>%
  filter(!map_lgl(ego, is.null)) %>%
  mutate(plot = map2(ego, trajectory, ~ dotplot(.x, font.size = 8,
                                                title = .y))) %>%
  pull(plot) %>%
  wrap_plots(ncol = 2)







trajectory_class <- zi_pred_resp %>%
  mutate(SE = 1 - estimate) %>%
  filter(time == "PreExc") %>%  # focus on baseline age trajectory
  group_by(target) %>%
  summarise(
    # early age effect (young to middle)
    early_slope = mean(SE[scaled_age >= 0.25 & scaled_age <= 0.5]) -
      mean(SE[scaled_age >= 0 & scaled_age < 0.25]),
    # late age effect (middle to old)
    late_slope  = mean(SE[scaled_age > 0.5 & scaled_age <= 0.75]) -
      mean(SE[scaled_age > 0.5 & scaled_age <= 1.0]),
    .groups = "drop"
  ) %>%
  mutate(
    trajectory = case_when(
      early_slope > 0  & late_slope < 0  ~ "Rise then Decline",
      early_slope < 0  & late_slope > 0  ~ "Decline then Rise",
      early_slope > 0  & late_slope > 0  ~ "Consistent Improvement",
      early_slope < 0  & late_slope < 0  ~ "Consistent Decline",
      TRUE                               ~ "Flat"
    )
  )




trajectory_class %>%
  count(trajectory) %>%
  ggplot(aes(x = trajectory, y = n, fill = trajectory)) +
  geom_col() +
  geom_text(aes(label = n), vjust = -0.3) +
  theme_minimal() +
  labs(title = "Intron trajectory patterns across age",
       x = NULL, y = "Number of introns") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))



zi_pred_resp %>%
  mutate(SE = 1 - estimate) %>%
  filter(time == "PreExc") %>%
  inner_join(trajectory_class, by = "target") %>%
  group_by(trajectory, scaled_age) %>%
  summarise(mean_SE = mean(SE, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = scaled_age, y = mean_SE, colour = trajectory)) +
  geom_line(linewidth = 0.9) +
  labs(
    title  = "Mean SE trajectory by intron class",
    x      = "Scaled age",
    y      = "Mean splicing efficiency",
    colour = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))



# Example for "Rise then Decline" introns
trajectory_go <- trajectory_class %>%
  inner_join(
    zi_age_df %>% 
      dplyr::select(target, external_gene_name) %>% distinct(),
    by = "target"
  ) %>%
  group_by(trajectory) %>%
  summarise(genes = list(unique(external_gene_name)), .groups = "drop")

# Run GO for each trajectory class
go_results <- trajectory_go %>%
  mutate(
    ego = map(genes, ~ tryCatch(
      enrichGO(
        gene          = .x,
        keyType       = "SYMBOL",
        universe      = gene_exp_df$gene_name,
        OrgDb         = org.Hs.eg.db,
        ont           = "BP",
        pAdjustMethod = "BH",
        qvalueCutoff  = 0.05,
        readable      = TRUE
      ),
      error = function(e) NULL
    ))
  )

# Plot GO results per class
go_results %>%
  filter(!map_lgl(ego, is.null)) %>%
  mutate(plot = map2(ego, trajectory, ~ dotplot(.x, font.size = 8,
                                                title = .y))) %>%
  pull(plot) %>%
  wrap_plots(ncol = 2)



# The zprob model that looks at how perfect splicing changes with age



zprob <- zi_pred_zprob %>%
  group_by(time) %>%
  mutate(adj.p = p.adjust(p.value, method = "fdr"),
         neg_log10_fdr = -log10(adj.p),
         effect = case_when(
           estimate > 0 & adj.p <= 0.05 ~ "Decreased Splicing Efficiency",
           estimate < 0 & adj.p <= 0.05 ~ "Improved Splicing Efficiency",
           TRUE ~ "No effect"),
         rank_score = -log10(adj.p) * abs(estimate),
         sig = adj.p <= 0.05) %>% 
  annotate_introns(gene_annotation , intron_length) %>%
  ungroup()
  
  
zprob %>%
  ggplot(aes(x = estimate, y = neg_log10_fdr, colour = effect)) +
  geom_point(alpha = 0.9, size = 2) +
  # geom_text_repel(
  #   data = top9_introns,
  #   aes(label = gene_intron),
  #   size = 6,
  #   max.overlaps = Inf
  # ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+ 
  scale_color_manual(values = effect_colors) +
  geom_text(
    data = summary,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.1,
    vjust = 1.0,
    size = 5,
    fontface = "italic"
  )
  
  
  
  
zprob <- zprob %>%
  filter(adj.p <= 0.05)

length(unique(zprob$target))















  group_by(scaled_age, time) %>%
  summarise(mean_zprob = mean(estimate, na.rm = TRUE),
            se         = sd(estimate, na.rm = TRUE) / sqrt(n()),
            .groups    = "drop") %>%
  ggplot(aes(x = scaled_age, y = mean_zprob, colour = time, fill = time)) +
  geom_ribbon(aes(ymin = mean_zprob - se, ymax = mean_zprob + se),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = c("PreExc" = colors[5], "PostExc" = colors[1])) +
  scale_fill_manual(values   = c("PreExc" = colors[5], "PostExc" = colors[1])) +
  labs(
    title = "Probability of perfect splicing across age",
    x     = "Scaled age",
    y     = "P(SE = 1)",
    colour = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))



sig_zprob_introns <- zi_pred_zprob %>%
  group_by(target) %>%
  mutate(adj.p = p.adjust(p.value, method = "fdr")) %>%
  filter(adj.p <= 0.05) %>%
  pull(target) %>%
  unique()




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

