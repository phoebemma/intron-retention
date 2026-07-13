
# Load the file with the model and data
source("R/Load_data_for_visualisation.R")


# Plot global SE trajectory from the secondary data

plot_global_trajectory <- zi_predictions %>%
  mutate(
    SE       = 1 - estimate,
    CI_low   = 1 - conf.high,
    CI_high  = 1 - conf.low,
    real_age = scaled_age * (age_max - age_min) + age_min
  ) %>%
  group_by(real_age, time) %>%
  summarise(
    mean_SE = mean(SE),
    CI_low  = mean(CI_low),
    CI_high = mean(CI_high),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = real_age, y = mean_SE, colour = time, fill = time)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_smooth(method = "lm", se = FALSE,
              linetype = "dashed", linewidth = 0.6) +
  scale_colour_manual(values = c("PreExc" = colors[5],
                                 "PostExc" = colors[1])) +
  scale_fill_manual(values   = c("PreExc" = colors[5],
                                 "PostExc" = colors[1])) +
  labs(
    title    = "Global splicing efficiency trajectory across age",
    subtitle = "Solid = spline fit, dashed = linear reference",
    x        = "Age (years)",
    y        = "Mean splicing efficiency",
    colour   = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
     plot.subtitle = element_text(hjust = 0.5, size = 13),
    strip.text = element_text(face= "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12), 
    # legend.position = ,
    axis.text.y = element_text(size = 16, face= "bold"),
    axis.text.x = element_text(size = 16, face= "bold"),
    axis.title = element_text(size = 16, face= "bold"),
    panel.grid = element_blank(), 
  )

 

unique(relief_contrasts$hypothesis)


# First the effect of aging before exercise intervention begins
baseline_age_effect <- relief_contrasts %>%
  filter(hypothesis == "Age effect at baseline") #%>%
#filter(sig)

# Top labels per facet
top_introns <- baseline_age_effect %>%
  filter(sig) %>%
  slice_max(abs(estimate), n = 9, with_ties = FALSE)




# Count summary per facet
facet_summary <- baseline_age_effect %>%
  summarise(
    n_total  = n(),
    n_sig    = sum(sig),
    perc_sig = round(100 * n_sig / n_total, 1),
    x        = Inf,
    y        = Inf,
    label    = paste0("ds: ", n_sig, "/", n_total,
                      " (", perc_sig, "%)")
  )

baseline_age_volcano <- ggplot(baseline_age_effect,
                               aes(estimate, neg_log10_fdr, colour = effect)) +
  geom_point(alpha = 0.6, size = 3) +
  geom_text_repel(
    data        = top_introns,
    aes(label   = gene_intron),
    size        = 4,
    max.overlaps = Inf
  ) +
  geom_text(
    data         = facet_summary,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes  = FALSE,
    hjust        = -0.5,
    vjust        = 1.0,
    size         = 5,
    fontface     = "bold"
  ) +
  coord_cartesian(ylim = c(0, 8))+
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0,
             linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = effect_colors) +
  labs(
    title    = "Aging-related effects on intron splicing efficiency at baseline",
    x        = "Effect size",
    y        = expression(-log[10](FDR)),
    colour   = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    # plot.subtitle = element_text(hjust = 0.5, size = 13),
    strip.text = element_text(face= "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12), 
    # legend.position = ,
    axis.text.y = element_text(size = 16, face= "bold"),
    axis.text.x = element_text(size = 16, face= "bold"),
    axis.title = element_text(size = 16, face= "bold"),
    panel.grid = element_blank(), 
  )



# Becuase not all top ds introns in the primary model are present in the secondary
# We plot for every ds intron
ds_intron_aging <- baseline_age_effect %>%
  filter(sig) %>%
  arrange(desc(abs(estimate))) %>%
  slice_head(n = 22) #%>%
#pull(target)
#

x <- zi_predictions %>%
  filter(target %in% ds_intron_aging$target) %>% # extract top 9
  mutate(
    SE       = 1 - estimate,
    CI_low   = 1 - conf.high,
    CI_high  = 1 - conf.low,
    real_age = scaled_age * (age_max - age_min) + age_min
  ) %>%
  left_join(
    (zi_age_slopes_fdr  %>%
       dplyr::select(target, gene_intron) %>% distinct()),
    by = "target"
  ) %>%
  filter(time == "PreExc")

x %>%
  ggplot(aes(x = real_age, y = SE, colour =  gene_intron)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ gene_intron, scales = "free_y", ncol = 3) +
  labs(
    title    = "Age-related trajectory of some ds introns relected in secondary analysis",
    x        = "Age (years)",
    y        = "Splicing efficiency",
    colour   = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle   = element_text(hjust = 0.5),
    strip.text      = element_text(face = "bold", size = 14)
  )



# Unique introns
introns <- unique(x$gene_intron)

# Okabe-Ito palette
okabe_ito <- c(
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#F0E442", # yellow
  "#0072B2", # blue
  "#D55E00", # vermillion
  "#CC79A7", # reddish purple
  "#999999", # grey
  "#000000"  # black
)

# Match colours to introns
cols <- setNames(okabe_ito[seq_along(introns)], introns)

# Generate one plot per intron
age_trajectory_plots <- lapply(introns, function(intron) {
  
  pdat <- x %>%
    filter(gene_intron == intron)
  
  ggplot(pdat, aes(x = real_age, y = SE)) +
    geom_ribbon(
      aes(ymin = CI_low, ymax = CI_high),
      fill = cols[intron],
      alpha = 0.20
    ) +
    geom_line(
      colour = cols[intron],
      linewidth = 1
    ) +
    labs(
      title = intron,
      x = "Age (years)",
      y = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 12
      ),
      axis.title = element_text(size = 14),
      strip.text = element_text(face = "bold")
    )
})

# Assemble using patchwork
combined_plot <- wrap_plots(
  age_trajectory_plots,
  ncol = 3
) +
  plot_annotation(
    title = "Age-related trajectories of some differentially spliced introns",
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 14
      )
    )
  )

combined_plot

# ggsave("SVG_files/Age_related_trajectories.svg", plot = combined_plot, width = 16, height = 15)

# ggsave("SVG_files/Age_related_volcano_plot.svg", plot = baseline_age_volcano, width = 16, height = 15)





# combine volcano + insets

RT_final_volcano <- ggdraw() +
  draw_plot(baseline_age_volcano) +
  draw_plot(
    combined_plot,
    x = 0.55,   # adjust horizontally
    y = 0.55,  # adjust vertically
    width = 0.38,
    height = 0.4
  )

RT_final_volcano


df <- baseline_age_effect %>%
  filter(sig)

# Functional annotation of the genes affected
ego_RT <- enrichGO(gene =  df$external_gene_name,
                   keyType = "SYMBOL",
                   universe = gene_exp_df$gene_name,
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = T)


## Output results from GO analysis to a table
cluster_RT <- data.frame(ego_RT)

# go_RT <- dotplot(ego_RT,
#                  showCategory= 5,
#                  font.size = 12, title = "Enriched biological processes in genes containing introns with age-associated SE") +
#   theme(axis.text = element_text(size = 12, face = "bold"),
#         plot.title = element_text(hjust = 0.5, face = "bold") )
# 
# print(go_RT)
# 




# 
# Figure_3 <-  RT_final_volcano / go_RT 
# Figure_3 +
#   plot_annotation(tag_levels = "A") +
#   plot_layout(guides = "collect" , heights = c(1.2, 1)) &
#   theme(
#     plot.tag = element_text(size = 14, face = "bold"),
#     plot.tag.position = c(0.02, 0.98)
#   )
# 
# ggsave( "Figures/Trainome_fig_4.png",Figure_3,  width = 17, height = 12, dpi = 300)



# training effects
Training_effects <- relief_contrasts %>%
  filter(hypothesis == "Exercise effect among young participants" | hypothesis == "Exercise effect among older participants")


top_labels_training <- Training_effects %>%
  filter(sig) %>%
  group_by(hypothesis) %>%
  slice_max(abs(estimate), n = 9, with_ties = FALSE) %>%
  ungroup()




# Count summary per facet
facet_summary_training <- Training_effects %>%
  group_by(hypothesis) %>%
  summarise(
    n_total  = n(),
    n_sig    = sum(sig),
    perc_sig = round(100 * n_sig / n_total, 1),
    x        = Inf,
    y        = Inf,
    label    = paste0("ds: ", n_sig, "/", n_total,
                      " (", perc_sig, "%)"),
    .groups  = "drop"
  )

# Training_effects_plot <- 

ggplot(Training_effects,
       aes(estimate, neg_log10_fdr, colour = effect)) +
  geom_point(alpha = 0.6, size = 3) +
  geom_text_repel(
    data        = top_labels_training,
    aes(label   = gene_intron),
    size        = 4,
    max.overlaps = Inf
  ) +
  geom_text(
    data         = facet_summary_training,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes  = FALSE,
    hjust        = -0.5,
    vjust        = 1.0,
    size         = 5,
    fontface     = "bold"
  ) +
  coord_cartesian(ylim = c(0, 8))+
  facet_wrap(~hypothesis) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0,
             linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = effect_colors) +
  labs(
    title    = "Training-related effects on intron splicing efficiency",
    x        = "Effect size",
    y        = expression(-log[10](FDR)),
    colour   = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    # plot.subtitle = element_text(hjust = 0.5, size = 13),
    strip.text = element_text(face= "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12), 
    # legend.position = ,
    axis.text.y = element_text(size = 16, face= "bold"),
    axis.text.x = element_text(size = 16, face= "bold"),
    axis.title = element_text(size = 16, face= "bold"),
    panel.grid = element_blank(), 
  )



train_sig <- Training_effects %>%
  filter(sig)

ego_training <- enrichGO(gene =  train_sig$external_gene_name,
                         keyType = "SYMBOL",
                         universe = gene_exp_df$gene_name,
                         OrgDb = org.Hs.eg.db, 
                         ont = "BP", 
                         pAdjustMethod = "BH", 
                         qvalueCutoff = 0.05, 
                         readable = T)


## Output results from GO analysis to a table
cluster_train <- data.frame(ego_training)

