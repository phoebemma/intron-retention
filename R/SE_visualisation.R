relief_contrasts_df <- relief_contrasts %>%
  filter(component != "Overall retention")
unique(relief_contrasts_df$hypothesis)

baseline_age_effect <- relief_contrasts_df %>%
  filter(hypothesis == "Age effect at baseline") %>%
  filter(sig)

table(baseline_age_effect$effect)


# SUPPLEMENATARY FIGURE 1A
ggplot(baseline_age_effect , aes(x = estimate, y = reorder(gene_intron, estimate), color = effect)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  #facet_wrap(~ hypothesis, scales = "free") +
  scale_color_manual(values = c("Improved SE" = colors[6],
                                "Reduced SE" = colors[1]),
                     name = "Effect") +
  labs(
    x = "Effect size",
    y = NULL,
    title = "Introns with age-related splicing efficiency at baseline"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 13),
    strip.text = element_text(face= "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 16, face = "bold"), 
    # legend.position = "none",
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16, face= "bold"),
    axis.title = element_text(size = 16, face= "bold")
  )

# Functional annotation of the genes affected
ego_aging <- enrichGO(gene =  unique(baseline_age_effect$external_gene_name),
                   keyType = "SYMBOL",
                   universe = gene_exp_df$gene_name,
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = T)


## Output results from GO analysis to a table
cluster_RT <- data.frame(ego_aging)

go_aging <- dotplot(ego_aging,
                 
                 font.size = 8, title = "Enriched biological processes in genes containing introns with aging-associated SE") +
  theme(axis.text = element_text(size = 10), axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 10),
        plot.title = element_text(hjust = 0) )

print(go_aging)


# extract the top 9 of these introns based on ranking and visualise their trajectory

top_10_baseline_ints <- baseline_age_effect %>%
  group_by(hypothesis) %>%
  slice_max(abs(estimate), n = 9, with_ties = FALSE) %>%
  ungroup()


zi_time_effects_fdr %>%
  filter(target %in% top_10_baseline_ints$target) %>%
  group_by(real_age, effect) %>%
  summarise(
    mean_estimate = mean(estimate, na.rm = TRUE),
    se            = sd(estimate, na.rm = TRUE) / sqrt(n()),
    .groups       = "drop"
  ) %>%
  ggplot(aes(x = real_age, y = mean_estimate,
             colour = effect, fill = effect)) +
  geom_ribbon(aes(ymin = mean_estimate - se, ymax = mean_estimate + se),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = effect_colors) +
  scale_fill_manual(values   = effect_colors) +
  labs(
    title    = "Exercise effect on SE across age significant introns",
    subtitle = "Split by direction of effect",
    x        = "Age (years)",
    y        = "Exercise effect on SE (PostExc - PreExc)",
    colour   = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle   = element_text(hjust = 0.5)#,
    # legend.position = "top"
  )


zi_predictions %>%
  filter(target %in% baseline_age_effect$target) %>%
  mutate(
    SE       = 1 - estimate,
    CI_low   = 1 - conf.high,
    CI_high  = 1 - conf.low,
    real_age = scaled_age * (age_max - age_min) + age_min
  ) %>%
  left_join(
    zi_age_slopes_fdr %>%
      dplyr::select(target, gene_intron) %>% distinct(),
    by = "target"
  ) %>%
  ggplot(aes(x = real_age, y = SE, colour = time, fill = time)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ gene_intron, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = c("PreExc" = colors[5],
                                 "PostExc" = colors[1])) +
  scale_fill_manual(values   = c("PreExc" = colors[5],
                                 "PostExc" = colors[1])) +
  labs(
    title    = "Top aging-associated introns SE trajectory",
    subtitle = "Pooled model with study-specific exercise effects removed",
    x        = "Age (years)",
    y        = "Splicing efficiency",
    colour   = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle   = element_text(hjust = 0.5),
    strip.text      = element_text(face = "bold", size = 8),
    legend.position = "top"
  )


# age effect at post exercise


post_age_effect <- relief_contrasts_df %>%
  filter(hypothesis == "Age effect postexercise") %>%
  filter(sig)

table(post_age_effect$effect)
# post_age_effect$gene_intron %in% baseline_age_effect$gene_intron
# intersect(post_age_effect$gene_intron, baseline_age_effect$gene_intron)

ggplot(post_age_effect , aes(x = estimate, y = reorder(gene_intron, estimate), color = effect)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  #facet_wrap(~ hypothesis, scales = "free") +
  scale_color_manual(values = c("Improved SE" = colors[6],
                                "Reduced SE" = colors[1]),
                     name = "Effect") +
  labs(
    x = "Effect size",
    y = NULL,
    title = "Introns with age-related splicing efficiency at postexercise"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 13),
    strip.text = element_text(face= "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 16, face = "bold"), 
    # legend.position = "none",
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16, face= "bold"),
    axis.title = element_text(size = 16, face= "bold")
  )





## PLOT THE INTRONS AFFECTED BY AGING AT POST AND PRE_EXERCISE
comparison <- baseline_age_effect %>%
  dplyr::select(gene_intron, baseline_estimate = estimate) %>%
  inner_join(
    post_age_effect %>%
      dplyr::select(gene_intron, post_estimate = estimate),
    by = "gene_intron"
  )

ggplot(comparison,
       aes(x = baseline_estimate,
           y = post_estimate,
           label = gene_intron)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed",
              color = "grey50") +
  ggrepel::geom_text_repel() +
  labs(
    x = "Baseline effect size",
    y = "Postexercise effect size",
    title = "Comparison of age effects across timepoints"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 13),
    strip.text = element_text(face= "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 16, face = "bold"), 
    # legend.position = "none",
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16, face= "bold"),
    axis.title = element_text(size = 16, face= "bold")
  )



# Looking at introns affected by Training in the old and young

Training_effect_old <- relief_contrasts_df %>%
  filter(hypothesis == "Exercise effect among older participants") %>%
  filter(sig)


ggplot(Training_effect_old , aes(x = estimate, y = reorder(gene_intron, estimate), color = effect)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  #facet_wrap(~ hypothesis, scales = "free") +
  scale_color_manual(values = c("Improved SE" = colors[6],
                                "Reduced SE" = colors[1]),
                     name = "Effect") +
  labs(
    x = "Effect size",
    y = NULL,
    title = "Introns with training-related splicing efficiency among the aged"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 13),
    strip.text = element_text(face= "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 16, face = "bold"), 
    # legend.position = "none",
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 16, face= "bold"),
    axis.title = element_text(size = 16, face= "bold")
  )



# determine  any  significant gene ontology 
# Functional annotation of the genes affected
ego_RT_old <- enrichGO(gene =  unique(Training_effect_old$external_gene_name),
                      keyType = "SYMBOL",
                      universe = gene_exp_df$gene_name,
                      OrgDb = org.Hs.eg.db, 
                      ont = "BP", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05, 
                      readable = T)


## Output results from GO analysis to a table
cluster_RT_old <- data.frame(ego_RT_old)

 dotplot(ego_RT_old,
                    
                    font.size = 8, title = "Enriched biological processes in genes containing introns with aging-associated SE") +
  theme(axis.text = element_text(size = 10), axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 10),
        plot.title = element_text(hjust = 0) )






Training_effect_young <- relief_contrasts_df %>%
  filter( hypothesis == "Exercise effect among young participants") %>%
  filter(sig)
table(Training_effect_young$effect)

ggplot(Training_effect_young , aes(x = estimate, y = reorder(gene_intron, estimate), color = effect)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  #facet_wrap(~ hypothesis, scales = "free") +
  scale_color_manual(values = c("Improved SE" = colors[6],
                                "Reduced SE" = colors[1]),
                     name = "Effect") +
  labs(
    x = "Effect size",
    y = NULL,
    title = "Introns with training-related splicing efficiency among the young"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 13),
    strip.text = element_text(face= "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 16, face = "bold"), 
    # legend.position = "none",
    axis.text.y = element_text(size = 13),
    axis.text.x = element_text(size = 16, face= "bold"),
    axis.title = element_text(size = 16, face= "bold")
  )


# Functional annotation of the genes affected
ego_RT_young <- enrichGO(gene =  unique(Training_effect_young$external_gene_name),
                       keyType = "SYMBOL",
                       universe = gene_exp_df$gene_name,
                       OrgDb = org.Hs.eg.db, 
                       ont = "BP", 
                       pAdjustMethod = "BH", 
                       qvalueCutoff = 0.05, 
                       readable = T)


## Output results from GO analysis to a table
cluster_RT_young <- data.frame(ego_RT_young)

dotplot(ego_RT_old,
        
        font.size = 8, title = "Enriched biological processes in genes containing introns with aging-associated SE") +
  theme(axis.text = element_text(size = 10), axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 10),
        plot.title = element_text(hjust = 0) )

# filter and plot the trajectory of introns affected by training among the old

zi_time_effects_fdr %>%
  filter(gene_intron %in% Training_effect_old$gene_intron) %>%
  group_by(real_age, effect) %>%
  summarise(
    mean_estimate = mean(estimate, na.rm = TRUE),
    se            = sd(estimate, na.rm = TRUE) / sqrt(n()),
    .groups       = "drop"
  ) %>%
  ggplot(aes(x = real_age, y = mean_estimate,
             colour = effect, fill = effect)) +
  geom_ribbon(aes(ymin = mean_estimate - se, ymax = mean_estimate + se),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = effect_colors) +
  scale_fill_manual(values   = effect_colors) +
  labs(
    title    = "Exercise effect on SE across age — significant introns",
    subtitle = "Split by direction of effect",
    x        = "Age (years)",
    y        = "Exercise effect on SE (PostExc - PreExc)",
    colour   = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle   = element_text(hjust = 0.5)#,
    # legend.position = "top"
  )
  

# Trajectory plots for top training-associated introns 
# top_Training_targets <- Training_effect_old %>%
#   filter(sig) %>%
#   slice_max(abs(estimate), n = 20) %>%
#   pull(target)


zi_predictions %>%
  filter(target %in% Training_effect_young$target) %>%
  mutate(
    SE       = 1 - estimate,
    CI_low   = 1 - conf.high,
    CI_high  = 1 - conf.low,
    real_age = scaled_age * (age_max - age_min) + age_min
  ) %>%
  left_join(
    zi_age_slopes_fdr %>%
      dplyr::select(target, gene_intron) %>% distinct(),
    by = "target"
  ) %>%
  ggplot(aes(x = real_age, y = SE, colour = time, fill = time)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ gene_intron, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = c("PreExc" = colors[5],
                                 "PostExc" = colors[1])) +
  scale_fill_manual(values   = c("PreExc" = colors[5],
                                 "PostExc" = colors[1])) +
  labs(
    title    = "Top aging-associated introns — SE trajectory",
    x        = "Age (years)",
    y        = "Splicing efficiency",
    colour   = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle   = element_text(hjust = 0.5),
    strip.text      = element_text(face = "bold", size = 8),
    legend.position = "top"
  )



# Visualise introns that are associated with the interaction between aging and training

inte_effect <- relief_contrasts_df %>%
  filter(hypothesis == "Interaction effect of training and aging") %>%
  filter(sig)


# top_inte_targets <- inte_effect %>%
#   filter(sig) %>%
#   slice_max(abs(estimate), n = 6) %>%
#   pull(target)


zi_predictions %>%
  filter(target %in% inte_effect$target) %>%
  mutate(
    SE       = 1 - estimate,
    CI_low   = 1 - conf.high,
    CI_high  = 1 - conf.low,
    real_age = scaled_age * (age_max - age_min) + age_min
  ) %>%
  left_join(
    zi_age_slopes_fdr %>%
      dplyr::select(target, gene_intron) %>% distinct(),
    by = "target"
  ) %>%
  ggplot(aes(x = real_age, y = SE, colour = time, fill = time)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ gene_intron, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = c("PreExc" = colors[5],
                                 "PostExc" = colors[1])) +
  scale_fill_manual(values   = c("PreExc" = colors[5],
                                 "PostExc" = colors[1])) +
  labs(
    title    = "Introns associated with interactive effect of aging and exercise SE trajectory",
    x        = "Age (years)",
    y        = "Splicing efficiency",
    colour   = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle   = element_text(hjust = 0.5),
    strip.text      = element_text(face = "bold", size = 8),
    legend.position = "top"
  )


