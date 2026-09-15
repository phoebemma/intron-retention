# Load the libraries and data
source("R/Load_data_for_visualisation.R") 


# Explore the Relief model

# The marginal effects predictions were built into thesummary function and thus
# Extracting the data will follow the traditional seqwrap method
Relief_results <- seqwrap_summarise(Relief_model)

Relief_results <- Relief_results$summaries




contrasts_of_interest <- c("train_young", "train_old", 
                           "age_effect_pre", "age_effect_post", 
                           "train_old_young")

relief_contrasts <- Relief_results %>%
  filter(hypothesis %in% contrasts_of_interest) %>%
  group_by( hypothesis, component) %>%
  mutate(adj.p = p.adjust(p.value, method = "fdr"),
         component = recode(component, 
                                   "response"    = "Overall retention",
                                   "zprob"       = "Probability of Perfect splicing",
                                   "conditional" = "Degree of retention"),
         hypothesis = recode(hypothesis, "train_young"     = "Exercise effect among young participants",
                                        "train_old"       = "Exercise effect among older participants",
                                        "age_effect_pre"  = "Age effect at baseline",
                                        "age_effect_post" = "Age effect postexercise",
                                        "train_old_young" = "Interaction effect of training and aging")) %>%
  ungroup() %>%
  mutate(
    sig           = adj.p <= 0.05,
    sig_ci        = conf.low > 0 | conf.high < 0,
    neg_log10_fdr = -log10(adj.p),
    rank_score    = -log10(adj.p) * abs(estimate),
    effect        = case_when(
      conf.high < 0 & sig ~ "Improved SE",
      conf.low  > 0 & sig ~ "Reduced SE",
      TRUE             ~ "No effect"
    )
  ) %>%
  annotate_introns(gene_annotation, intron_length, flip = TRUE) #%>%
 # filter(component == "Overall retention")




# extract top 9 introns per component and hypothesis

# Top labels per facet
top_labels <- relief_contrasts %>%
  filter(sig) %>%
  group_by(hypothesis, component) %>%
  slice_max(abs(estimate), n = 9, with_ties = FALSE) %>%
  ungroup()




# Count summary per facet
facet_summary <- relief_contrasts %>%
  group_by(hypothesis, component) %>%
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

ggplot(relief_contrasts,
       aes(estimate, neg_log10_fdr, colour = effect)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text_repel(
    data        = top_labels,
    aes(label   = gene_intron),
    size        = 3,
    max.overlaps = Inf
  ) +
  geom_text(
    data         = facet_summary,
    aes(x = x, y = y, label = label),
    inherit.aes  = FALSE,
    hjust        = 1.1,
    vjust        = 1.5,
    size         = 3,
    fontface     = "bold"
  ) +
  coord_cartesian(ylim = c(0, 10))+
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0,
             linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = effect_colors) +
  facet_grid(component ~ hypothesis, scales = "free") +
  labs(
    title    = "ReLiEf: exercise and age effects on splicing efficiency",
   # subtitle = "Rows = model component, Columns = contrast",
    x        = "Effect size (SE scale)",
    y        = expression(-log[10](FDR)),
    colour   = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 13),
    strip.text = element_text(face= "bold", size = 13),
    legend.title = element_blank(),
    legend.text = element_text(size = 12, face = "bold"), 
    # legend.position = "none",
    axis.text.y = element_text(size = 12, face= "bold"),
    axis.text.x = element_text(size = 12, face= "bold"),
    axis.title = element_text(size = 12, face= "bold")
  )

# subset and plot for each component

overall_retention <- relief_contrasts %>%
  filter(component == "Overall retention")

top_labels_overall <- overall_retention  %>%
  filter(sig) %>%
  group_by( hypothesis) %>%
  slice_max(abs(estimate), n = 9, with_ties = FALSE) %>%
  ungroup()


facet_summary_overall <- overall_retention %>%
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

ggplot(overall_retention,
       aes(estimate, neg_log10_fdr, colour = effect)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text_repel(
    data        = top_labels_overall,
    aes(label   = gene_intron),
    size        = 5,
    max.overlaps = Inf
  ) +
  geom_text(
    data         = facet_summary_overall,
    aes(x = x, y = y, label = label),
    inherit.aes  = FALSE,
    hjust        = 1.1,
    vjust        = 1.5,
    size         = 3,
    fontface     = "bold"
  ) +
  coord_cartesian(ylim = c(0, 8))+
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0,
             linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = effect_colors) +
  facet_grid(~hypothesis, scales = "free") +
  labs(
    title    = "Exercise and age effects on overall intron retention",
   # subtitle = "Rows = model component, Columns = contrast",
    x        = "Effect size (SE scale)",
    y        = expression(-log[10](FDR)),
    colour   = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 13),
    strip.text = element_text(face= "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 16, face = "bold"), 
    # legend.position = "none",
    axis.text.y = element_text(size = 16, face= "bold"),
    axis.text.x = element_text(size = 16, face= "bold"),
    axis.title = element_text(size = 16, face= "bold")
  )








# Compare this to the zi_conditional in the full model

zi_conditional


# Compare aging effect at baseline and training effect among the aged

# Which introns were differentially spliced among the aged at baseline
aging_baseline <- relief_contrasts %>%
  filter(hypothesis ==  "Age effect at baseline") %>%
  filter(component != "Overall retention")



sig_baseline <- aging_baseline %>%
  filter(sig)
length(unique(sig_baseline$gene_label))

table(sig_baseline$gene_label)

# how many genes contained more than 1 intron

repeated_genes <- sig_baseline %>%
  count(gene_label) %>%
  filter(n > 1) %>%
  pull(gene_label)



repeated_gene_rows <- sig_baseline %>%
  filter(gene_label %in% repeated_genes)

# Filyter the top nine for annotation 

# 
# top_9_labels <- aging_baseline %>%
#   filter(sig) %>%
#   group_by( component) %>%
#   slice_max(abs(estimate), n = 9) %>%
#   ungroup()


top_9_labels <- sig_baseline  %>%
  group_by(component) %>%
  mutate(rank_score = abs(estimate) * -log10(adj.p)) %>% 
  arrange(desc(rank_score), adj.p, desc(abs(estimate))) %>%
  slice_max(rank_score, n = 9, with_ties = FALSE) %>%
  ungroup()

# Count summary per facet
facet_sum <- aging_baseline %>%
  group_by( component) %>%
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

volcano_aging_baseline <- ggplot(aging_baseline,
       aes(estimate, neg_log10_fdr, colour = effect)) +
  geom_point(alpha = 0.6, size = 2.5) +
  geom_text_repel(
    data        = top_9_labels,
    aes(label   = gene_intron),
    size        = 4,
    max.overlaps = Inf
  ) +
  geom_text(
    data         = facet_sum,
    aes(x = x, y = y, label = label),
    inherit.aes  = FALSE,
    hjust        = 2.9,
    vjust        = 2.9,
    size         = 4,
    fontface     = "bold"
  ) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0,
             linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = effect_colors) +
  scale_y_continuous(limits = c(0, 16)) +  # adjust 20 to a sensible cap
  coord_cartesian(clip = "off")+
  facet_wrap(~ component, scales = "free") +
  labs(
    title    = "Introns affected by aging at baseline",
      subtitle = "Top 9 differentially spliced introns annotated",
    x        = "Effect size (SE scale)",
    y        = expression(-log[10](FDR)),
    colour   = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    strip.text = element_text(face= "bold", size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 12, face = "bold"), 
    # legend.position = "none",
    axis.text.y = element_text(size = 12, face= "bold"),
    axis.text.x = element_text(size = 12, face= "bold"),
    axis.title = element_text(size = 12, face= "bold")
  )

#ggsave("Figures/volc_baseline.png",  width = 26, height = 10)



ds_intron_aging <- sig_baseline %>%
  #arrange(desc(abs(estimate))) %>%
  mutate(rank_score = abs(estimate) * -log10(adj.p)) %>% 
  arrange(desc(rank_score), adj.p, desc(abs(estimate))) %>%
  slice_head(n = 59) 



aged_intersect <- zi_predictions %>%
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

length(unique(aged_intersect$gene_intron))

traj_aged <- aged_intersect %>%
  ggplot(aes(x = real_age, y = SE, colour =  gene_intron)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ gene_intron, scales = "free_y", ncol = 3) +
  labs(
    title    = "Age-related trajectory of the top 12 ds introns reflected in secondary analysis",
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

#ggsave("Figures/trajectories_of_12.png",  width = 26, height = 15)




volcano_aging_baseline + traj_aged +  plot_annotation(tag_levels = 'A')

ggsave("Figures/Figure_2.png",  width = 26, height = 15)





## EXPLORING THE EFFECT OF EXERCISE TRAINING

# explore effects of exercise among the old and young
Training_old <- relief_contrasts %>%
  filter(hypothesis == "Exercise effect among older participants"| hypothesis == "Exercise effect among young participants") %>%
  filter(component != "Overall retention")



# How many introns were significantly ds
sig_train_old <- Training_old %>%
  filter(hypothesis == "Exercise effect among older participants") %>%
  filter(sig)
length(unique(sig_train_old$gene_label))

table(sig_train_old$effect)


sig_train_young <- Training_old %>%
  filter(hypothesis == "Exercise effect among young participants") %>%
  filter(sig)


table(sig_train_young$effect)


duplicated_rows <- sig_train_old %>%
  group_by(gene_label) %>%
  filter(n() > 1) %>%
  ungroup()
table(duplicated_rows$gene_label)
length(unique(duplicated_rows$gene_label))


common_gene_introns <- intersect(sig_train_old$gene_intron, sig_train_young$gene_intron)
# Full rows from the "old" data frame, restricted to overlapping gene_introns
overlap_old <- sig_train_old %>%
  filter(gene_intron %in% common_gene_introns)

# Same for "young"
overlap_young <- sig_train_young %>%
  filter(gene_intron %in% common_gene_introns)



overlap_comparison <- sig_train_old %>%
  inner_join(sig_train_young, by = "gene_intron", suffix = c("_old", "_young")) %>%
  dplyr::select(gene_intron, estimate_old, estimate_young, adj.p_old, adj.p_young) %>%
  mutate(
    estimate_diff = estimate_old - estimate_young,
    same_direction = sign(estimate_old) == sign(estimate_young)
  )
table(overlap_comparison$same_direction)

top_9_labels_train <- Training_old %>%
  filter(sig) %>%
  group_by(hypothesis) %>%
  mutate(rank_score = abs(estimate) * -log10(adj.p)) %>% 
  arrange(desc(rank_score), adj.p, desc(abs(estimate))) %>%
  slice_max(rank_score, n = 9, with_ties = FALSE) %>%
  ungroup()


# Count summary per facet
facet_sum <- Training_old   %>%
  group_by( hypothesis) %>%
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

ggplot(Training_old ,
       aes(estimate, neg_log10_fdr, colour = effect)) +
  geom_point(alpha = 0.6, size = 2.5) +
  geom_text_repel(
    data        = top_9_labels_train,
    aes(label   = gene_intron),
    size        = 5,
    max.overlaps = Inf
  ) +
  geom_text(
    data         = facet_sum,
    aes(x = x, y = y, label = label),
    inherit.aes  = FALSE,
    hjust        = 2.9,
    vjust        = 2.9,
    size         = 4,
    fontface     = "bold"
  ) +
  coord_cartesian(ylim = c(0, 10))+
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0,
             linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = effect_colors) +
  facet_wrap(~ hypothesis, scales = "free") +
  labs(
    title    = "Effect of exercise in aged participants",
    #  subtitle = "Top 9 differentially spliced introns annotated",
    x        = "Effect size (SE scale)",
    y        = expression(-log[10](FDR)),
    colour   = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 13),
    strip.text = element_text(face= "bold", size = 13),
    legend.title = element_blank(),
    legend.text = element_text(size = 13, face = "bold"), 
    # legend.position = "none",
    axis.text.y = element_text(size = 13, face= "bold"),
    axis.text.x = element_text(size = 13, face= "bold"),
    axis.title = element_text(size = 13, face= "bold")
  )







# In the Relief data, exercise did not alter the probability of perfect splicing

# Take a look at the pooled model

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
  theme_minimal(base_size = 12) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle   = element_text(hjust = 0.5),
    #legend.position = "top"
  )

print(plot_global_trajectory)

#  Conditional trajectory (partially retained introns only)
plot_conditional <- zi_conditional %>%
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
  scale_colour_manual(values = c("PreExc" = colors[5],
                                 "PostExc" = colors[1])) +
  scale_fill_manual(values   = c("PreExc" = colors[5],
                                 "PostExc" = colors[1])) +
  labs(
    title    = "Conditional SE trajectory across age",
    subtitle = "Among introns that are not perfectly spliced",
    x        = "Age (years)",
    y        = "Mean splicing efficiency",
    colour   = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle   = element_text(hjust = 0.5),
   # legend.position = "top"
  )

print(plot_conditional)

#  Volcano age slopes per timepoint 
plot_volcano_age <- zi_age_slopes_fdr %>%
  ggplot(aes(x = estimate, y = neg_log10_fdr, colour = effect)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0,
             linetype = "dashed", colour = "grey50") +
  geom_text_repel(
    data = . %>% filter(sig) %>%
      group_by(time) %>%
      slice_max(abs(estimate), n = 9),
    aes(label = gene_intron),
    size = 3, max.overlaps = Inf
  ) +
  scale_colour_manual(values = effect_colors) +
  facet_wrap(~ time, ncol = 2) +
  labs(
    title    = "Age effect on splicing efficiency",
    subtitle = "After accounting for study-specific exercise responses",
    x        = "Age slope (effect on SE)",
    y        = expression(-log[10](FDR)),
    colour   = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle   = element_text(hjust = 0.5),
    strip.text      = element_text(face = "bold"),
    legend.position = "top"
  )

print(plot_volcano_age)

# Age slopes FDR within each timepoint ---
zi_age_slopes_fdr <- zi_age_slopes %>%
  group_by(time) %>%
  mutate(adj.p = p.adjust(p.value, method = "fdr")) %>%
  ungroup() %>%
  mutate(
    sig           = adj.p <= 0.05,
    sig_ci        = conf.low > 0 | conf.high < 0,
    neg_log10_fdr = -log10(adj.p),
    rank_score    = -log10(adj.p) * abs(estimate),
    effect        = case_when(
      conf.high < 0 & sig ~ "Improved SE",
      conf.low  > 0 & sig ~ "Reduced SE",
      TRUE                 ~ "No effect"
    )
  ) %>%
  annotate_introns(gene_annotation, intron_length, flip = TRUE)

# Exercise effects  FDR within each age anchor 
zi_time_effects_fdr <- zi_time_effects %>%
  group_by(scaled_age) %>%
  mutate(adj.p = p.adjust(p.value, method = "fdr")) %>%
  ungroup() %>%
  mutate(
    real_age      = scaled_age * (age_max - age_min) + age_min,
    sig           = adj.p <= 0.05,
    sig_ci        = conf.low > 0 | conf.high < 0,
    neg_log10_fdr = -log10(adj.p),
    rank_score    = -log10(adj.p) * abs(estimate),
    effect        = case_when(
      conf.high < 0 & sig ~ "Improved SE",
      conf.low  > 0 & sig ~ "Reduced SE",
      TRUE                 ~ "No effect"
    )
  ) %>%
  annotate_introns(gene_annotation, intron_length, flip = TRUE)


# Exercise effect trajectory across age (continuous) 
plot_exercise_trajectory <- zi_time_effects_fdr %>%
  group_by(real_age) %>%
  summarise(
    mean_estimate = mean(estimate, na.rm = TRUE),
    se            = sd(estimate, na.rm = TRUE) / sqrt(n()),
    .groups       = "drop"
  ) %>%
  ggplot(aes(x = real_age, y = mean_estimate)) +
  geom_ribbon(aes(ymin = mean_estimate - se, ymax = mean_estimate + se),
              alpha = 0.15, fill = colors[5], colour = NA) +
  geom_line(linewidth = 0.9, colour = colors[5]) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  labs(
    title    = "Global exercise effect on SE across age",
    subtitle = "After removing study-specific exercise responses",
    x        = "Age (years)",
    y        = "Exercise effect on SE (PostExc - PreExc)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(plot_exercise_trajectory)




#  Exercise effect trajectory for significant introns only 
sig_exercise_targets <- zi_time_effects_fdr %>%
  filter(sig) %>%
  pull(target) %>%
  unique()

plot_exercise_traj_sig <- zi_time_effects_fdr %>%
  filter(target %in% sig_exercise_targets) %>%
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

print(plot_exercise_traj_sig)



# Trajectory plots for top aging-associated introns 
top_aging_targets <- zi_age_slopes_fdr %>%
  filter(sig, time == "PostExc") %>%
  slice_max(abs(estimate), n = 6) %>%
  pull(target)

plot_traj_top_aging <- zi_predictions %>%
  filter(target %in% top_aging_targets) %>%
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

print(plot_traj_top_aging)


# training effect in he young 

Training_young <- relief_contrasts %>%
  filter(hypothesis == "Exercise effect among young participants") %>%
  filter(component != "Overall retention")



# How many introns were significantly ds
sig_train_young <-Training_young %>%
  filter(sig)

table(sig_train_young$effect)


# HOW MANY IN YOUNG Participants


# Trajectory plots for top introns in ReLiEf 
  Relief_sig_age_baseline <- Training_young%>%
  filter(sig) %>%
 # slice_max(abs(estimate), n = 10) %>%
  pull(target)

plot_traj_relief <- zi_predictions %>%
  filter(target %in% Relief_sig_age_baseline) %>%
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
    title    = "Top exercise-associated introns ReLiEf Young vs Old",
    x        = "Age group",
    y        = "Splicing efficiency",
    colour   = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    strip.text      = element_text(face = "bold", size = 8),
    legend.position = "top"
  )

print(plot_traj_relief)




#Interaction effect of aging
