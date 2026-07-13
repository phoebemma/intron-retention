# load the file with the data
source("R/Load_data_for_visualisation.R")

library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggrepel)


# This investigates interaction effect
zi_interaction_time <- readRDS("data/zi_main_time_effects.RDS")
# which introns were signififcant across the different ages
sig_any_age <- zi_interaction_time %>%
  group_by(scaled_age) %>%
  mutate(adj.p = p.adjust(p.value, method = "fdr")) %>%
  ungroup() %>%
  filter(adj.p <= 0.05) #%>%
  #pull(target) %>%
  #unique()



zi_interaction_time %>%
  group_by(scaled_age) %>%
  mutate(
    # CI does not cross zero = confident effect exists
    sig_ci = conf.low > 0 | conf.high < 0,
    # Direction of effect
    effect = case_when(
      conf.high < 0 ~ "Improved SE",   # entire CI below zero = confident improvement
      conf.low  > 0 ~ "Reduced SE",    # entire CI above zero = confident reduction
      TRUE          ~ "Uncertain"       # CI crosses zero = no confident direction
    )
  )

head(zi_interaction_time, 5)
zi_time_effects_fdr <- zi_interaction_time %>%
  group_by(scaled_age) %>%
  mutate(adj.p = p.adjust(p.value, method = "fdr")) %>%
  ungroup() %>%
  mutate(
    sig    = adj.p <= 0.05,
    rank_score = -log10(adj.p) * abs(estimate),
    sig_ci = conf.low > 0 | conf.high < 0,
    effect = case_when(
      conf.high < 0 & sig ~ "Improved SE",
      conf.low  > 0 & sig ~ "Reduced SE",
      TRUE                 ~ "No effect"
    ),
    # age_group = cut(
    #   scaled_age,
    #   breaks = c(-Inf, 0.20, 0.40, 0.60, 0.80, Inf),
    #   labels = c("Young", "Young-Mid", "Middle", "Mid-Old", "Old"),
    #   include.lowest = TRUE
    # ),
    # age_group = factor(age_group, 
    #                    levels = c("Young", "Young-Mid", 
    #                               "Middle", "Mid-Old", "Old"))
  ) %>%
  annotate_introns(gene_annotation, intron_length)


# count significant introns per group
count_plot <- zi_time_effects_fdr %>%
  filter(sig) %>%
  count(scaled_age, effect) %>%
  ggplot(aes(x = scaled_age, y = n, fill = effect)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = n), 
            position = position_dodge(width = 0.9), 
            vjust = -0.3, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("Improved SE" = colors[6],
                               "Reduced SE"  = colors[1])) +
  labs(
    title = "Number of introns with significant exercise effect per age group",
    x     = NULL,
    y     = "Number of introns",
    fill  = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top"
  )


# Identify introns significant at each age anchor
sig_targets <- zi_time_effects_fdr %>%
  filter(sig) %>%
  pull(target) %>%
  unique()

# plot them 
heatmap_plot <- zi_time_effects_fdr %>%
  filter(target %in% sig_targets) %>%
  mutate(
    # flip estimate back to SE scale
    estimate_SE = -estimate
  ) %>%
  ggplot(aes(x = scaled_age, y = gene_intron, fill = estimate_SE)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low      = colors[1],   # reduced SE
    mid      = "white",
    high     = colors[6],   # improved SE
    midpoint = 0
  ) +
  labs(
    title = "Exercise effect on splicing efficiency across age groups",
    x     = NULL,
    y     = NULL,
    fill  = "Effect on SE\n(Post - Pre)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(face = "bold")
  )


# Volcano plot per age group

volcano_age <- zi_time_effects_fdr %>%
  mutate(neg_log10_fdr = -log10(adj.p),
         estimate_SE   = -estimate) %>%
  ggplot(aes(x = estimate_SE, y = neg_log10_fdr, colour = effect)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), 
             linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0, 
             linetype = "dashed", colour = "grey50") +
  geom_text_repel(
    data = . %>% filter(sig) %>% 
      group_by(age_group) %>% 
      slice_max(abs(estimate_SE), n = 5),
    aes(label = gene_intron),
    size = 3, max.overlaps = Inf
  ) +
  scale_colour_manual(values = c(
    "Improved SE" = colors[6],
    "Reduced SE"  = colors[1],
    "No effect"   = "grey70"
  )) +
 # facet_wrap(~ age_group, ncol = 2) +
  coord_cartesian(xlim = c(-0.05, 0.1)) +
  labs(
    title = "Exercise effect on splicing efficiency across age groups",
    x     = "Effect size (Post - Pre SE)",
    y     = expression(-log[10](FDR)),
    colour = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    strip.text      = element_text(face = "bold"),
    legend.position = "top"
  )


# Which age anchor shows the strongest exercise effect
zi_time_effects_fdr %>%
  filter(sig) %>%
  group_by(scaled_age) %>%
  summarise(
    n_sig        = n_distinct(target),
    mean_effect  = mean(abs(estimate), na.rm = TRUE),
    median_effect = median(abs(estimate), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sig))


# Get the original age range
age_min <- min(metadata$age, na.rm = TRUE)
age_max <- max(metadata$age, na.rm = TRUE)

# Convert scaled_age back to real age
zi_pred_cond <- readRDS("data/zi_main_conditional.RDS")
zi_pred_cond <- zi_pred_cond %>%
  mutate(real_age = scaled_age * (age_max - age_min) + age_min)



# global conditional trajectory across age
# The investigates/predicts the effect of age and exercise on intron retention among imperfectly spliced introns
zi_pred_cond <- zi_pred_cond %>%
  mutate(SE = 1 - estimate) %>%
  group_by(real_age, time) %>%
  summarise(
    mean_SE = mean(SE, na.rm = TRUE),
    se      = sd(SE, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) 

zi_pred_cond %>%
  ggplot(aes(x = real_age, y = mean_SE, colour = time, fill = time)) +
  geom_ribbon(aes(ymin = mean_SE - se, ymax = mean_SE + se),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = c("PreExc" = colors[5], 
                                 "PostExc" = colors[1])) +
  scale_fill_manual(values   = c("PreExc" = colors[5], 
                                 "PostExc" = colors[1])) +
  labs(
    title    = "Conditional SE trajectory across age",
    subtitle = "Among introns that are not perfectly spliced",
    x        = "Scaled age",
    y        = "Mean splicing efficiency",
    colour   = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "top"
  )

# Compare conditional versus response trajectories
# This reveals how much the overall signal is driven by perfectly spliced pool versus the partially retained pool''
zi_pred_resp <- readRDS( "data/zi_main_predictions.RDS")

# Combine both prediction types
combined_pred <- bind_rows(
  zi_pred_resp %>% mutate(type = "Overall (response)"),
  zi_pred_cond %>% mutate(type = "Conditional (partial retention only)")
) %>%
  mutate(SE = 1 - estimate) %>%
  group_by(type, scaled_age, time) %>%
  summarise(
    mean_SE = mean(SE, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(combined_pred, aes(x = scaled_age, y = mean_SE, 
                          colour = time, linetype = type)) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = c("PreExc" = colors[5],
                                 "PostExc" = colors[1])) +
  labs(
    title    = "Overall vs conditional SE trajectory across age",
    subtitle = "Solid = overall response, dashed = conditional on partial retention",
    x        = "Scaled age",
    y        = "Mean splicing efficiency",
    colour   = NULL,
    linetype = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle   = element_text(hjust = 0.5),
    legend.position = "top"
  )


  divergent_introns <- zi_pred_resp %>%
  inner_join(zi_pred_cond, 
             by = c("target", "scaled_age", "time"),
             suffix = c("_response", "_conditional")) %>%
  mutate(
    divergence = abs(estimate_response - estimate_conditional)
  ) %>%
  group_by(target) %>%
  summarise(
    mean_divergence = mean(divergence, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_divergence)) %>%
  slice_head(n = 20) %>%
  pull(target)
  
  

  