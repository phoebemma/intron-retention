# 
# TRAINOME MODEL 2REFINED ZERO-INFLATED BETA REGRESSION PIPELINE
# Impact of Aging and Resistance Training on Intron Splicing Efficiency
#
# Key improvements over Trainome_model.R:
#   - study:time interaction added to account for study-specific exercise effects
#   - Age effect estimated after removing study-specific exercise responses
#   - ReLiEf sensitivity restructured as Young (18-30) vs Old (60-93) comparison
#   - Marginalisation over study AND study:time for clean population inference
#
# Models fitted:
#   1. ZI-Beta pooled  age + time + sex + study + study:time + (1|participant)
#   2. ZI-Beta ReLiEf Young vs Old, age_group * time + sex + (1|participant)
# 

library(dplyr)
library(tidyverse)
library(scales)
library(seqwrap)
library(glmmTMB)
library(marginaleffects)
library(purrr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(cowplot)
library(ggVennDiagram)
library(clusterProfiler)
library(org.Hs.eg.db)


# 
# COLOUR PALETTE
#

colors <- c("#d7191c", "#fdae61", "#ffffbf", "#abd9e9", "#2c7bb6",
            "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")

effect_colors <- c(
  "Improved SE" = colors[6],
  "Reduced SE"  = colors[1],
  "No effect"   = "grey70"
)


#
# LOAD AND COMPILE METADATA
# 

copd_metadata <- readRDS("data/copd_metadata.RDS") %>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)

vol_metadata <- readRDS("data/volume_metadata.RDS") %>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)

ct_metadata <- readRDS("data/contratrain_metadata.RDS") %>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)

ao_metadata <- readRDS("data/Alpha_Omega_metadata.RDS") %>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)

relief_metadata <- readRDS("data/Relief_metadata.RDS") %>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)

# Compile metadata — study as fixed factor, vol as reference level
metadata <- bind_rows(
  copd_metadata, vol_metadata, ct_metadata, ao_metadata, relief_metadata
) %>%
  mutate(
    age         = round(age, 0),
    sex         = factor(sex, levels = c("female", "male")),
    time        = factor(time, levels = c("PreExc", "PostExc")),
    scaled_age  = round(rescale(age), digits = 2),
    participant = paste0(study, "_", participant),
    study       = factor(study, levels = c("vol", "ct", "copd",
                                           "Alpha/Omega", "ReLiEf"))
  )

# Store age range for back-transformation
age_min <- min(metadata$age, na.rm = TRUE)
age_max <- max(metadata$age, na.rm = TRUE)

# Diagnose age-study overlap — important context for interpretation
age_study_overlap <- metadata %>%
  mutate(
    age_group = cut(scaled_age,
                    breaks = c(-Inf, 0.20, 0.40, 0.60, 0.80, Inf),
                    labels = c("Young", "Young-Mid", "Middle", "Mid-Old", "Old"))
  ) %>%
  count(study, age_group) %>%
  pivot_wider(names_from = age_group, values_from = n, values_fill = 0)

print(age_study_overlap)


#
# LOAD AND COMPILE SPLICING DATA
#

all_splice_df <- readRDS("data/copd_splicing_data.RDS") %>%
  inner_join(readRDS("data/volume_splicing_data.RDS"),      by = "transcript_ID") %>%
  inner_join(readRDS("data/contratrain_splicing_data.RDS"), by = "transcript_ID") %>%
  inner_join(readRDS("data/Alpha_Omega_splicing_data.RDS"), by = "transcript_ID") %>%
  inner_join(readRDS("data/Relief_splicing_data.RDS"),      by = "transcript_ID") %>%
  drop_na()

# Retain only samples present in both splice data and metadata
intersect_ids <- intersect(colnames(all_splice_df), metadata$seq_sample_id)

all_splice_df <- all_splice_df %>%
  subset(select = c("transcript_ID", intersect_ids)) %>%
  drop_na()

# Subset metadata to samples present in splice data
metadata <- metadata %>%
  filter(seq_sample_id %in% colnames(all_splice_df))

# Reorder splice columns to match metadata row order
all_splice_reordered <- all_splice_df[, c("transcript_ID", metadata$seq_sample_id)]


#
# PREPARE MODEL MATRICES
#

# Flip SE to retention scale: 0 = perfectly spliced (structural zero for zi)
zi_mat        <- all_splice_reordered
zi_mat[-1]    <- 1 - all_splice_reordered[-1]


# 
# ANNOTATION AND HELPER FUNCTIONS
#

gene_annotation <- readRDS("data/ensembl_gene_annotation.RDS")

gene_exp_df <- readRDS("data_new/gene_counts/batch_corrected_genecounts.RDS")

intron_length <- readr::read_tsv("data_new/Alpha_Omega_SpliceQ_outputs/A_102.tsv") %>%
  distinct(across(6:ncol(.)), .keep_all = TRUE) %>%
  mutate(
    transcript_ID  = paste0(transcript_ID, "_", intron_ID, "_", chr),
    intron_length  = abs((sj3start - sj5end) + 1),
    intron_ID      = as.character(intron_ID),
    number_introns = n()
  ) %>%
  dplyr::select(transcript_ID, intron_ID, gene_ID, intron_length, number_introns)

# Annotates outputs with gene names and intron metadata
# flip = TRUE reverses estimate from retention to SE scale
annotate_introns <- function(df, gene_annotation, intron_length, flip = TRUE) {
  df %>%
    mutate(
      estimate       = if (flip) -estimate else estimate,
      conf.low_se    = if (flip) -conf.high else conf.low,
      conf.high_se   = if (flip) -conf.low else conf.high,
      transcript_ID  = str_split(target, "_", simplify = TRUE)[, 1]
    ) %>%
    inner_join(gene_annotation,
               by = c("transcript_ID" = "ensembl_transcript_id_version")) %>%
    separate(target, into = c(NA, "intron_ID", NA),
             sep = "_", remove = FALSE) %>%
    inner_join(intron_length %>%
                 mutate(intron_ID = as.character(intron_ID)),
               by = c("intron_ID",
                      "ensembl_gene_id_version" = "gene_ID",
                      "target" = "transcript_ID")) %>%
    mutate(
      gene_label  = ifelse(
        is.na(external_gene_name) | external_gene_name == "",
        ensembl_gene_id, external_gene_name
      ),
      gene_intron = paste(gene_label, intron_ID, sep = " : ")
    ) %>%
    {if ("rank_score" %in% colnames(.)) arrange(., desc(rank_score)) else .}
}


# 
# POOLED MODEL  ZI-BETA WITH STUDY:TIME INTERACTION
#
# Key change from Trainome_model.R:
#   study:time added to both formula and ziformula
#   This removes study-specific exercise effects before estimating age slope
#   Prevents study-age confounding from masking the true age trajectory
# 

zi_main_container <- seqwrap_compose(
  data       = zi_mat,
  metadata   = metadata,
  samplename = "seq_sample_id",
  modelfun   = glmmTMB::glmmTMB,
  arguments  = alist(
    formula   = y ~ splines::ns(scaled_age, df = 4) + time + sex +
                    study + study:time +
                    (1 | participant),
    ziformula = ~ splines::ns(scaled_age, df = 4) + time +
                    study + study:time,
    family    = glmmTMB::beta_family(link = "logit")
  )
)

zi_main_results <- seqwrap(
  zi_main_container,
  return_models = TRUE,
  cores         = 10
)

saveRDS(zi_main_results, "data/zi_main_results.RDS")
# zi_main_results <- readRDS("data/zi_main_results.RDS")


# 
# MARGINAL EFFECTS POOLED MODEL
#
# Marginalise over ALL studies AND study:time combinations
# by including study in datagrid and averaging  this gives a
# population-level estimate free from any single study's exercise response
#

valid_zi_main  <- compact(zi_main_results@models)
age_anchors    <- seq(0, 1, by = 0.10)

#  Predicted SE trajectories (response scale) 
zi_predictions <- map_dfr(
  valid_zi_main,
  function(mod) {
    tryCatch({
      predictions(
        mod,
        type       = "response",
        re_formula = NA,
        vcov       = TRUE,
        newdata    = datagrid(
          scaled_age = seq(0, 1, length.out = 100),
          time       = c("PreExc", "PostExc"),
          study      = unique(metadata$study)
        )
      ) %>%
        # Average over studies to get population-level estimate
        group_by(scaled_age, time) %>%
        summarise(
          estimate  = mean(estimate),
          conf.low  = mean(conf.low),
          conf.high = mean(conf.high),
          .groups   = "drop"
        )
    }, error = function(e) NULL)
  },
  .id = "target"
)


saveRDS(zi_predictions,  "data/zi_predictions.RDS")


#  Probability of perfect splicing (zi component) 
zi_zprob <- map_dfr(
  valid_zi_main,
  function(mod) {
    tryCatch({
      predictions(
        mod,
        type       = "zprob",
        re_formula = NA,
        vcov       = TRUE,
        newdata    = datagrid(
          scaled_age = seq(0, 1, length.out = 100),
          time       = c("PreExc", "PostExc"),
          study      = unique(metadata$study)
        )
      ) %>%
        group_by(scaled_age, time) %>%
        summarise(
          estimate  = mean(estimate),
          conf.low  = mean(conf.low),
          conf.high = mean(conf.high),
          .groups   = "drop"
        )
    }, error = function(e) NULL)
  },
  .id = "target"
)

saveRDS(zi_zprob, "data/zi_zprob.RDS")

#  Degree of retention among partially retained introns (beta component)
zi_conditional <- map_dfr(
  valid_zi_main,
  function(mod) {
    tryCatch({
      predictions(
        mod,
        type       = "conditional",
        re_formula = NA,
        vcov       = TRUE,
        newdata    = datagrid(
          scaled_age = seq(0, 1, length.out = 100),
          time       = c("PreExc", "PostExc"),
          study      = unique(metadata$study)
        )
      ) %>%
        group_by(scaled_age, time) %>%
        summarise(
          estimate  = mean(estimate),
          conf.low  = mean(conf.low),
          conf.high = mean(conf.high),
          .groups   = "drop"
        )
    }, error = function(e) NULL)
  },
  .id = "target"
)


saveRDS(zi_conditional,  "data/zi_conditional.RDS")

# Age slopes at 0.10 increments, separately per timepoint 
zi_age_slopes <- map_dfr(
  valid_zi_main,
  function(mod) {
    tryCatch({
      avg_slopes(
        mod,
        variables  = "scaled_age",
        by         = "time",
        newdata    = datagrid(
          scaled_age = age_anchors,
          time       = c("PreExc", "PostExc"),
          study      = unique(metadata$study)
        ),
        type       = "response",
        re_formula = NA
      )
    }, error = function(e) NULL)
  },
  .id = "target"
)

saveRDS(zi_age_slopes,   "data/zi_age_slopes.RDS")

#  Exercise effect (PostExc - PreExc) at each age anchor 
zi_time_effects <- map_dfr(
  valid_zi_main,
  function(mod) {
    tryCatch({
      avg_comparisons(
        mod,
        variables  = "time",
        by         = "scaled_age",
        newdata    = datagrid(
          scaled_age = age_anchors,
          study      = unique(metadata$study)
        ),
        type       = "response",
        re_formula = NA
      )
    }, error = function(e) NULL)
  },
  .id = "target"
)


saveRDS(zi_time_effects, "data/zi_time_effects.RDS")


# 
# SECTION 8: FDR FILTERING  POOLED MODEL
#

# Age slopes FDR within each timepoint ---
zi_age_slopes_fdr <- zi_age_slopes %>%
  group_by(time) %>%
  mutate(adj.p = p.adjust(p.value, method = "fdr")) %>%
  ungroup() %>%
  mutate(
   # real_age      = scaled_age * (age_max - age_min) + age_min,
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


# =============================================================================
# SECTION 9: ReLiEf SENSITIVITY — YOUNG VS OLD
#
# ReLiEf is the only study with participants at both ends of the age spectrum
# Young: 18-30, Old: 60-93
# Participants aged 30-60 excluded to create clean young vs old comparison
# This avoids the spline extrapolation problem across the age gap
# =============================================================================

# --- 9a. Subset ReLiEf metadata ---
metadata_relief <- metadata %>%
  filter(study == "ReLiEf") %>%
  filter(age <= 30 | age >= 60) %>%
  mutate(
    age_group = factor(
      ifelse(age <= 30, "Young", "Old"),
      levels = c("Young", "Old")
    )
  )

cat("ReLiEf participants:\n")
print(table(metadata_relief$age_group, metadata_relief$time))

# --- 9b. Subset splice data to ReLiEf samples ---
relief_ids    <- metadata_relief$seq_sample_id

splice_relief <- all_splice_reordered %>%
  dplyr::select(transcript_ID, any_of(relief_ids))

splice_relief <- splice_relief[, c("transcript_ID",
                                    metadata_relief$seq_sample_id)]

zi_mat_relief      <- splice_relief
zi_mat_relief[-1]  <- 1 - splice_relief[-1]

# --- 9c. Fit ReLiEf model ---
# age_group * time tests whether exercise response differs by age group
# No study random effect needed — single study
zi_relief_container <- seqwrap_compose(
  data       = zi_mat_relief,
  metadata   = metadata_relief,
  samplename = "seq_sample_id",
  modelfun   = glmmTMB::glmmTMB,
  arguments  = alist(
    formula   = y ~ age_group * time + sex +
                    (1 | participant),
    ziformula = ~ age_group * time,
    family    = glmmTMB::beta_family(link = "logit")
  )
)

zi_relief_results <- seqwrap(
  zi_relief_container,
  return_models = TRUE,
  cores         = 10
)

saveRDS(zi_relief_results, "data/zi_relief_young_old_results.RDS")
# zi_relief_results <- readRDS("data/zi_relief_young_old_results.RDS")


# =============================================================================
# SECTION 10: MARGINAL EFFECTS — ReLiEf MODEL
# =============================================================================

valid_relief <- compact(zi_relief_results@models)

# --- 10a. Predicted SE per age group and timepoint ---
zi_relief_predictions <- map_dfr(
  valid_relief,
  function(mod) {
    tryCatch({
      predictions(
        mod,
        type       = "response",
        re_formula = NA,
        vcov       = TRUE,
        newdata    = datagrid(
          age_group = c("Young", "Old"),
          time      = c("PreExc", "PostExc")
        )
      )
    }, error = function(e) NULL)
  },
  .id = "target"
)

# --- 10b. Exercise effect within each age group ---
zi_relief_time_effects <- map_dfr(
  valid_relief,
  function(mod) {
    tryCatch({
      avg_comparisons(
        mod,
        variables  = "time",
        by         = "age_group",
        newdata    = datagrid(
          age_group = c("Young", "Old")
        ),
        type       = "response",
        re_formula = NA
      )
    }, error = function(e) NULL)
  },
  .id = "target"
)

# --- 10c. Age group effect within each timepoint ---
zi_relief_age_effects <- map_dfr(
  valid_relief,
  function(mod) {
    tryCatch({
      avg_comparisons(
        mod,
        variables  = "age_group",
        by         = "time",
        newdata    = datagrid(
          time = c("PreExc", "PostExc")
        ),
        type       = "response",
        re_formula = NA
      )
    }, error = function(e) NULL)
  },
  .id = "target"
)

# --- 10d. zprob predictions for ReLiEf ---
zi_relief_zprob <- map_dfr(
  valid_relief,
  function(mod) {
    tryCatch({
      predictions(
        mod,
        type       = "zprob",
        re_formula = NA,
        vcov       = TRUE,
        newdata    = datagrid(
          age_group = c("Young", "Old"),
          time      = c("PreExc", "PostExc")
        )
      )
    }, error = function(e) NULL)
  },
  .id = "target"
)

saveRDS(zi_relief_predictions,  "data/zi_relief_predictions.RDS")
saveRDS(zi_relief_time_effects, "data/zi_relief_time_effects.RDS")
saveRDS(zi_relief_age_effects,  "data/zi_relief_age_effects.RDS")
saveRDS(zi_relief_zprob,        "data/zi_relief_zprob.RDS")


# =============================================================================
# SECTION 11: FDR FILTERING — ReLiEf MODEL
# =============================================================================

# Exercise effect within each age group
zi_relief_time_fdr <- zi_relief_time_effects %>%
  group_by(age_group) %>%
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

# Age group effect within each timepoint
zi_relief_age_fdr <- zi_relief_age_effects %>%
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


# 
# VISUALISATIONS POOLED MODEL
#

#  Global SE trajectory across age 
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
    legend.position = "top"
  )

print(plot_global_trajectory)
# ggsave("Figures/global_trajectory_v2.svg", plot_global_trajectory,
#         width = 10, height = 6)


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
    legend.position = "top"
  )

print(plot_conditional)
# ggsave("Figures/conditional_trajectory_v2.svg", plot_conditional,
#         width = 10, height = 6)


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
      slice_max(abs(estimate), n = 5),
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
# ggsave("Figures/volcano_age_v2.svg", plot_volcano_age,
#         width = 12, height = 6)


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
# ggsave("Figures/exercise_trajectory_v2.svg", plot_exercise_trajectory,
#         width = 10, height = 6)


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
    plot.subtitle   = element_text(hjust = 0.5),
    legend.position = "top"
  )

print(plot_exercise_traj_sig)
# ggsave("Figures/exercise_traj_sig_v2.svg", plot_exercise_traj_sig,
#         width = 10, height = 6)


#  Overlap between PreExc and PostExc age-significant introns 
sig_preexc <- zi_age_slopes_fdr %>%
  filter(time == "PreExc", sig) %>%
  pull(target) %>% unique()

sig_postexc <- zi_age_slopes_fdr %>%
  filter(time == "PostExc", sig) %>%
  pull(target) %>% unique()

plot_venn_age <- ggVennDiagram(
  list(PreExc = sig_preexc, PostExc = sig_postexc)
) +
  labs(title = "Age-associated introns: rest vs post-exercise") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(plot_venn_age)
# ggsave("Figures/venn_age_prepost_v2.svg", plot_venn_age,
#         width = 8, height = 6)


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
# ggsave("Figures/traj_top_aging_v2.svg", plot_traj_top_aging,
#         width = 12, height = 8)


# =============================================================================
# SECTION 13: VISUALISATIONS — ReLiEf YOUNG VS OLD
# =============================================================================

# --- 13a. Mean SE per age group and timepoint ---
plot_relief_bar <- zi_relief_predictions %>%
  mutate(SE = 1 - estimate) %>%
  group_by(age_group, time) %>%
  summarise(
    mean_SE = mean(SE, na.rm = TRUE),
    se      = sd(SE, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = age_group, y = mean_SE, fill = time)) +
  geom_col(position = "dodge", colour = "white") +
  geom_errorbar(
    aes(ymin = mean_SE - se, ymax = mean_SE + se),
    position = position_dodge(width = 0.9),
    width = 0.25
  ) +
  scale_fill_manual(values = c("PreExc" = colors[5],
                               "PostExc" = colors[1])) +
  labs(
    title = "Mean SE by age group and timepoint — ReLiEf",
    x     = NULL,
    y     = "Mean splicing efficiency",
    fill  = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top"
  )

print(plot_relief_bar)
# ggsave("Figures/relief_bar_v2.svg", plot_relief_bar,
#         width = 8, height = 6)


# --- 13b. Volcano — exercise effect within each age group ---
plot_relief_volcano_exercise <- zi_relief_time_fdr %>%
  ggplot(aes(x = estimate, y = neg_log10_fdr, colour = effect)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0,
             linetype = "dashed", colour = "grey50") +
  geom_text_repel(
    data = . %>% filter(sig) %>%
      group_by(age_group) %>%
      slice_max(abs(estimate), n = 5),
    aes(label = gene_intron),
    size = 3, max.overlaps = Inf
  ) +
  scale_colour_manual(values = effect_colors) +
  facet_wrap(~ age_group, ncol = 2) +
  labs(
    title    = "ReLiEf: exercise effect on SE by age group",
    subtitle = "Young (18-30) vs Old (60-93)",
    x        = "Effect size (PostExc - PreExc SE)",
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

print(plot_relief_volcano_exercise)
# ggsave("Figures/relief_volcano_exercise_v2.svg", plot_relief_volcano_exercise,
#         width = 12, height = 6)


# --- 13c. Volcano — age group effect at rest and post-exercise ---
plot_relief_volcano_age <- zi_relief_age_fdr %>%
  ggplot(aes(x = estimate, y = neg_log10_fdr, colour = effect)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0,
             linetype = "dashed", colour = "grey50") +
  geom_text_repel(
    data = . %>% filter(sig) %>%
      group_by(time) %>%
      slice_max(abs(estimate), n = 5),
    aes(label = gene_intron),
    size = 3, max.overlaps = Inf
  ) +
  scale_colour_manual(values = effect_colors) +
  facet_wrap(~ time, ncol = 2) +
  labs(
    title    = "ReLiEf: age group effect on SE at rest and after exercise",
    subtitle = "Old vs Young",
    x        = "Effect size (Old - Young SE)",
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

print(plot_relief_volcano_age)
# ggsave("Figures/relief_volcano_age_v2.svg", plot_relief_volcano_age,
#         width = 12, height = 6)


# --- 13d. Overlap: exercise-responsive introns young vs old ---
sig_young <- zi_relief_time_fdr %>%
  filter(age_group == "Young", sig) %>%
  pull(target) %>% unique()

sig_old <- zi_relief_time_fdr %>%
  filter(age_group == "Old", sig) %>%
  pull(target) %>% unique()

plot_relief_venn_exercise <- ggVennDiagram(
  list(Young = sig_young, Old = sig_old)
) +
  labs(title = "Exercise-responsive introns: Young vs Old") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(plot_relief_venn_exercise)
# ggsave("Figures/relief_venn_exercise_v2.svg", plot_relief_venn_exercise,
#         width = 8, height = 6)


# --- 13e. Overlap: age-associated introns vs pooled model ---
sig_pooled <- zi_age_slopes_fdr %>%
  filter(sig) %>%
  pull(target) %>% unique()

sig_relief_age <- zi_relief_age_fdr %>%
  filter(sig) %>%
  pull(target) %>% unique()

plot_venn_validation <- ggVennDiagram(
  list(`Pooled model` = sig_pooled, `ReLiEf Young vs Old` = sig_relief_age)
) +
  labs(title = "Validation: pooled model vs ReLiEf sensitivity") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(plot_venn_validation)
# ggsave("Figures/venn_validation_v2.svg", plot_venn_validation,
#         width = 8, height = 6)


# --- 13f. Trajectory plots for top introns in ReLiEf ---
top_relief_targets <- zi_relief_time_fdr %>%
  filter(sig) %>%
  slice_max(abs(estimate), n = 6) %>%
  pull(target)

plot_traj_relief <- zi_relief_predictions %>%
  filter(target %in% top_relief_targets) %>%
  mutate(SE = 1 - estimate) %>%
  left_join(
    zi_relief_time_fdr %>%
      dplyr::select(target, gene_intron) %>% distinct(),
    by = "target"
  ) %>%
  ggplot(aes(x = age_group, y = SE, colour = time,
             group = time, fill = time)) +
  geom_point(size = 3) +
  geom_line(linewidth = 0.8) +
  geom_errorbar(
    aes(ymin = 1 - conf.high, ymax = 1 - conf.low),
    width = 0.1
  ) +
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
# ggsave("Figures/traj_relief_top_v2.svg", plot_traj_relief,
#         width = 12, height = 8)


# =============================================================================
# SECTION 14: COMBINED PUBLICATION FIGURES
# =============================================================================

# --- Figure 1: Pooled model main results ---
figure_pooled <- (plot_global_trajectory | plot_venn_age) /
  plot_volcano_age /
  (plot_exercise_trajectory | plot_exercise_traj_sig) +
  plot_annotation(
    title    = "Pooled analysis: age and exercise effects on splicing efficiency",
    tag_levels = "A"
  ) &
  theme(
    plot.tag   = element_text(face = "bold", size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  )

print(figure_pooled)
# ggsave("Figures/figure_pooled_v2.svg", figure_pooled,
#         width = 16, height = 18)


# --- Figure 2: ReLiEf sensitivity ---
figure_relief <- (plot_relief_bar | plot_relief_venn_exercise) /
  plot_relief_volcano_exercise /
  plot_relief_volcano_age /
  (plot_traj_relief | plot_venn_validation) +
  plot_annotation(
    title      = "ReLiEf sensitivity: Young vs Old",
    tag_levels = "A"
  ) &
  theme(
    plot.tag   = element_text(face = "bold", size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  )

print(figure_relief)
# ggsave("Figures/figure_relief_v2.svg", figure_relief,
#         width = 16, height = 22)
