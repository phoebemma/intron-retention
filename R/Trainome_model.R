# 
# Impact of Aging and Resistance Training on Intron Splicing Efficiency
# Models fitted:
#   1. ZI-Beta main effects   age + time + sex + study (fixed) + (1|participant)
#   2. ZI-Beta interaction    age * time + sex + study (fixed) + (1|participant)
#   3. ZI-Beta ReLiEf only   sensitivity analysis within a single study
#
# Key design decisions:
#   - Study included as FIXED effect to account for age-study confounding
#   - SE scores flipped (1 - SE) so 0 = perfect splicing (structural zero)
#   - Natural splines (df=4) for non-linear age modelling
#   - Marginal effects marginalised over all studies for population inference
#   - Age increments of 0.20 for avg_slopes 


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

# DEFINE COLOUR PALETTE (for visualisation)

colors <- c("#d7191c", "#fdae61", "#ffffbf", "#abd9e9", "#2c7bb6",
            "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")

effect_colors <- c(
  "Improved SE"  = colors[6],
  "Reduced SE"   = colors[1],
  "No effect"    = "grey70"
)


# LOAD AND COMPILE METADATA

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


# Compile all metadata and derive scaled age across pooled sample
metadata <- bind_rows(
  copd_metadata, vol_metadata, ct_metadata, ao_metadata, relief_metadata
) %>%
  mutate(
    age         = round(age, 0),
    sex         = factor(sex, levels = c("female", "male")),
    time        = factor(time, levels = c("PreExc", "PostExc")),
    scaled_age  = round(rescale(age), digits = 2),
    participant = paste0(study, "_", participant),
    # Study as factor with vol as reference level (largest young group)
    study       = factor(study, levels = c("vol", "ct", "copd",
                                           "Alpha/Omega", "ReLiEf"))
  )


# Store age range for back-transformation later
age_min <- min(metadata$age, na.rm = TRUE)
age_max <- max(metadata$age, na.rm = TRUE)


# 
#  LOAD AND COMPILE SPLICING DATA
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

# Flip SE scores: retention = 1 - SE
# Perfect splicing (SE=1) becomes structural zero (retention=0)
# Complete retention (SE=0) becomes retention=1
zi_mat <- all_splice_reordered
zi_mat[-1] <- 1 - all_splice_reordered[-1]


# 
#  LOAD FILES AND ANNOTATION HELPER FUNCTION
# 

# Load annotation files
gene_annotation <- readRDS("data/ensembl_gene_annotation.RDS")

intron_length <- readr::read_tsv("data_new/Alpha_Omega_SpliceQ_outputs/A_102.tsv") %>%
  distinct(across(6:ncol(.)), .keep_all = T) %>%
  mutate(transcript_ID = paste0(transcript_ID, "_", intron_ID, "_", chr),
         intron_length = abs((sj3start - sj5end) + 1)) %>%
  group_by(gene_ID) %>%
  mutate(number_introns = n()) %>%
  ungroup()

# Annotates model outputs with gene names, intron IDs, and intron length
# flip = TRUE reverses estimate direction back to SE scale for interpretation
annotate_introns <- function(df, gene_annotation, intron_length, flip = TRUE) {
  df %>%
    mutate(
      estimate      = if (flip) -estimate else estimate,
      conf.low_flip = if (flip) -conf.high else conf.low,
      conf.high_flip = if (flip) -conf.low else conf.high,
      transcript_ID = str_split(target, "_", simplify = TRUE)[, 1]
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
# MODEL 1  ZI-BETA MAIN EFFECTS (STUDY AS FIXED EFFECT)
# Biological question: Does aging and/or resistance training change
# intron retention, after accounting for study-specific baselines?
# 

zi_main_container <- seqwrap_compose(
  data       = zi_mat,
  metadata   = metadata,
  samplename = "seq_sample_id",
  modelfun   = glmmTMB::glmmTMB,
  arguments  = alist(
    formula   = y ~ splines::ns(scaled_age, df = 4) + time + sex + study +
                    (1 | participant),
    ziformula = ~ splines::ns(scaled_age, df = 4) + time + study,
    family    = glmmTMB::beta_family(link = "logit")
  )
)

zi_main_results <- seqwrap(
  zi_main_container,
  return_models = TRUE,
  cores         = 10
)

# saveRDS(zi_main_results, "data/zi_main_results.RDS")
# zi_main_results <- readRDS("data/zi_main_results.RDS")


# 
# MODEL 2 ZI-BETA INTERACTION (AGE x TIME)
# Biological question: Does the splicing response to exercise differ
# by age, after accounting for study-specific baselines?
# 

zi_interaction_container <- seqwrap_compose(
  data       = zi_mat,
  metadata   = metadata,
  samplename = "seq_sample_id",
  modelfun   = glmmTMB::glmmTMB,
  arguments  = alist(
    formula   = y ~ splines::ns(scaled_age, df = 4) * time + sex + study +
                    (1 | participant),
    ziformula = ~ splines::ns(scaled_age, df = 4) * time + study,
    family    = glmmTMB::beta_family(link = "logit")
  )
)

zi_interaction_results <- seqwrap(
  zi_interaction_container,
  return_models = TRUE,
  cores         = 10
)

# saveRDS(zi_interaction_results, "data/zi_interaction_results.RDS")
# zi_interaction_results <- readRDS("data/zi_interaction_results.RDS")


# 
#  MARGINAL EFFECTS MAIN EFFECTS MODEL
# All predictions marginalized over studies by including all studies
# in data grid and averaging, ensuring population-level inference
# 

valid_zi_main <- compact(zi_main_results@models)

# Age increments of 0.20 for finer resolution than previous 0.25
age_anchors <- seq(0, 1, by = 0.10)

# Predicted retention trajectories (overall response)
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

#  Probability of perfect splicing across age (zi component) ---
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

# Degree of retention among partially retained introns (beta component) ---
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

#  Age slopes at 0.20 increments, separately at each timepoint 
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

#  Exercise effect (PostExc vs PreExc) at each age anchor 
zi_time_effects_main <- map_dfr(
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

saveRDS(zi_predictions,        "data/zi_main_predictions.RDS")
saveRDS(zi_zprob,              "data/zi_main_zprob.RDS")
saveRDS(zi_conditional,        "data/zi_main_conditional.RDS")
saveRDS(zi_age_slopes,         "data/zi_main_age_slopes.RDS")
saveRDS(zi_time_effects_main,  "data/zi_main_time_effects.RDS")


#
# MARGINAL EFFECTS FOR INTERACTION MODEL
# 

valid_zi_interaction <- compact(zi_interaction_results@models)

# Full predicted trajectories 
zi_int_predictions <- map_dfr(
  valid_zi_interaction,
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

# Age-dependent exercise effect at each age anchor ---
zi_time_effects_interaction <- map_dfr(
  valid_zi_interaction,
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

saveRDS(zi_int_predictions,          "data/zi_interaction_predictions.RDS")
saveRDS(zi_time_effects_interaction, "data/zi_interaction_time_effects.RDS")


# 
# FDR FILTERING  MAIN EFFECTS MODEL
# 

# --- 10a. Age slopes: FDR within each timepoint ---
zi_age_slopes_fdr <- zi_age_slopes %>%
  group_by(time) %>%
  mutate(adj.p = p.adjust(p.value, method = "fdr")) %>%
  ungroup() %>%
  mutate(
    real_age   = scaled_age * (age_max - age_min) + age_min,
    sig        = adj.p <= 0.05,
    sig_ci     = conf.low > 0 | conf.high < 0,
    rank_score = -log10(adj.p) * abs(estimate),
    effect     = case_when(
      conf.high < 0 & sig & sig_ci ~ "Improved SE",
      conf.low  > 0 & sig & sig_ci ~ "Reduced SE",
      conf.high < 0 & sig          ~ "Improved SE",
      conf.low  > 0 & sig          ~ "Reduced SE",
      TRUE                          ~ "No effect"
    )
  ) %>%
  annotate_introns(gene_annotation, intron_length, flip = TRUE)

# --- 10b. Exercise effects: FDR within each age anchor ---
zi_time_effects_fdr <- zi_time_effects_main %>%
  group_by(scaled_age) %>%
  mutate(adj.p = p.adjust(p.value, method = "fdr")) %>%
  ungroup() %>%
  mutate(
    real_age      = scaled_age * (age_max - age_min) + age_min,
    sig           = adj.p <= 0.05,
    sig_ci        = conf.low > 0 | conf.high < 0,
    neg_log10_fdr = -log10(adj.p),
    rank_score    = -log10(adj.p) * abs(estimate),
    age_group     = case_when(
      scaled_age == 0.0 ~ "Young",
      scaled_age == 0.2 ~ "Young-Mid",
      scaled_age == 0.4 ~ "Middle",
      scaled_age == 0.6 ~ "Mid-Old",
      scaled_age == 0.8 ~ "Old",
      scaled_age == 1.0 ~ "Oldest"
    ),
    age_group = factor(age_group,
                       levels = c("Young", "Young-Mid", "Middle",
                                  "Mid-Old", "Old", "Oldest")),
    effect = case_when(
      conf.high < 0 & sig ~ "Improved SE",
      conf.low  > 0 & sig ~ "Reduced SE",
      TRUE                 ~ "No effect"
    )
  ) %>%
  annotate_introns(gene_annotation, intron_length, flip = TRUE)

# --- 10c. Interaction model time effects: FDR within each age anchor ---
zi_interaction_fdr <- zi_time_effects_interaction %>%
  group_by(scaled_age) %>%
  mutate(adj.p = p.adjust(p.value, method = "fdr")) %>%
  ungroup() %>%
  mutate(
    real_age      = scaled_age * (age_max - age_min) + age_min,
    sig           = adj.p <= 0.05,
    sig_ci        = conf.low > 0 | conf.high < 0,
    neg_log10_fdr = -log10(adj.p),
    rank_score    = -log10(adj.p) * abs(estimate),
    age_group     = case_when(
      scaled_age == 0.0 ~ "Young",
      scaled_age == 0.2 ~ "Young-Mid",
      scaled_age == 0.4 ~ "Middle",
      scaled_age == 0.6 ~ "Mid-Old",
      scaled_age == 0.8 ~ "Old",
      scaled_age == 1.0 ~ "Oldest"
    ),
    age_group = factor(age_group,
                       levels = c("Young", "Young-Mid", "Middle",
                                  "Mid-Old", "Old", "Oldest")),
    effect = case_when(
      conf.high < 0 & sig ~ "Improved SE",
      conf.low  > 0 & sig ~ "Reduced SE",
      TRUE                 ~ "No effect"
    )
  ) %>%
  annotate_introns(gene_annotation, intron_length, flip = TRUE)


# =============================================================================
# SECTION 11: SENSITIVITY ANALYSIS — ReLiEf ONLY
# ReLiEf is the only study spanning sufficient age range for within-study
# age effect estimation, free from age-study confounding
# =============================================================================

# --- 11a. Subset metadata to ReLiEf ---
metadata_relief <- metadata %>%
  filter(study == "ReLiEf") %>%
  # Recompute scaled_age within ReLiEf age range
  mutate(scaled_age = round(rescale(age), digits = 2))

age_min_relief <- min(metadata_relief$age, na.rm = TRUE)
age_max_relief <- max(metadata_relief$age, na.rm = TRUE)

# --- 11b. Subset splice data to ReLiEf samples ---
relief_ids    <- metadata_relief$seq_sample_id

splice_relief <- all_splice_reordered %>%
  dplyr::select(transcript_ID, any_of(relief_ids))

# Reorder columns to match metadata row order
splice_relief <- splice_relief[, c("transcript_ID",
                                    metadata_relief$seq_sample_id)]

# Flip to retention scale
zi_mat_relief        <- splice_relief
zi_mat_relief[-1]    <- 1 - splice_relief[-1]

# --- 11c. Fit ZI-Beta model within ReLiEf ---
# No study random or fixed effect needed — single study
# df reduced to 3 to account for smaller sample and narrower age range
zi_relief_container <- seqwrap_compose(
  data       = zi_mat_relief,
  metadata   = metadata_relief,
  samplename = "seq_sample_id",
  modelfun   = glmmTMB::glmmTMB,
  arguments  = alist(
    formula   = y ~ splines::ns(scaled_age, df = 3) + time + sex +
                    (1 | participant),
    ziformula = ~ splines::ns(scaled_age, df = 3) + time,
    family    = glmmTMB::beta_family(link = "logit")
  )
)

zi_relief_results <- seqwrap(
  zi_relief_container,
  return_models = TRUE,
  cores         = 10
)

saveRDS(zi_relief_results, "data/zi_relief_results.RDS")
# zi_relief_results <- readRDS("data/zi_relief_results.RDS")

# --- 11d. Marginal effects — ReLiEf only ---
valid_relief <- compact(zi_relief_results@models)

# Predicted trajectories
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
          scaled_age = seq(0, 1, length.out = 100),
          time       = c("PreExc", "PostExc")
        )
      )
    }, error = function(e) NULL)
  },
  .id = "target"
)

# Age slopes at 0.20 increments
zi_relief_slopes <- map_dfr(
  valid_relief,
  function(mod) {
    tryCatch({
      avg_slopes(
        mod,
        variables  = "scaled_age",
        by         = "time",
        newdata    = datagrid(
          scaled_age = age_anchors,
          time       = c("PreExc", "PostExc")
        ),
        type       = "response",
        re_formula = NA
      )
    }, error = function(e) NULL)
  },
  .id = "target"
)

# Exercise effect at each age anchor
zi_relief_time_effects <- map_dfr(
  valid_relief,
  function(mod) {
    tryCatch({
      avg_comparisons(
        mod,
        variables  = "time",
        by         = "scaled_age",
        newdata    = datagrid(
          scaled_age = age_anchors
        ),
        type       = "response",
        re_formula = NA
      )
    }, error = function(e) NULL)
  },
  .id = "target"
)

saveRDS(zi_relief_predictions,   "data/zi_relief_predictions.RDS")
saveRDS(zi_relief_slopes,        "data/zi_relief_slopes.RDS")
saveRDS(zi_relief_time_effects,  "data/zi_relief_time_effects.RDS")

# --- 11e. FDR filtering — ReLiEf ---
zi_relief_slopes_fdr <- zi_relief_slopes %>%
  group_by(time) %>%
  mutate(adj.p = p.adjust(p.value, method = "fdr")) %>%
  ungroup() %>%
  mutate(
    real_age   = scaled_age * (age_max_relief - age_min_relief) + age_min_relief,
    sig        = adj.p <= 0.05,
    sig_ci     = conf.low > 0 | conf.high < 0,
    rank_score = -log10(adj.p) * abs(estimate),
    effect     = case_when(
      conf.high < 0 & sig ~ "Improved SE",
      conf.low  > 0 & sig ~ "Reduced SE",
      TRUE                 ~ "No effect"
    )
  ) %>%
  annotate_introns(gene_annotation, intron_length, flip = TRUE)

zi_relief_time_fdr <- zi_relief_time_effects %>%
  group_by(scaled_age) %>%
  mutate(adj.p = p.adjust(p.value, method = "fdr")) %>%
  ungroup() %>%
  mutate(
    real_age      = scaled_age * (age_max_relief - age_min_relief) + age_min_relief,
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
# SECTION 12: VISUALISATIONS — MAIN EFFECTS MODEL
# =============================================================================

# --- 12a. Global SE trajectory across age (main model) ---
plot_global_trajectory <- zi_predictions %>%
  mutate(SE = 1 - estimate,
         CI_low  = 1 - conf.high,
         CI_high = 1 - conf.low,
         real_age = scaled_age * (age_max - age_min) + age_min) %>%
  group_by(real_age, time) %>%
  summarise(
    mean_SE  = mean(SE),
    CI_low   = mean(CI_low),
    CI_high  = mean(CI_high),
    .groups  = "drop"
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
    subtitle = "Solid = spline fit, dashed = linear fit",
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
# ggsave("Figures/global_trajectory_main.svg", plot_global_trajectory,
#         width = 10, height = 6)


# --- 12b. Volcano plot — age slopes per timepoint ---
plot_volcano_age <- zi_age_slopes_fdr %>%
  mutate(neg_log10_fdr = -log10(adj.p)) %>%
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
    title    = "Age effect on splicing efficiency — PreExc vs PostExc",
    subtitle = "Positive = improved SE with age",
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
# ggsave("Figures/volcano_age_slopes.svg", plot_volcano_age,
#         width = 12, height = 6)


# --- 12c. Count of significant introns per age group and timepoint ---
plot_sig_counts <- zi_time_effects_fdr %>%
  filter(sig) %>%
  count(age_group, effect) %>%
  ggplot(aes(x = age_group, y = n, fill = effect)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = n),
            position = position_dodge(width = 0.9),
            vjust = -0.3, size = 4, fontface = "bold") +
  scale_fill_manual(values = effect_colors) +
  labs(
    title = "Significant exercise-associated introns per age group",
    x     = NULL,
    y     = "Number of introns",
    fill  = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top"
  )

print(plot_sig_counts)
# ggsave("Figures/sig_counts_age_group.svg", plot_sig_counts,
#         width = 10, height = 6)


# --- 12d. Volcano per age group — exercise effect ---
plot_volcano_exercise <- zi_time_effects_fdr %>%
  ggplot(aes(x = estimate, y = neg_log10_fdr, colour = effect)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0,
             linetype = "dashed", colour = "grey50") +
  geom_text_repel(
    data = . %>% filter(sig) %>%
      group_by(age_group) %>%
      slice_max(abs(estimate), n = 3),
    aes(label = gene_intron),
    size = 3, max.overlaps = Inf
  ) +
  scale_colour_manual(values = effect_colors) +
  facet_wrap(~ age_group, ncol = 3) +
  labs(
    title    = "Exercise effect on splicing efficiency across age groups",
    subtitle = "Positive = improved SE after exercise",
    x        = "Effect size (PostExc - PreExc SE)",
    y        = expression(-log[10](FDR)),
    colour   = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle   = element_text(hjust = 0.5),
    strip.text      = element_text(face = "bold"),
    legend.position = "top"
  )

print(plot_volcano_exercise)
# ggsave("Figures/volcano_exercise_age_groups.svg", plot_volcano_exercise,
#         width = 14, height = 10)


# --- 12e. Heatmap — exercise effect across age groups ---
sig_exercise_targets <- zi_time_effects_fdr %>%
  filter(sig) %>%
  pull(target) %>%
  unique()

plot_heatmap <- zi_time_effects_fdr %>%
  filter(target %in% sig_exercise_targets) %>%
  ggplot(aes(x = age_group, y = gene_intron, fill = estimate)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low      = colors[1],
    mid      = "white",
    high     = colors[6],
    midpoint = 0
  ) +
  labs(
    title = "Exercise effect on SE across age groups",
    x     = NULL,
    y     = NULL,
    fill  = "Effect on SE\n(Post - Pre)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(face = "bold")
  )

print(plot_heatmap)
# ggsave("Figures/heatmap_exercise_age.svg", plot_heatmap,
#         width = 10, height = 12)


# --- 12f. Trajectory plots for top aging-associated introns ---
top_aging_targets <- zi_age_slopes_fdr %>%
  filter(sig, time == "PostExc") %>%
  slice_max(abs(estimate), n = 6) %>%
  pull(target)

plot_traj_aging <- zi_predictions %>%
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
    x        = "Age (years)",
    y        = "Splicing efficiency",
    colour   = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    strip.text      = element_text(face = "bold", size = 8),
    legend.position = "top"
  )

print(plot_traj_aging)
# ggsave("Figures/traj_top_aging.svg", plot_traj_aging,
#         width = 12, height = 8)


# --- 12g. Overlap between PreExc and PostExc age-significant introns ---
sig_preexc <- zi_age_slopes_fdr %>%
  filter(time == "PreExc", sig) %>%
  pull(target) %>% unique()

sig_postexc <- zi_age_slopes_fdr %>%
  filter(time == "PostExc", sig) %>%
  pull(target) %>% unique()

plot_venn <- ggVennDiagram(
  list(PreExc = sig_preexc, PostExc = sig_postexc)
) +
  labs(title = "Overlap of age-associated introns at rest vs post-exercise") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(plot_venn)
# ggsave("Figures/venn_preexc_postexc.svg", plot_venn,
#         width = 8, height = 6)


# =============================================================================
# SECTION 13: VISUALISATIONS — ReLiEf SENSITIVITY ANALYSIS
# =============================================================================

# --- 13a. Global trajectory — ReLiEf only ---
plot_relief_trajectory <- zi_relief_predictions %>%
  mutate(
    SE       = 1 - estimate,
    CI_low   = 1 - conf.high,
    CI_high  = 1 - conf.low,
    real_age = scaled_age * (age_max_relief - age_min_relief) + age_min_relief
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
    title    = "ReLiEf sensitivity: splicing efficiency trajectory across age",
    subtitle = "Within-study age effect free from study-age confounding",
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

print(plot_relief_trajectory)
# ggsave("Figures/relief_global_trajectory.svg", plot_relief_trajectory,
#         width = 10, height = 6)


# --- 13b. Volcano — ReLiEf age slopes ---
plot_relief_volcano_age <- zi_relief_slopes_fdr %>%
  mutate(neg_log10_fdr = -log10(adj.p)) %>%
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
    title    = "ReLiEf sensitivity: age effect on splicing efficiency",
    subtitle = "Within-study age effect",
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

print(plot_relief_volcano_age)
# ggsave("Figures/relief_volcano_age.svg", plot_relief_volcano_age,
#         width = 12, height = 6)


# --- 13c. Compare main vs ReLiEf significant introns ---
sig_main   <- zi_age_slopes_fdr %>% filter(sig) %>% pull(target) %>% unique()
sig_relief <- zi_relief_slopes_fdr %>% filter(sig) %>% pull(target) %>% unique()

plot_venn_sensitivity <- ggVennDiagram(
  list(`All studies` = sig_main, `ReLiEf only` = sig_relief)
) +
  labs(title = "Overlap: main analysis vs ReLiEf sensitivity") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(plot_venn_sensitivity)
# ggsave("Figures/venn_main_vs_relief.svg", plot_venn_sensitivity,
#         width = 8, height = 6)


# --- 13d. Side by side global trajectories: main vs ReLiEf ---
traj_main <- zi_predictions %>%
  mutate(
    SE       = 1 - estimate,
    real_age = scaled_age * (age_max - age_min) + age_min,
    source   = "All studies"
  ) %>%
  group_by(real_age, time, source) %>%
  summarise(mean_SE = mean(SE), .groups = "drop")

traj_relief <- zi_relief_predictions %>%
  mutate(
    SE       = 1 - estimate,
    real_age = scaled_age * (age_max_relief - age_min_relief) + age_min_relief,
    source   = "ReLiEf only"
  ) %>%
  group_by(real_age, time, source) %>%
  summarise(mean_SE = mean(SE), .groups = "drop")

plot_comparison <- bind_rows(traj_main, traj_relief) %>%
  ggplot(aes(x = real_age, y = mean_SE, colour = time, linetype = source)) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = c("PreExc" = colors[5],
                                 "PostExc" = colors[1])) +
  labs(
    title    = "Main analysis vs ReLiEf sensitivity: global SE trajectory",
    x        = "Age (years)",
    y        = "Mean splicing efficiency",
    colour   = NULL,
    linetype = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top"
  )

print(plot_comparison)
# ggsave("Figures/comparison_main_vs_relief.svg", plot_comparison,
#         width = 10, height = 6)


# =============================================================================
# SECTION 14: COMBINED FIGURE FOR PUBLICATION
# =============================================================================

combined_figure <- (plot_global_trajectory | plot_venn) /
  (plot_volcano_age) /
  (plot_relief_trajectory | plot_venn_sensitivity) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

print(combined_figure)
# ggsave("Figures/combined_main_figure.svg", combined_figure,
#         width = 16, height = 18)
