# =============================================================================
# RECOMMENDED ANALYSIS PIPELINE
# Impact of Aging and Resistance Training on Intron Splicing Efficiency
# =============================================================================
# Models:
#   1. Binomial      — probability of perfect splicing (SE = 1 vs SE < 1)
#   2. ZI-Beta main  — degree of retention, main effects of age and time
#   3. ZI-Beta interaction — age-dependent exercise effect on retention
# =============================================================================

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


# =============================================================================
# SECTION 1: LOAD AND COMPILE METADATA
# =============================================================================

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

metadata <- bind_rows(
  copd_metadata, vol_metadata, ct_metadata, ao_metadata, relief_metadata
) %>%
  mutate(
    age          = round(age, 0),
    sex          = factor(sex, levels = c("female", "male")),
    time         = factor(time, levels = c("PreExc", "PostExc")),
    scaled_age   = round(rescale(age), digits = 2),
    participant  = paste0(study, "_", participant)
  )


# =============================================================================
# SECTION 2: LOAD AND COMPILE SPLICING DATA
# =============================================================================

all_splice_df <- readRDS("data/copd_splicing_data.RDS") %>%
  inner_join(readRDS("data/volume_splicing_data.RDS"),    by = "transcript_ID") %>%
  inner_join(readRDS("data/contratrain_splicing_data.RDS"), by = "transcript_ID") %>%
  inner_join(readRDS("data/Alpha_Omega_splicing_data.RDS"), by = "transcript_ID") %>%
  inner_join(readRDS("data/Relief_splicing_data.RDS"),    by = "transcript_ID") %>%
  drop_na()

# Keep only samples present in both splice data and metadata
intersect_ids <- intersect(colnames(all_splice_df), metadata$seq_sample_id)

all_splice_df <- all_splice_df %>%
  subset(select = c("transcript_ID", intersect_ids)) %>%
  drop_na()

# Subset metadata to only samples present in splice data
metadata <- metadata %>%
  filter(seq_sample_id %in% colnames(all_splice_df))

# Reorder splice columns to match metadata row order
all_splice_reordered <- all_splice_df[, c("transcript_ID", metadata$seq_sample_id)]


# =============================================================================
# SECTION 3: PREPARE MODEL MATRICES
# =============================================================================

# --- Binomial matrix: 1 = perfectly spliced, 0 = any retention ---
binom_mat <- all_splice_reordered
binom_mat[-1] <- lapply(binom_mat[-1], function(x) as.integer(x == 1))

# --- ZI-Beta matrix: flip SE so 0 = perfectly spliced (structural zero) ---
zi_mat <- all_splice_reordered
zi_mat[-1] <- 1 - all_splice_reordered[-1]


# =============================================================================
# SECTION 4: MODEL 1 — BINOMIAL
# Biological question: Does aging/exercise change the probability of
# perfect splicing?
# =============================================================================

binom_container <- seqwrap_compose(
  data       = binom_mat,
  metadata   = metadata,
  samplename = "seq_sample_id",
  modelfun   = glmmTMB::glmmTMB,
  arguments  = alist(
    formula = y ~ splines::ns(scaled_age, df = 4) + time + sex +
                  (1 | study) + (1 | participant),
    family  = binomial(link = "logit")
  )
)

binom_results <- seqwrap(
  binom_container,
  return_models = TRUE,
  cores         = 10
)

saveRDS(binom_results, "data/binom_results.RDS")
# binom_results <- readRDS("data/binom_results.RDS")


# =============================================================================
# SECTION 5: MODEL 2 — ZERO-INFLATED BETA (MAIN EFFECTS)
# Biological question: Does aging/exercise change the degree of intron
# retention, independent of whether splicing is perfect?
# =============================================================================

zi_main_container <- seqwrap_compose(
  data       = zi_mat,
  metadata   = metadata,
  samplename = "seq_sample_id",
  modelfun   = glmmTMB::glmmTMB,
  arguments  = alist(
    formula   = y ~ splines::ns(scaled_age, df = 4) + time + sex +
                    (1 | study) + (1 | participant),
    ziformula = ~ splines::ns(scaled_age, df = 4) + time,
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


# =============================================================================
# SECTION 6: MODEL 3 — ZERO-INFLATED BETA (AGE x TIME INTERACTION)
# Biological question: Does the splicing response to exercise differ
# depending on the age of the participant?
# =============================================================================

zi_interaction_container <- seqwrap_compose(
  data       = zi_mat,
  metadata   = metadata,
  samplename = "seq_sample_id",
  modelfun   = glmmTMB::glmmTMB,
  arguments  = alist(
    formula   = y ~ splines::ns(scaled_age, df = 4) * time + sex +
                    (1 | study) + (1 | participant),
    ziformula = ~ splines::ns(scaled_age, df = 4) * time,
    family    = glmmTMB::beta_family(link = "logit")
  )
)

zi_interaction_results <- seqwrap(
  zi_interaction_container,
  return_models = TRUE,
  cores         = 10
)

saveRDS(zi_interaction_results, "data/zi_interaction_results.RDS")
# zi_interaction_results <- readRDS("data/zi_interaction_results.RDS")


# =============================================================================
# SECTION 7: MARGINAL EFFECTS — MAIN EFFECTS MODEL
# Extract predicted retention trajectories and main effects of age and time
# =============================================================================

valid_zi_main <- compact(zi_main_results@models)

# 7a. Predicted retention across age for both time points (overall response)
zi_predictions <- map_dfr(
  valid_zi_main,
  function(mod) {
    tryCatch({
      predictions(
        mod,
        type      = "response",
        re_formula = NA,
        vcov      = TRUE,
        newdata   = datagrid(
          scaled_age = seq(0, 1, length.out = 100),
          time       = c("PreExc", "PostExc")
        )
      )
    }, error = function(e) NULL)
  },
  .id = "target"
)

# 7b. Probability of perfect splicing across age (zi component)
zi_zprob <- map_dfr(
  valid_zi_main,
  function(mod) {
    tryCatch({
      predictions(
        mod,
        type      = "zprob",
        re_formula = NA,
        vcov      = TRUE,
        newdata   = datagrid(
          scaled_age = seq(0, 1, length.out = 100),
          time       = c("PreExc", "PostExc")
        )
      )
    }, error = function(e) NULL)
  },
  .id = "target"
)

# 7c. Degree of retention among partially retained introns (conditional)
zi_conditional <- map_dfr(
  valid_zi_main,
  function(mod) {
    tryCatch({
      predictions(
        mod,
        type      = "conditional",
        re_formula = NA,
        vcov      = TRUE,
        newdata   = datagrid(
          scaled_age = seq(0, 1, length.out = 100),
          time       = c("PreExc", "PostExc")
        )
      )
    }, error = function(e) NULL)
  },
  .id = "target"
)

# 7d. Age slopes at both time points
zi_age_slopes <- map_dfr(
  valid_zi_main,
  function(mod) {
    tryCatch({
      avg_slopes(
        mod,
        variables  = "scaled_age",
        by         = "time",
        newdata    = datagrid(scaled_age = seq(0, 1, by = 0.25)),
        type       = "response",
        re_formula = NA
      )
    }, error = function(e) NULL)
  },
  .id = "target"
)

saveRDS(zi_predictions,  "data/zi_main_predictions.RDS")
saveRDS(zi_zprob,        "data/zi_main_zprob.RDS")
saveRDS(zi_conditional,  "data/zi_main_conditional.RDS")
saveRDS(zi_age_slopes,   "data/zi_main_age_slopes.RDS")


# =============================================================================
# SECTION 8: MARGINAL EFFECTS — INTERACTION MODEL
# Extract age-dependent exercise effects
# =============================================================================

valid_zi_interaction <- compact(zi_interaction_results@models)

# 8a. Full predicted trajectories (PreExc vs PostExc across age)
zi_int_predictions <- map_dfr(
  valid_zi_interaction,
  function(mod) {
    tryCatch({
      predictions(
        mod,
        type      = "response",
        re_formula = NA,
        vcov      = TRUE,
        newdata   = datagrid(
          scaled_age = seq(0, 1, length.out = 100),
          time       = c("PreExc", "PostExc")
        )
      )
    }, error = function(e) NULL)
  },
  .id = "target"
)

# 8b. Exercise effect at each age anchor — key output of interaction model
zi_time_effects <- map_dfr(
  valid_zi_interaction,
  function(mod) {
    tryCatch({
      avg_comparisons(
        mod,
        variables  = "time",
        by         = "scaled_age",
        newdata    = datagrid(scaled_age = seq(0, 1, by = 0.25)),
        type       = "response",
        re_formula = NA
      )
    }, error = function(e) NULL)
  },
  .id = "target"
)

saveRDS(zi_int_predictions, "data/zi_interaction_predictions.RDS")
saveRDS(zi_time_effects,    "data/zi_time_effects.RDS")


# =============================================================================
# SECTION 9: FDR FILTERING
# =============================================================================

# --- From main effects model summaries ---
zi_main_summary <- seqwrap_summarise(zi_main_results)

zi_main_outputs <- zi_main_summary$summaries %>%
  dplyr::select(-group) %>%
  filter(term != "(Intercept)", term != "sexmale") %>%
  drop_na() %>%
  group_by(term) %>%
  mutate(
    adj.p = p.adjust(p.value, method = "fdr"),
    term  = recode(term,
                   "scaled_age"   = "Aging",
                   "timePostExc"  = "Resistance Training")
  ) %>%
  ungroup() %>%
  mutate(
    sig           = adj.p <= 0.05,
    neg_log10_fdr = -log10(adj.p),
    effect        = case_when(
      estimate > 0 & adj.p <= 0.05 ~ "Increased Retention",
      estimate < 0 & adj.p <= 0.05 ~ "Decreased Retention",
      TRUE                          ~ "No effect"
    )
  )

# --- From interaction model: age-dependent exercise effect ---
sig_interaction_introns <- zi_time_effects %>%
  group_by(target) %>%
  mutate(adj.p = p.adjust(p.value, method = "fdr")) %>%
  filter(adj.p <= 0.05) %>%
  pull(target) %>%
  unique()

# Introns where exercise effect changes direction across age
age_dependent_introns <- zi_time_effects %>%
  group_by(target) %>%
  mutate(adj.p = p.adjust(p.value, method = "fdr")) %>%
  filter(adj.p <= 0.05) %>%
  filter(any(estimate > 0) & any(estimate < 0)) %>%
  pull(target) %>%
  unique()

# Subset interaction predictions to significant introns only
zi_int_predictions_sig <- zi_int_predictions %>%
  filter(target %in% sig_interaction_introns)


# =============================================================================
# SECTION 10: LOAD ANNOTATION FOR DOWNSTREAM VISUALISATION
# =============================================================================

gene_annotation <- readRDS("data/ensembl_gene_annotation.RDS")

intron_length <- readr::read_tsv("data_new/Alpha_Omega_SpliceQ_outputs/A_102.tsv") %>%
  distinct(across(6:ncol(.)), .keep_all = TRUE) %>%
  mutate(
    transcript_ID = paste0(transcript_ID, "_", intron_ID, "_", chr),
    intron_length = abs((sj3start - sj5end) + 1)
  ) %>%
  dplyr::select(transcript_ID, intron_length)

# Annotate main outputs
zi_main_annotated <- zi_main_outputs %>%
  inner_join(intron_length, by = c("target" = "transcript_ID")) %>%
  mutate(transcript_ID = str_split(target, "_", simplify = TRUE)[, 1]) %>%
  inner_join(gene_annotation, by = c("transcript_ID" = "ensembl_transcript_id_version")) %>%
  separate(target, into = c(NA, "intron_ID", NA), sep = "_", remove = FALSE) %>%
  mutate(
    gene_label  = ifelse(
      is.na(external_gene_name) | external_gene_name == "",
      ensembl_gene_id, external_gene_name
    ),
    gene_intron = paste(gene_label, intron_ID, sep = " : ")
  )


# =============================================================================
# SECTION 11: VISUALISATIONS
# =============================================================================

colors <- c("#d7191c", "#fdae61", "#ffffbf", "#abd9e9", "#2c7bb6",
            "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")

effect_colors <- c(
  "Increased Retention" = colors[1],
  "Decreased Retention" = colors[6],
  "No effect"           = "grey70"
)


# --- 11a. Volcano plot: Aging effect (main effects model) ---
volcano_aging <- zi_main_annotated %>%
  filter(term == "Aging") %>%
  ggplot(aes(estimate, neg_log10_fdr, colour = effect)) +
  geom_point(alpha = 0.6, size = 1.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = . %>% filter(sig) %>% slice_max(abs(estimate), n = 10),
    aes(label = gene_intron), size = 3, max.overlaps = Inf
  ) +
  scale_colour_manual(values = effect_colors) +
  labs(
    title    = "Effect of Aging on Intron Retention",
    subtitle = "ZI-Beta model — response scale",
    x        = "Effect size (log-odds retention)",
    y        = expression(-log[10](FDR))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title  = element_blank()
  )


# --- 11b. Volcano plot: Resistance Training effect ---
volcano_RT <- zi_main_annotated %>%
  filter(term == "Resistance Training") %>%
  ggplot(aes(estimate, neg_log10_fdr, colour = effect)) +
  geom_point(alpha = 0.6, size = 1.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = . %>% filter(sig) %>% slice_max(abs(estimate), n = 10),
    aes(label = gene_intron), size = 3, max.overlaps = Inf
  ) +
  scale_colour_manual(values = effect_colors) +
  labs(
    title    = "Effect of Resistance Training on Intron Retention",
    subtitle = "ZI-Beta model — response scale",
    x        = "Effect size (log-odds retention)",
    y        = expression(-log[10](FDR))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title  = element_blank()
  )

volcano_aging + volcano_RT +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

# ggsave("Figures/volcano_main_effects.png", width = 14, height = 6, dpi = 400)


# --- 11c. Age trajectory plots for top aging-associated introns ---
top_aging <- zi_main_annotated %>%
  filter(term == "Aging", sig) %>%
  slice_max(abs(estimate), n = 6) %>%
  pull(target)

traj_aging <- zi_predictions %>%
  filter(target %in% top_aging) %>%
  left_join(
    zi_main_annotated %>% dplyr::select(target, gene_intron) %>% distinct(),
    by = "target"
  )

ggplot(traj_aging, aes(x = scaled_age, y = estimate, colour = time, fill = time)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, colour = NA) +
  geom_line(size = 0.8) +
  facet_wrap(~ gene_intron, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = c("PreExc" = colors[5], "PostExc" = colors[1])) +
  scale_fill_manual(values   = c("PreExc" = colors[5], "PostExc" = colors[1])) +
  labs(
    title  = "Retention trajectories across age — top aging-associated introns",
    x      = "Scaled age",
    y      = "Predicted retention (1 - SE)",
    colour = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    strip.text  = element_text(face = "bold"),
    legend.position = "top"
  )

# ggsave("Figures/traj_aging_top6.png", width = 12, height = 8, dpi = 400)


# --- 11d. Age trajectory plots for top RT-associated introns ---
top_RT <- zi_main_annotated %>%
  filter(term == "Resistance Training", sig) %>%
  slice_max(abs(estimate), n = 6) %>%
  pull(target)

traj_RT <- zi_predictions %>%
  filter(target %in% top_RT) %>%
  left_join(
    zi_main_annotated %>% dplyr::select(target, gene_intron) %>% distinct(),
    by = "target"
  )

ggplot(traj_RT, aes(x = scaled_age, y = estimate, colour = time, fill = time)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, colour = NA) +
  geom_line(size = 0.8) +
  facet_wrap(~ gene_intron, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = c("PreExc" = colors[5], "PostExc" = colors[1])) +
  scale_fill_manual(values   = c("PreExc" = colors[5], "PostExc" = colors[1])) +
  labs(
    title  = "Retention trajectories — top RT-associated introns",
    x      = "Scaled age",
    y      = "Predicted retention (1 - SE)",
    colour = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    strip.text  = element_text(face = "bold"),
    legend.position = "top"
  )

# ggsave("Figures/traj_RT_top6.png", width = 12, height = 8, dpi = 400)


# --- 11e. Age-dependent exercise effect — interaction model ---
# Shows how the PreExc vs PostExc difference changes across age
zi_time_effects_sig <- zi_time_effects %>%
  group_by(target) %>%
  mutate(adj.p = p.adjust(p.value, method = "fdr")) %>%
  filter(target %in% age_dependent_introns) %>%
  left_join(
    zi_main_annotated %>% dplyr::select(target, gene_intron) %>% distinct(),
    by = "target"
  )

ggplot(zi_time_effects_sig, aes(x = scaled_age, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, fill = colors[5]) +
  geom_line(colour = colors[5], size = 0.8) +
  facet_wrap(~ gene_intron, scales = "free_y") +
  labs(
    title    = "Age-dependent exercise effect on intron retention",
    subtitle = "Estimate = PostExc - PreExc retention at each age point",
    x        = "Scaled age",
    y        = "Exercise effect (PostExc - PreExc)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text    = element_text(face = "bold")
  )

# ggsave("Figures/age_dependent_exercise_effect.png", width = 14, height = 10, dpi = 400)


# --- 11f. Crossing trajectories for age-dependent introns ---
# Most compelling visualisation — shows PreExc/PostExc lines diverging or crossing
traj_interaction <- zi_int_predictions_sig %>%
  left_join(
    zi_main_annotated %>% dplyr::select(target, gene_intron) %>% distinct(),
    by = "target"
  ) %>%
  filter(target %in% age_dependent_introns)

ggplot(traj_interaction, aes(x = scaled_age, y = estimate, colour = time, fill = time)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.12, colour = NA) +
  geom_line(size = 0.9) +
  facet_wrap(~ gene_intron, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = c("PreExc" = colors[5], "PostExc" = colors[1])) +
  scale_fill_manual(values   = c("PreExc" = colors[5], "PostExc" = colors[1])) +
  labs(
    title    = "Age-dependent divergence of exercise effect on intron retention",
    subtitle = "Interaction model — trajectories that cross or diverge with age",
    x        = "Scaled age",
    y        = "Predicted retention (1 - SE)",
    colour   = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle   = element_text(hjust = 0.5),
    strip.text      = element_text(face = "bold"),
    legend.position = "top"
  )

# ggsave("Figures/crossing_trajectories.png", width = 14, height = 10, dpi = 400)


# --- 11g. Summary count plot — how many introns per effect category ---
zi_main_annotated %>%
  filter(sig) %>%
  count(term, effect) %>%
  ggplot(aes(x = term, y = n, fill = effect)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.3) +
  scale_fill_manual(values = effect_colors) +
  labs(
    title = "Number of significantly affected introns",
    x     = NULL,
    y     = "Number of introns",
    fill  = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# ggsave("Figures/effect_counts.png", width = 8, height = 6, dpi = 400)
