# This script shows the script for investigatng intron retention across age.
# First the Relief model is built. Then the pooled data is used as validation of the Relief Results

# In the ReLiEf analysis, SE is investigated between Young (18-30) vs Old (60-93) comparison
# Models fitted:
#   1. ZI-Beta ReLiEf Young vs Old, age_group * time + sex + (1|participant)
#   2. ZI-Beta pooled  age + time + sex + study + study:time + (1|participant)

library(dplyr)
library(tidyverse)
library(seqwrap)
library(reliefdata)
library(marginaleffects)




# extract the Relief SPLICE-q data and metadata



Relief_full_meta <- readRDS("data/Relief_metadata.RDS")%>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age) %>%
  mutate(sex = factor(sex, levels = c("female", "male")),
         time = factor(time, levels = c("PreExc", "PostExc")),
         age_group = case_when(age < 40 ~ "Young",
                               age > 40 ~ "Old"),
         age_group = factor(age_group, levels = c("Young", "Old")))


Relief_full_splice <- readRDS("data/Relief_splicing_data.RDS")  %>%
  drop_na()



# REORDER THE SEQUENCE ID TO MATCH BOTH DATAFRAMMES
splice_reordered <- Relief_full_splice[, c("transcript_ID",Relief_full_meta$seq_sample_id)] 

# Check if everything matches except the transcript_id
match(colnames(splice_reordered), Relief_meta$seq_sample_id)


# Build model using zero-inflated beta-binomial since 57% of data are 1s

# To achieve this, the data has to be flipped


zi_splice_df <- splice_reordered 

zi_splice_df[-1] <- 1 - splice_reordered[-1] # Flip all except first column

# write the summary function
sum_fun_relief <- function(m) {
  options(marginaleffects_safe = FALSE)
  
  grid <- marginaleffects::datagrid(
    model     = m,
    age_group = unique,
    time      = unique
  )
  
  hyp <- c(
    "young_pre"       = "b1 = 0",
    "young_post"      = "b2 = 0",
    "old_pre"         = "b3 = 0",
    "old_post"        = "b4 = 0",
    "train_young"     = "b2 - b1 = 0",
    "train_old"       = "b4 - b3 = 0",
    "age_effect_pre"  = "b3 - b1 = 0",
    "age_effect_post" = "b4 - b2 = 0",
    "train_old_young" = "(b4 - b3) - (b2 - b1) = 0"
  )
  
  # Overall predicted retention (zi + beta combined)
  response <- marginaleffects::predictions(
    m, re.form = NA, type = "response",
    newdata = grid, hypothesis = hyp
  ) |> data.frame() |> dplyr::mutate(component = "response")
  
  # Probability of perfect splicing (zi component only) 
  zprob <- marginaleffects::predictions(
    m, re.form = NA, type = "zprob",
    newdata = grid, hypothesis = hyp
  ) |> data.frame() |>
    dplyr::mutate(component = "zprob")
  
  # Degree of retention among partially retained introns (beta component only)
  conditional <- marginaleffects::predictions(
    m, re.form = NA, type = "conditional",
    newdata = grid, hypothesis = hyp
  ) |> data.frame() |> dplyr::mutate(component = "conditional")
  
  dplyr::bind_rows(response, zprob, conditional)
}

#
# Fit the model
#
zi_relief_container <- seqwrap_compose(
  data        = zi_splice_df,
  metadata    = Relief_meta,
  samplename  = "seq_sample_id",
  modelfun    = glmmTMB::glmmTMB,
  arguments   = alist(
    formula   = y ~ age_group * time + sex + (1 | participant),
    ziformula = ~ age_group * time,
    family    = glmmTMB::beta_family(link = "logit")
  ),
  summary_fun = sum_fun_relief
)

zi_relief_results <- seqwrap(
  zi_relief_container,
  return_models = FALSE,  
  cores         = 10
)



saveRDS(zi_relief_results, "data/Relief_zi_model.RDS") 




# MODEL TWO Validation of Relief results
# Key design decisions:
#   - Study included as FIXED effect to account for age-study confounding
#   - study:time interaction added to account for study-specific exercise effects
#   - Age effect estimated after removing study-specific exercise responses
#   - SE scores flipped (1 - SE) so 0 = perfect splicing (structural zero)
#   - Natural splines (df=4) for non-linear age modelling
#   - Marginalisation over study AND study:time for clean population inferencee
#   - Age increments of 0.10 for avg_slopes 
#   study:time added to both formula and ziformula
#   This removes study-specific exercise effects before estimating age slope
#   Prevents study-age confounding from masking the true age trajectory



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

# Compile metadata study as fixed factor, vol as reference level
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

# Diagnose age-study overlap important context for interpretation
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

