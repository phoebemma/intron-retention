library(glmmTMB)
library(seqwrap)
library(splines)
library(marginaleffects)
library(dplyr)
library(ggplot2)
library(broom.mixed)
library(purrr)
# This analysis intends to capture the non-linear realtionships 
# between splicing efficiency , aging and resistance training

# load splicing efficiency data
all_splice_df <- readRDS("data/Trainome_all_splice_df.RDS")
# Load metadata
metadata <- readRDS("data/Trainome_metadata.RDS")



# REORDER THE SEQUENCE ID TO MATCH BOTH DATAFRAMMES
all_splice_reordered <- all_splice_df[, c("transcript_ID",metadata$seq_sample_id)] 

# Check if everything matches except the transcript_id
match(colnames(all_splice_reordered), metadata$seq_sample_id)


# Data preparation

# First the binary model

perfect_mat <- all_splice_reordered

perfect_mat[-1] <- lapply(
  perfect_mat[-1],
  function(x) as.integer(x == 1)
)



binom_container <- seqwrap_compose(
  data       = perfect_mat,
  metadata   = metadata,
  samplename = "seq_sample_id",
  modelfun   = glmmTMB::glmmTMB,
  arguments  =  alist(
    formula = y ~ splines::ns(scaled_age, df = 4) + time + sex +
      (1 | study) + (1 | participant),
    family  = binomial(link = "logit")
  )
)


binom_results <- seqwrap(
  binom_container,
  return_models = TRUE,
  cores = 10,
  # subset = 1:200
)

saveRDS(binom_results, "data/splined_binom_results.RDS")

# 1. Filter out the NULL models while preserving the names of valid ones
valid_models1 <- compact(binom_results@models) 
# Alternative in base R: valid_models <- Filter(Negate(is.null), binom_results@models)

# 2. Run  predictions only on the valid models
# df_predictions <- map_dfr(
#   valid_models,
#   ~ predictions(
#     .x, 
#     type = "response", 
#     re_formula = NA,
#     vcov = FALSE, # Skips standard error math to prevent the Jacobian crash
#     newdata = datagrid(scaled_age = seq(0, 1, length.out = 100))
#   ), 
#   .id = "target"
# )


df_predictions1 <- map_dfr(
  valid_models1,
  function(mod) {
    tryCatch({
      # Force standard error calculation using vcov = TRUE
      predictions(
        mod, 
        type = "response", 
        re_formula = NA,
        vcov = TRUE, 
        newdata = datagrid(scaled_age = seq(0, 1, length.out = 100),
                           time = c("PreExc", "PostExc"))
      )
    }, error = function(e) {
      # If a model fails the matrix math, return NULL and continue the loop
      return(NULL) 
    })
  },
  .id = "target"
)


saveRDS(df_predictions1, "data/splined_binomial_predictions.RDS")



#  The Slope of Age Conditional on Time
df_age_slopes1 <- map_dfr(
  valid_models1,
  function(mod) {
    tryCatch({
      # Calculate age slopes explicitly at both time points
      avg_slopes(
        mod, 
        variables = "scaled_age", 
        by = "time",
        type = "response", 
        re_formula = NA
      )
    }, error = function(e) return(NULL))
  },
  .id = "target"
)

saveRDS(df_age_slopes1, "data/splined_binomial_age_slopes.RDS")

# The Marginal Effect of Time Across Age
df_time_effects1 <- map_dfr(
  valid_models1,
  function(mod) {
    tryCatch({
      # Calculate the discrete change of time (post vs pre) across age steps
      avg_comparisons(
        mod,
        variables = "time",
        by = "scaled_age",
        # Test at specific age anchors: 0, 0.25, 0.5, 0.75, and 1.0
        newdata = datagrid(scaled_age = seq(0, 1, by = 0.25)),
        type = "response",
        re_formula = NA
      )
    }, error = function(e) return(NULL))
  },
  .id = "target"
)

saveRDS(df_time_effects1, "data/splined_binomial_time_effect.RDS")

# Extract the Significant Outcomes

# 1. Targets where the Age Slope is Significant (Pre or Post):

significant_age_targets <- df_age_slopes |>
  group_by(time) |> 
  mutate(p.adj = p.adjust(p.value, method = "BH")) |> 
  ungroup() |> 
  filter(p.adj < 0.05) |> 
  select(target, time, estimate, p.value, p.adj)
# Interpretation: Look at the time column. If a target is significant in post 
# but not pre, training fundamentally altered how age relates to the outcome.




# Targets where the Training Effect (time) is Significant:

significant_time_targets <- df_time_effects |>
  group_by(scaled_age) |> 
  mutate(p.adj = p.adjust(p.value, method = "BH")) |> 
  ungroup() |> 
  filter(p.adj < 0.05) |> 
  select(target, scaled_age, contrast, estimate, p.value, p.adj)
# Interpretation: Look at the scaled_age column. The estimate tells you the 
# exact probability jump caused by training at that specific age checkpoint




# THE ULTIMATE INTERACTION TEST (the difference in slopes)
# To find out targets where theinteraction is significant, that is the slope of age
# in post is significatly different from the slope of age in pre

df_interaction_test <- map_dfr(
  valid_models,
  function(mod) {
    tryCatch({
      # Generate the slopes for both times
      s <- avg_slopes(mod, variables = "scaled_age", by = "time", type = "response", re_formula = NA)
      # Test if the difference between the two slopes equals zero
      hypotheses(s, hypothesis = "b1 - b2 = 0")
    }, error = function(e) return(NULL))
  },
  .id = "target"
)

# Filter for targets with a true, shifting interaction
true_interactions <- df_interaction_test |> 
  mutate(p.adj = p.adjust(p.value, method = "BH")) |> 
  filter(p.adj < 0.05)


# identify targets significant for both factors


# 1. Get unique targets where age has a significant effect (pre or post)
targets_aging <- significant_age_targets |> 
  pull(target) |> 
  unique()

# 2. Get unique targets where training (time) has a significant effect
targets_training <- significant_time_targets |> 
  pull(target) |> 
  unique()

# 3. Find the overlap (targets significant for BOTH aging AND training)
dual_significant_targets <- intersect(targets_aging, targets_training)

print(paste("Found", length(dual_significant_targets), "targets significant for both factors."))





# regenerate of filter predictions with time split

# Generate dense prediction curves specifically for the targets that matter
df_plot_data <- map_dfr(
  valid_models[dual_significant_targets], # Only loops through the dual-significant subset
  ~ predictions(
    .x, 
    type = "response", 
    re_formula = NA,
    # Generate age sequence across BOTH time points simultaneously
    newdata = datagrid(
      scaled_age = seq(0, 1, length.out = 100),
      time = c("pre", "post")
    )
  ), 
  .id = "target"
)



# Create the interaction plot
ggplot(df_plot_data, aes(x = scaled_age, y = estimate, color = time, fill = time)) +
  # Main trend lines for pre and post training
  geom_line(linewidth = 1.2) +
  # 95% Confidence interval ribbons to show uncertainty
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, color = NA) +
  # Separate panel for each dual-significant target
  facet_wrap(~ target, scales = "free_y") + 
  # Aesthetics and styling
  scale_color_manual(values = c("pre" = "#E69F00", "post" = "#56B4E9")) +
  scale_fill_manual(values = c("pre" = "#E69F00", "post" = "#56B4E9")) +
  labs(
    x = "Scaled Age (0 to 1)",
    y = "Predicted Probability (y)",
    title = "Targets with Significant Main/Interaction Effects",
    subtitle = "Showing altered aging trajectories pre vs. post training",
    color = "Training Status",
    fill = "Training Status"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"), # Bold panel titles (target names)
    legend.position = "bottom"
  )

