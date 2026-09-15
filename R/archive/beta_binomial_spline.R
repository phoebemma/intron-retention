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



# The dataset for the contiouns SE values
cont_mat <- all_splice_reordered
cont_mat[cont_mat == 1] <- 0.999




cont_container <- seqwrap_compose(
  data       = cont_mat,
  metadata   = metadata,
  samplename = "seq_sample_id",
  modelfun   = glmmTMB,
  arguments  = alist(
    formula = y ~ splines::ns(scaled_age, df = 4) + time + sex +
      (1 | study) + (1 | participant),
    family = glmmTMB::beta_family(link = "logit")
  )
)



cont_results <- seqwrap(
  cont_container,
  return_models = TRUE,
  cores = 10,
  # subset = 1:20
)

# saveRDS(cont_results, "data/splined_beta_binomial_model.RDS")


# Filter out the NULL models while preserving the names of valid ones
valid_models <- compact(cont_results@models)




df_predictions <- map_dfr(
  valid_models,
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


# saveRDS(df_predictions, "data/splined_beta_binomial_predictions.RDS")

#  The Slope of Age Conditional on Time
df_age_slopes <- map_dfr(
  valid_models,
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


# saveRDS(df_age_slopes, "data/splined_beta_binomial_age_slopes.RDS")

# The Marginal Effect of Time Across Age
df_time_effects <- map_dfr(
  valid_models,
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



# saveRDS(df_time_effects, "data/splined_beta_binomial_time_effect.RDS")
# Extract the Significant Outcomes

# 1. Targets where the Age Slope is Significant (Pre or Post):

significant_age_targets <- df_age_slopes |>
  group_by(time) |> 
  mutate(p.adj = p.adjust(p.value, method = "BH")) |> 
  ungroup() |> 
  filter(p.adj < 0.05) # |> 
#  select(target, time, estimate, p.value, p.adj)



significant_time_targets <- df_time_effects |>
  group_by(scaled_age) |> 
  mutate(p.adj = p.adjust(p.value, method = "BH")) |> 
  ungroup() |> 
  filter(p.adj < 0.05) # |> 
 # select(target, scaled_age, contrast, estimate, p.value, p.adj)


library(ggplot2)

ggplot(df_predictions, aes(x = scaled_age, y = estimate)) +
  # Shaded confidence ribbons
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = target), alpha = 0.15) +
  # Smooth prediction lines
  geom_line(aes(color = target), linewidth = 1) +
  labs(
    title = "Model-Predicted Trajectories Across Age",
    x = "Scaled Age (0 to 1)",
    y = "Predicted Response Value",
    fill = "Model Target",
    color = "Model Target"
  ) +
  theme_minimal()
