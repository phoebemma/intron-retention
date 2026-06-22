# This presents the predictions using marginal effects on the binomial and betabinomial models
library(glmmTMB)
library(seqwrap)
library(splines)
library(marginaleffects)
library(dplyr)
library(ggplot2)
library(broom.mixed)
library(purrr)
library(tidyr)


# Load the zero-inflated betabinomial model

zi_results <- readRDS("data/splined_zi_results.RDS")

valid_models <- compact(zi_results@models)

# this has type = "response"
 # checks the net effect of aging and exercise on retention
df_pred_resp <- map_dfr(
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
                           time = c("PreExc", "PostExc")) # Forces predictions for BOTH states
      )
    }, error = function(e) {
      # If a model fails the matrix math, return NULL and continue the loop
      return(NULL) 
    })
  },
  .id = "target"
)

saveRDS(df_pred_resp, "data/splined_pred_resp.RDS")

# This calculates the probability of perfect splicing
# Answers the question " Does aging/ exercise change whether an intron is perfectly spliced or not"
df_pred_zero_inf <- map_dfr(
  valid_models,
  function(mod) {
    tryCatch({
      # Force standard error calculation using vcov = TRUE
      predictions(
        mod, 
        type = "zprob", 
        re_formula = NA,
        vcov = TRUE, 
        newdata = datagrid(scaled_age = seq(0, 1, length.out = 100),
                           time = c("PreExc", "PostExc")) # Forces predictions for BOTH states
      )
    }, error = function(e) {
      # If a model fails the matrix math, return NULL and continue the loop
      return(NULL) 
    })
  },
  .id = "target"
)

saveRDS(df_pred_zero_inf, "data/splined_pred_zero_inf.RDS")

# This is sensitive to small changes. 
# Answeres the question "Among introns that are already partially retained, does aging/exercise make retention worse or better?"
# It is sensitive to small changes in partially retained introns


df_pred_cond <- map_dfr(
  valid_models,
  function(mod) {
    tryCatch({
      # Force standard error calculation using vcov = TRUE
      predictions(
        mod, 
        type = "conditional", 
        re_formula = NA,
        vcov = TRUE, 
        newdata = datagrid(scaled_age = seq(0, 1, length.out = 100),
                           time = c("PreExc", "PostExc")) # Forces predictions for BOTH states
      )
    }, error = function(e) {
      # If a model fails the matrix math, return NULL and continue the loop
      return(NULL) 
    })
  },
  .id = "target"
)

saveRDS(df_pred_cond, "data/splined_pred_cond.RDS")




df_age_slopes <- map_dfr(
  valid_models,
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

saveRDS(df_age_slopes, "data/splined_age_slopes.RDS")


# The interaction model 


zi_interaction_results <- readRDS("data/splined_zi_interaction_results.RDS")

valid_zi_interaction <- compact(zi_interaction_results@models)

#  Full predicted trajectories (PreExc vs PostExc across age)
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

#  Exercise effect at each age anchor a key output of interaction model
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


