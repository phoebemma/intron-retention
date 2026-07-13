library(glmmTMB)
library(seqwrap)
library(splines)
library(marginaleffects)
library(dplyr)
library(ggplot2)

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

# The dataset for the contiouns SE values
cont_mat <- all_splice_reordered
cont_mat[cont_mat == 1] <- 0.999


# Define the arguments for the binary model

# args_binom <- list(
#   formula = y ~ splines::ns(scaled_age, df = 4) * time + sex +
#     (1 | study) + (1 | participant),
#   family  = binomial(link = "logit")
# )



# 
# 
# sum_fun <- function(model, data, target_name = NULL, ...) {
#   
#   # Ensure marginaleffects doesn't block complex hypotheses
#   options(marginaleffects_safe = FALSE)
#   # Define grids for marginal effects
#   grid <- marginaleffects::datagrid(
#     scaled_age = seq(0, 1, length.out = 50),
#     time = c("pre", "post")
#   )
#   
#   
#   # Predictions (trajectory across age)
#   pred <- tryCatch({
#     predictions(
#       model,
#       newdata = grid,
#       type = "response",
#       re.form = NA   # population-level
#     ) %>%
#       as.data.frame() %>%
#       transform(type = "prediction")
#   }, error = function(e) NULL)
#   
#   # ----- Pre vs Post contrast (training effect) -----
#   contr <- tryCatch({
#     comparisons(
#       model,
#       newdata = grid,
#       variables = "time"
#     ) %>%
#       as.data.frame() %>%
#       transform(type = "contrast")
#   }, error = function(e) NULL)
#   
#   # ----- Combine -----
#   out <- dplyr::bind_rows(pred, contr)
#   
#   # Add feature ID if available
#   if (!is.null(target_name)) {
#     out$feature <- target_name
#   }
#   
#   return(out)
# }
# 
# sum_with_mfx <- function(model, data, target_name = NULL, ...) {
#   
#   options(marginaleffects_safe = FALSE)
#   
#   # Define grid (dense, spline-friendly)
#   grid <- marginaleffects::datagrid(
#     scaled_age = seq(0, 1, length.out = 50),
#     time = levels(model.frame(model)$time),
#     model = model
#   )
#   
#   # ---- Predictions (like allEffects) ----
#   pred <- tryCatch({
#     predictions(
#       model,
#       newdata = grid,
#       type = "response",
#       re.form = NA
#     ) %>%
#       as.data.frame() %>%
#       transform(type = "prediction")
#   }, error = function(e) NULL)
#   
#   # ---- Contrasts (new, powerful addition) ----
#   contr <- tryCatch({
#     comparisons(
#       model,
#       newdata = grid,
#       variables = "time"
#     ) %>%
#       as.data.frame() %>%
#       transform(type = "contrast")
#   }, error = function(e) NULL)
#   
#   out <- dplyr::bind_rows(pred, contr)
#   
#   # Attach feature ID
#   if (!is.null(target_name)) {
#     out$feature <- target_name
#   }
#   
#   return(out)
# }
# sum_fun <- function(model, data, target = NULL, ...) {
#   
#   # Expand original data across age
#   ages <- seq(0, 1, length.out = 50)
#   newdata <- data[rep(seq_len(nrow(data)), each = length(ages)), ]
#   newdata$scaled_age <- rep(ages, times = nrow(data))
#   
#   pred <- marginaleffects::predictions(
#     model,
#     newdata = newdata,
#     type = "response"
#   )
#   
#   contr <- marginaleffects::comparisons(
#     model,
#     newdata = newdata,
#     variables = "time"
#   )
#   
#   out <- dplyr::bind_rows(
#     as.data.frame(pred),
#     as.data.frame(contr)
#   )
#   
#   if (!is.null(target)) {
#     out$feature <- target
#   }
#   
#   return(out)
# }
# 



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

binom_results@models
names(binom_results@models)

x <- binom_results@summaries


mod_sum <- bind_rows(within(binom_results@summaries, rm(names(which(binom_results@summaries == "NULL"))))) %>%
  mutate(target = rep(names(which(binom_results@summaries != "NULL")), each = 13))


mod_sum <- binom_results@summaries %>%
  discard(is.null) %>%              # remove NULL elements
  mutate(target = rep(., each = 13))  # bind rows + add name as column


mod_sum <- binom_results@summaries %>%
  discard(is.null) %>%
  imap_dfr(~ mutate(.x, target = .y))


x <- x[[1]]
str(binom_results$fits[[1]])
models <- binom_results$fits

x <- seqwrap_summarise(binom_results)
x$evaluations
x <- x$summaries
binom_results@errors$warn_eval



colnames(x)




# define arguments for beta-binomial model

# args_cont <- list(
#   formula = y ~ splines::ns(scaled_age, df = 4) * time + sex +
#     (1 | study) + (1 | participant),
#   family = beta_family(link = "logit")
# )


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

saveRDS(cont_results, "data/splined_beta_binomial_model.RDS")

# extract one model with the hope of knowing how to build the summary function
m1 <- cont_results@models[[1]]

x <- cont_results@models
names(x)
# 
# dat <- data.frame(m1)
# 
# ages <- seq(0, 1, length.out = 50)
# 
# pred <- marginaleffects::predictions(
#   m1,
#   newdata = as.data.frame(model.frame(m1)),
#   type = "response",
#   re.form = NA
# )



dat <- metadata

# keep only rows used in the model
dat <- dat[complete.cases(dat), ]

ages <- seq(0, 1, length.out = 50)

newdata <- dat[rep(seq_len(nrow(dat)), each = length(ages)), ]
newdata$scaled_age <- rep(ages, times = nrow(dat))


pred <- marginaleffects::predictions(
  m1,
  newdata = newdata,
  type = "response",
  re.form = NA
)




names(cont_results@models)


all_pred <- purrr::imap_dfr(models, function(m, id) {
  
  dat <- as.data.frame(model.frame(m))
  
  pred <- marginaleffects::predictions(
    m,
    newdata = dat,
    type = "response"
  )
  
  out <- as.data.frame(pred)
  out$feature <- id   # ??? THIS is the target_name
  
  return(out)
})
cont_results@errors$err_sum

y <- seqwrap_summarise(cont_results)
y$evaluations
y <- y$summaries
unique(y$term)







# extract results
res_binom <- binom_results$summary %>%
  mutate(model_type = "perfect_splicing")

res_cont <- cont_results$summary %>%
  mutate(model_type = "overall_efficiency")

res_all <- bind_rows(res_binom, res_cont)


# clean variables
# Convert scaled_age ??? real age
res_all <- res_all %>%
  mutate(age = 18 + scaled_age * (93 - 18))

# summarise across introns 
summary_pred <- res_all %>%
  filter(type == "prediction") %>%
  group_by(model_type, age, time) %>%
  summarise(
    mean = mean(estimate, na.rm = TRUE),
    p25 = quantile(estimate, 0.25, na.rm = TRUE),
    p75 = quantile(estimate, 0.75, na.rm = TRUE),
    .groups = "drop"
  )
# Training effect
summary_contrast <- res_all %>%
  filter(type == "contrast") %>%
  group_by(model_type, age) %>%
  summarise(
    mean = mean(estimate, na.rm = TRUE),
    p25 = quantile(estimate, 0.25, na.rm = TRUE),
    p75 = quantile(estimate, 0.75, na.rm = TRUE),
    .groups = "drop"
  )



# Visualisation
ggplot(summary_pred,
       aes(x = age, y = mean, color = time)) +
  geom_line(size = 1) +
  facet_wrap(~ model_type, scales = "free_y") +
  labs(
    x = "Age (years)",
    y = "Predicted splicing outcome",
    title = "Age-related splicing patterns"
  )

# Training versus age
ggplot(summary_contrast,
       aes(x = age, y = mean)) +
  geom_line(size = 1) +
  facet_wrap(~ model_type) +
  labs(
    x = "Age (years)",
    y = "Training effect (post vs pre)",
    title = "Exercise effect across age"
  )


# feature level effects