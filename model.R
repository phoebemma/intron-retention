library(dplyr)
library(tidyverse)
library(seqwrap)
library(effects)
library(scales)
library(glmmTMB)
library(DHARMa)

source("R/Trainome_functions.R")

sigma_summary2 <- function(x) {
  
  ## ---- Convergence information ----
  conv <- if (is.null(x$fit$convergence)) 1 else x$fit$convergence[[1]]
  
  ## ---- DHARMa residual diagnostics ----
  # Simulated scaled residuals for uniformity and dispersion tests
  residObj <- DHARMa::simulateResiduals(x, plot = FALSE)
  
  unif <- DHARMa::testUniformity(residObj, plot = FALSE)
  disp <- DHARMa::testDispersion(residObj, plot = FALSE)
  
  ## ---- Random effects (observation-level RE, if present) ----
  tidy_coef <- broom.mixed::tidy(x)
  
  if (any(tidy_coef$group %in% "seq_sample_id")) {
    olre.sd <- tidy_coef |>
      dplyr::filter(group == "seq_sample_id") |>
      dplyr::pull(estimate)
  } else {
    olre.sd <- NA
  }
  
  ## ---- Extract beta precision (phi) ----
  sdr_tab <- as.data.frame(summary(x$sdr))
  
  phi_est <- sdr_tab["betadisp", "Estimate"]
  phi_se  <- sdr_tab["betadisp", "Std. Error"]
  
  ## ---- Extract marginal mean (accounts for zero inflation) ----
  # Conditional mean on link scale
  mu_link <- predict(x, type = "link")
  
  # Zero-inflation probability on link scale
  zi_link <- predict(x, type = "zprob")
  
  # Convert to response scale
  mu  <- plogis(mu_link)
  zi  <- plogis(zi_link)
  
  # Marginal mean E[Y] = (1 - zi) * mu
  mu_marginal <- (1 - zi) * mu
  
  ## ---- Combine outputs ----
  out <- data.frame(
    phi         = phi_est,
    phi_se      = phi_se,
    olre.sd     = olre.sd,
    log_mu      = mean(log(mu_marginal)),
    zi_prob     = mean(zi),
    convergence = conv,
    pdHess      = x$sdr$pdHess,
    aic         = AIC(x),
    unif.p      = unif$p.value,
    unif.stat   = unif$statistic,
    disp.p      = disp$p.value,
    disp.stat   = disp$statistic
  )
  
  return(out)
}

all_splice_df <- readRDS("data/all_splice.RDS") %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) 



all_full_metadata <- readRDS("data/all_full_metadata.RDS")
colnames(all_full_metadata)


# REORDER THE SEQUENCE ID TO MATCH BOTH DATAFRAMMES
all_splice_reordered <- all_splice_df[, c("transcript_ID",all_full_metadata$seq_sample_id)] 

# Check if everything matches except the transcript_id
match(colnames(all_splice_reordered), all_full_metadata$seq_sample_id)




# Convert intron matrix to long format


# splice_long <- all_splice_reordered |>
#   pivot_longer(
#     cols = where(is.numeric),
#     names_to  = "seq_sample_id",
#     values_to = "y"
#   )
# 
# hist(splice_long$y)
# 
# 
# # define full response to be where SE is = 1
# splice_long <- splice_long |>
#   mutate(full_response = as.integer(y == 1))
# 
# 
# 
# splice_long <- splice_long |>
#   left_join(
#     all_full_metadata,
#     by = "seq_sample_id"
#   )
# derive a matrix that indicates 0 if SE is not 1
one_inflated_mat <- all_splice_reordered

one_inflated_mat[-1] <- lapply(
  one_inflated_mat[-1],
  function(x) as.integer(x == 1)
)

# One inflated model
# This models the probability of SE being 1

## Binary outcome: full responder vs not
# It answers, given an intron, what is the probability of perfect splicing as a function of age and exercise

args_partA <- list(
  formula = y ~ scaled_age + time + sex +
    (1 | study) + (1 | participant),
  family  = binomial()
)

mA <- seqwrap_compose(
  data       = one_inflated_mat,
  metadata   = all_full_metadata,
  samplename = "seq_sample_id",
  modelfun   = glmmTMB::glmmTMB,
  arguments  = args_partA
)

mA_results <- seqwrap(
  mA,
  return_models = FALSE,
  cores = 6
)

mA_sums <- seqwrap_summarise(mA_results)


mod_sum <- mA_sums$summaries


mod_eval <- mA_sums$evaluations





# 
model1_pred <- mod_sum %>%
  dplyr::select(-group) %>%
  dplyr::filter(term != "(Intercept)") %>%
  drop_na() %>%
  group_by(term) %>%                                  
  mutate(adj.p = p.adjust(p.value, method = "fdr")) %>% 
  ungroup() %>%                                       
  mutate(odds_ratio = exp(estimate)) %>%
  filter(adj.p<= 0.05)







# convert the 1.0 to 0.999. This is becasue beta-model accepts only values between 0 and one
all_splice_reordered[all_splice_reordered == 1 ] <- 0.999



# initialise the argument
args_full <-list(formula = y ~  scaled_age * time + sex + (1|study) +(1|participant), 
                 family = glmmTMB::beta_family(link = "logit"))




# The old model


# check the functions and datasets
container <- seqwrap_compose(data = all_splice_reordered,
                             metadata = all_full_metadata,
                             samplename = "seq_sample_id",
                             modelfun = glmmTMB::glmmTMB,
                             arguments = args_full,
                            # summary_fun = sum_with_pred,
                             eval_fun = sigma_summary2 
                            )


# build model
model <- seqwrap(container,
                # summary_fun = sum_with_pred,
                 #eval_fun = eval_mod,
                 return_models = F,
                 # subset = 1:15,
                 cores = ncores-2)






mod_results <- seqwrap_summarise(model)


mod_sumB <- mod_results$summaries


mod_evalB <- mod_results$evaluations





# 
modelB_pred <- mod_sumB %>%
  dplyr::select(-group) %>%
  dplyr::filter(term != "(Intercept)") %>%
  drop_na() %>%
  group_by(term) %>%                                  
  mutate(adj.p = p.adjust(p.value, method = "fdr")) %>% 
  ungroup() %>%                                       
  mutate(odds_ratio = exp(estimate)) %>%
  filter(adj.p<= 0.05)


age_sig <- modelB_pred %>%
  filter(term == "scaled_age")

# Part B- beta model for sub-maxila responses
# Prepare submaximal data
# Exclude those with perfect SE scores

# beta_mat <- all_splice_reordered
# 
# beta_mat[-1] <- lapply(
#   beta_mat[-1],
#   function(x) ifelse(x == 1, NA_real_, x)
# )
# 
# beta_mat <- beta_mat |>
#   tidyr::drop_na()


# Uncomment this if you want to retain introns where over 30% are perfectly spliced
# keep <- beta_mat[-1] |>
#   summarise(across(everything(), ~ mean(. == 1, na.rm = TRUE))) |>
#   pivot_longer(everything()) |>
#   pull(value) < 0.3
# 
# beta_mat <- beta_mat[keep, ]
# Part B would check the magnitude of sub-optimal SE scores

# args_partB <- list(
#   formula = y ~ scaled_age * time + sex +
#     (1 | study) + (1 | participant),
#   family  = glmmTMB::beta_family(link = "logit")
# )
# 
# mB <- seqwrap_compose(
#   data       = beta_mat,
#   metadata   = all_full_metadata,
#   samplename = "seq_sample_id",
#   modelfun   = glmmTMB::glmmTMB,
#   eval_fun   = sigma_summary2,
#   arguments  = args_full
# )
# 
# mB_results <- seqwrap(
#   mB,
#   return_models = FALSE,
#   cores = 6
# )

#mB_sums <- seqwrap_summarise(mB_results)




targets <- mod_results$evaluations |>
  filter(!is.na(phi_se)) |>      # precision identified
  distinct(target) |>
  pull(target)


# Estimate mean on response scale
log_mu_df <- all_splice_reordered |>
  mutate(log_mu = log(rowMeans(across(where(is.numeric))))) |>
  select(target = transcript_ID, log_mu)

phi_df <- mod_results$evaluations |>
  filter(target %in% targets) |>
  select(target, phi, phi_se) |>
  left_join(log_mu_df, by = "target")

# LOESS trend for beta precision
trend_model <- loess(
  phi ~ log_mu,
  data = phi_df,
  weights = 1 / (phi_se^2),
  span = 0.7
)
# Construct intron-specific EB penalties for phi

phi_prior_df <- log_mu_df |>
  mutate(
    pred_phi = pmax(predict(trend_model, newdata = log_mu_df), 2),
    prior = paste0("normal(", round(pred_phi, 2), ",1)"))

# EB regularisation for fixed and random effects

fixef_dist <- mod_results$summaries |>
  filter(term %in% c("scaled_age", "time", "sex")) |>
  summarise(.by = term,
            m = mean(estimate),
            s = sd(estimate))

fixef_priors <- data.frame(
  prior = paste0("normal(", round(fixef_dist$m, 2), ",",
                 round(fixef_dist$s, 2), ")"),
  class = "fixef",
  coef  = fixef_dist$term
)

ranef_sd <- mod_results$summaries |>
  filter(term == "sd__(Intercept)") |>
  pull(estimate)

ranef_prior <- data.frame(
  prior = paste0("gamma(", round(mean(ranef_sd), 2), ",2)"),
  class = "ranef",
  coef  = c("study", "participant")
)


# Assemble intron-specific EB penalty tests
Priors_list <- lapply(seq_len(nrow(log_mu_df)), function(j) {
  
  intron <- log_mu_df$target[j]
  
  bind_rows(
    fixef_priors,
    ranef_prior,
    data.frame(
      prior = phi_prior_df |>
        filter(target == intron) |>
        pull(prior),
      class = "fixef_phi",   # precision term
      coef  = "1"
    )
  )
})


# Informed model
mB_EB <- seqwrap_compose(
  data       = all_splice_reordered,
  metadata   = all_full_metadata,
  samplename = "seq_sample_id",
  modelfun   = glmmTMB::glmmTMB,
  eval_fun   = sigma_summary2,
  targetdata = Priors_list,   # <<< THIS enables EB
  arguments  = alist(
    formula = y ~ scaled_age * time + sex +
      (1 | study) + (1 | participant),
    family  = glmmTMB::beta_family(link = "logit"),
    priors  = data.frame(
      prior = prior,
      class = class,
      coef  = coef
    )
  )
)

mB_EB_results <- seqwrap(
  mB_EB,
  return_models = FALSE,
  cores = 6
)

mB_EB_sums <- seqwrap_summarise(mB_EB_results)



sum_mB    <- mB_sums$summaries
sum_mBEB  <- mB_EB_sums$summaries

eval_mB   <- mB_sums$evaluations
eval_mBEB <- mB_EB_sums$evaluations


# precision and shrinkage
phi_compare <- eval_mB |>
  select(target, phi, phi_se) |>
  rename(phi_mle = phi) |>
  inner_join(
    eval_mBEB |>
      select(target, phi) |>
      rename(phi_eb = phi),
    by = "target"
  )



# Plot: MLE ?? vs EB ??

ggplot(phi_compare, aes(x = phi_mle, y = phi_eb)) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Precision (??) - unregularised",
    y = "Precision (??) - EB regularised",
    title = "Empirical Bayes shrinkage of beta precision (??)"
  )


# Mean precision relationship
phi_trend <- phi_compare |>
  inner_join(
    eval_mB |>
      select(target, log_mu),
    by = "target"
  )

ggplot(phi_trend, aes(x = log_mu)) +
  geom_point(aes(y = phi_mle), alpha = 0.3, colour = "grey50") +
  geom_point(aes(y = phi_eb), alpha = 0.3, colour = "steelblue") +
  labs(
    y = "Precision (??)",
    x = "log(mean sub-maximal splicing efficiency)",
    title = "Mean-precision relationship before and after EB"
  )


# Fixed effect concordance (age by time)
fixef_compare <- sum_mB |>
  filter(term == "scaled_age:time") |>
  select(target, estimate) |>
  rename(est_mle = estimate) |>
  inner_join(
    sum_mBEB |>
      filter(term == "scaled_age:time") |>
      select(target, estimate) |>
      rename(est_eb = estimate),
    by = "target"
  )

ggplot(fixef_compare, aes(est_mle, est_eb)) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    x = "Age × time effect (MLE)",
    y = "Age × time effect (EB)",
    title = "Stability of age-dependent training effects under EB"
  )


# Distribution of standard errors
se_compare <- sum_mB |>
  filter(term == "scaled_age:time") |>
  select(target, std.error) |>
  rename(se_mle = std.error) |>
  inner_join(
    sum_mBEB |>
      filter(term == "scaled_age:time") |>
      select(target, std.error) |>
      rename(se_eb = std.error),
    by = "target"
  )

ggplot(se_compare, aes(x = se_mle, y = se_eb)) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    x = "SE (MLE)",
    y = "SE (EB)",
    title = "Reduction in uncertainty under EB regularisation"
  )


