library(dplyr)
library(tidyverse)
library(seqwrap)
library(effects)
library(scales)
library(glmmTMB)
#source("R/Trainome_functions.R")



all_splice_df <- readRDS("data/all_splice.RDS") %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) 



all_full_metadata <- readRDS("data/all_full_metadata.RDS")
colnames(all_full_metadata)


# REORDER THE SEQUENCE ID TO MATCH BOTH DATAFRAMMES
all_splice_reordered <- all_splice_df[, c("transcript_ID",all_full_metadata$seq_sample_id)] 

# Check if everything matches except the transcript_id
match(colnames(all_splice_reordered), all_full_metadata$seq_sample_id)



# First, a one inflated model
# This models the probability of SE being 1

## Binary outcome: full responder vs not
# It answers, given an intron, what is the probability of perfect splicing as a function of age and exercise

# derive a matrix that indicates 0 if SE is not 1
one_inflated_mat <- all_splice_reordered

one_inflated_mat[-1] <- lapply(
  one_inflated_mat[-1],
  function(x) as.integer(x == 1)
)

# Intialise argument
args_binom <- list( formula = y ~ scaled_age + time + sex +
                      (1 | study) + (1 | participant), family  = binomial)

# containerise
binom <- seqwrap_compose(data       = one_inflated_mat,
                        metadata   = all_full_metadata,
                        samplename = "seq_sample_id",
                        modelfun   = glmmTMB::glmmTMB,
                        arguments  = args_binom)

binom_results <- seqwrap(binom,
                          return_models = FALSE,
                          cores = 10)

saveRDS(binom_results, "data/binom_model.RDS")

binom_sum <- seqwrap_summarise(binom_results)


binom_sum$summaries %>% 
  dplyr::select(-group)%>%
  dplyr::filter(term != "(Intercept)") %>%
  drop_na() %>%
  group_by(term) %>%                                  
  mutate(adj.p = p.adjust(p.value, method = "fdr")) %>% 
  ungroup() %>%
  filter(adj.p <= 0.05)
  

# To use the estimate from the first model as priors for the second, extract parameters
params <- binom_sum$summaries %>% 
  summarise(.by = term, 
            m = mean(estimate), 
            s = sd(estimate))  %>% 
  filter(term %in% c("scaled_age", "timePostExc", "sexmale") )

# Extract the random effects distribution
random_sd_estimate <- binom_sum$summaries %>%
        filter(term == "sd__(Intercept)") %>%
        select(target, term, estimate) %>%
        pull(estimate)

mean_sd <- mean(random_sd_estimate)
var_sd <- var(random_sd_estimate)


# Combine parameters in a data frame...
Priors_df <- #bind_rows(
  data.frame(prior = paste0("normal(",
                            round(params$m, ),
                            ",",
                            round(params$s,4 ), 
                            ")"),
             class = rep("fixef", 3),
             coef = params$term)
# ,
#                data.frame(prior = paste0(
#                        "gamma(",
#                        round(mean_sd ,4),
#                        ",",
#                        2,
#                        ")"),
#                        class = "ranef",
#                        coef = "id")
#               )  %>%
# print()



Priors_list <- list()
for( j in 1:nrow( one_inflated_mat )) {
  Priors_list[[j]] <- Priors_df
}


# containerise

informed_binom <- seqwrap_compose(data = one_inflated_mat,
                          metadata = all_full_metadata,
                          samplename = "seq_sample_id",
                          modelfun = glmmTMB::glmmTMB,
                          arguments = alist(
                            formula= y ~ scaled_age + time + sex +
                              (1 | study) + (1 | participant),
                            family = binomial,
                            priors = data.frame(
                              prior = prior,
                              class = class,
                              coef = coef)), 
                          targetdata = Priors_list)


informed_binom_results <- seqwrap(informed_binom,
                          return_models = FALSE,
                          #  subset = 1:2000,
                          cores = 10)


saveRDS(informed_binom_results, "data/informed_binom_model.RDS")

informed_binom_sum <- seqwrap_summarise(informed_binom_results)


informed_binom_sum$summaries %>% 
  dplyr::select(-group)%>%
  dplyr::filter(term != "(Intercept)") %>%
  drop_na() %>%
  group_by(term) %>%                                  
  mutate(adj.p = p.adjust(p.value, method = "fdr")) %>% 
  ungroup() %>%
  filter(adj.p <= 0.05)

informed_outputs <- informed_binom_sum$summaries %>% 
  dplyr::select(-group)%>%
  dplyr::filter(term != "(Intercept)") %>%
  drop_na() %>%
  group_by(term) %>%                                  
  mutate(adj.p = p.adjust(p.value, method = "fdr")) %>% 
  ungroup()


binom_outputs <- binom_sum$summaries %>% 
  dplyr::select(-group)%>%
  dplyr::filter(term != "(Intercept)") %>%
  drop_na() %>%
  group_by(term) %>%                                  
  mutate(adj.p = p.adjust(p.value, method = "fdr")) %>% 
  ungroup()

x <- informed_outputs %>%
  filter(adj.p <= 0.05)






# A model to check even the slightest variations in SE

# convert the 1.0 to 0.999. This is becasue beta-model accepts only values between 0 and one
all_splice_reordered[all_splice_reordered == 1 ] <- 0.999



# initialise the argument. This time we check the interaction of age and time
args_full <-list(formula = y ~  scaled_age + time + sex + (1|study) +(1|participant), 
                 family = glmmTMB::beta_family(link = "logit"))




# check the functions and datasets
container <- seqwrap_compose(data = all_splice_reordered,
                             metadata = all_full_metadata,
                             samplename = "seq_sample_id",
                             modelfun = glmmTMB::glmmTMB,
                             arguments = args_full)


# build model
full_model <- seqwrap(container,
                 # summary_fun = sum_with_pred,
                 #eval_fun = eval_mod,
                 return_models = F,
                 # subset = 1:150,
                 cores = 10)


# full_model<- readRDS("data/full_model.RDS")

 saveRDS(full_model, "data/full_model.RDS")

full_model_sum <- seqwrap_summarise(full_model)


full_model_sum$summaries %>% 
  dplyr::select(-group) %>%
  dplyr::filter(term != "(Intercept)") %>%
  drop_na() %>%
  group_by(term) %>%                                  
  mutate(adj.p = p.adjust(p.value, method = "fdr")) %>% 
  ungroup() %>%
  filter(adj.p <= 0.05)


# To use the estimate from the first model as priors for the second, extract parameters
full_params <- full_model_sum$summaries %>% 
  summarise(.by = term, 
            m = mean(estimate), 
            s = sd(estimate))  %>% 
  filter(term %in% c("scaled_age", "timePostExc", "sexmale", "scaled_age:timePostExc") ) %>%
  mutate(
    m_adj = m * 0.5,
    s_adj = pmax(ifelse(grepl(":", term), s * 3, s * 2), 0.5)
  )





# Combine parameters in a data frame...
full_Priors_df <-  data.frame(prior = paste0("normal(",
                            round(full_params$m_adj, 2),
                            ",",
                            round(full_params$s_adj,4 ), 
                            ")"),
             class = rep("beta", 4),
             coef = full_params$term) 




full_Priors_list <- list()
for( j in 1:nrow( all_splice_reordered )) {
  full_Priors_list[[j]] <- full_Priors_df
}


# containerise

informed_full_model <- seqwrap_compose(data = all_splice_reordered,
                                  metadata = all_full_metadata,
                                  samplename = "seq_sample_id",
                                  modelfun = glmmTMB::glmmTMB,
                                  arguments = alist(
                                    formula= y ~ scaled_age * time + sex +
                                      (1 | study) + (1 | participant),
                                    family = beta_family(link = "logit"),
                                    priors = data.frame(
                                      prior = prior,
                                      class = class,
                                      coef = coef)), 
                                  targetdata = full_Priors_list)


informed_full_model_results <- seqwrap(informed_full_model,
                                  return_models = FALSE,
                                  #  subset = 1:150,
                                  cores = 10)


saveRDS(informed_full_model_results, "data/informed_full_model.RDS")


