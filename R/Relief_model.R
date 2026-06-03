library(dplyr)
library(tidyverse)
library(seqwrap)



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
all_splice_reordered <- Relief_full_splice[, c("transcript_ID",Relief_full_meta$seq_sample_id)] 

# Check if everything matches except the transcript_id
match(colnames(all_splice_reordered), Relief_full_meta$seq_sample_id)


# visualise the data

# Build binomial model
# This model investigates the question, "given an intron,
# what is the probability of perfect splicing as a function of age and resistance exercise training"

# derive a matrix that indicates 0 if SE is not 1
one_inflated_mat <- all_splice_reordered

one_inflated_mat[-1] <- lapply(
  one_inflated_mat[-1],
  function(x) as.integer(x == 1)
)




# Intialise argument
args_binom <- list( formula = y ~ age_group + time + sex +
                      (1 | participant), family  = binomial)

# containerise using seqwrap_compose
binom <- seqwrap_compose(data       = one_inflated_mat,
                         metadata   = Relief_full_meta,
                         samplename = "seq_sample_id",
                         modelfun   = glmmTMB::glmmTMB,
                         arguments  = args_binom)

# build model
binom_results <- seqwrap(binom,
                         return_models = FALSE,
                         cores = 10)

saveRDS(binom_results, "data/Relief_binom_model.RDS")




# The second model
# This model accepts as input the full spectrum of SE values. 
# It investigates the impact of resistance training and aging 
# on the slightest SE variations of introns.

# convert the 1.0 to 0.999. This is becasue beta-model accepts only values between 0 and one
all_splice_reordered[all_splice_reordered == 1 ] <- 0.999



# initialise the argument. This time we check the interaction of age and time
args_full <-list(formula = y ~ age_group + time + sex +
                   (1 | participant), 
                 family = glmmTMB::beta_family(link = "logit"))




# check the functions and datasets
container <- seqwrap_compose(data = all_splice_reordered,
                             metadata = Relief_full_meta,
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

saveRDS(full_model, "data/Relief_full_model.RDS")


