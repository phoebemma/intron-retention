library(dplyr)
library(tidyverse)
library(seqwrap)
library(glmmTMB)
library(ggplot2)

# Load the dataframe containing the 6 different RT intervention studies
# Check the /Data_extraction subfolder for how each individual study's data extraction was done

# Check "data_compilation.R for how data was compiled



all_splice_df <- readRDS("data/all_splice.RDS") %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) 



all_full_metadata <- readRDS("data/all_full_metadata.RDS")
colnames(all_full_metadata)


# REORDER THE SEQUENCE ID TO MATCH BOTH DATAFRAMMES
all_splice_reordered <- all_splice_df[, c("transcript_ID",all_full_metadata$seq_sample_id)] 

# Check if everything matches except the transcript_id
match(colnames(all_splice_reordered), all_full_metadata$seq_sample_id)


# visualise the data

## Color scale ##
colors <-  c("#d7191c", "#fdae61", "#ffffbf", "#abd9e9", "#2c7bb6", "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")


# plot a distribution of the participants in each study

# In the metadata, most participants contributed two samples

# subset the data to only count one participant once
meta_unique <- all_full_metadata %>%
  distinct(participant, .keep_all = TRUE)%>%
  mutate(study = recode(study,
                        "ReLiEf" = "RELIEF",
                        "copd" = "COPD",
                        "ct" = "ContraTRAIN",
                        "vol" = "VOLUME"))

# create counts per study
# Aim is to add number of participants  in image

counts_df <- meta_unique %>%
  group_by(study) %>%
  summarise(
    n = n(),
    male = sum(sex == "male"),
    female = sum(sex == "female"),
    .groups = "drop"
  )


# for the summarised image

sum_df <- meta_unique %>%
  summarise(
    n = n(),
    male = sum(sex == "male"),
    female = sum(sex == "female")
  )


# To help set in-image text
max_count <- ggplot_build(
  ggplot(meta_unique, aes(x = age)) +
    geom_histogram(binwidth = 5)
)$data[[1]]$count %>% max()




# Plot the distribution of all participants in one image
ggplot(meta_unique, aes(x = age, fill = sex)) +
  geom_histogram(position = position_dodge(width = 5),
                 alpha = 0.5,
                 binwidth = 5,
                 color = "black") +
  scale_fill_manual(
    values = c(
      "male" = colors[7],
      "female" = colors[3]
    )
  ) +
  geom_text(
    data = sum_df, aes(x = min(meta_unique$age), y= 0.5 * max_count,
                       label = paste0("n = ", n,
                                      "\nMales = ", male,
                                      "\nFemales = ", female)),
    inherit.aes = F,
    hjust = 0.5,
    vjust = 0.5,
    size = 4,
    fontface= "bold.italic"
  )+
  # facet_wrap(~ study, ncol = 2) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Age distribution of all participants",
    x = "Age",
    y = "Number of participants"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank()
  )



# ggsave("Figures/Dist_all_participants.png", bg = colors[4], scale = 2.5, dpi = 400)



ggplot(meta_unique, aes(x = age, fill = sex)) +
  geom_histogram(position = "dodge",
                 alpha = 0.5,
                 binwidth = 5,
                 color = "black") +
  scale_fill_manual(
    values = c(
      "male" = colors[7],
      "female" = colors[3]
    )
  ) +
  facet_wrap(~ study, ncol = 2, scales = "fixed") +
  geom_text(
    data = counts_df, aes(x = min(meta_unique$age), y= 0.5 * max_count,
                          label = paste0("n = ", n,
                                         "\nMales = ", male,
                                         "\nFemales = ", female)),
    inherit.aes = F,
    hjust = 0.1,
    vjust = 1.5,
    size = 4,
    # fontface= "bold.italic"
  )+
  theme_minimal(base_size = 12) +
  labs(
    title = " Distribution of Participants across all studies",
    x = "Age",
    y = "Count"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank()
  )

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
args_binom <- list( formula = y ~ scaled_age + time + sex +
                      (1 | study) + (1 | participant), family  = binomial)

# containerise using seqwrap_compose
binom <- seqwrap_compose(data       = one_inflated_mat,
                         metadata   = all_full_metadata,
                         samplename = "seq_sample_id",
                         modelfun   = glmmTMB::glmmTMB,
                         arguments  = args_binom)

# build model
binom_results <- seqwrap(binom,
                         return_models = FALSE,
                         cores = 10)

#saveRDS(binom_results, "data/binom_model.RDS")




# The second model
# This model accepts as input the full spectrum of SE values. 
# It investigates the impact of resistance training and aging 
# on the slightest SE variations of introns.

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



# saveRDS(full_model, "data/full_model.RDS")
