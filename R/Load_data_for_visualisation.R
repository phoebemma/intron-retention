library(dplyr)
library(tidyverse)
library(seqwrap)

# define color scale for images
## Color scale ##
colors <-  c("#d7191c", "#fdae61", "#ffffbf", "#abd9e9", "#2c7bb6", "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")

all_splice_df <- readRDS("data/all_splice.RDS") %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) 



all_full_metadata <- readRDS("data/all_full_metadata.RDS")




# Load the binary model
binom_model <- readRDS("data/binom_model.RDS") %>%
  seqwrap_summarise()


# Load the beta-binomial model and extract its summary
beta_binom_model<- readRDS("data/full_model.RDS") %>%
  seqwrap_summarise()


