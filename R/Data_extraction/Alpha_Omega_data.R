library(AOData)
library(dplyr)

# Load the participant details
ids <- idkeys %>%
  select(participant, treat, age,  sex )

Sequenced_samples <- seq_samples %>%
  select(participant, condition, time) %>%
  inner_join(ids, by = "participant") %>%
  mutate(age_group = ifelse(age <=40,  "Young","Old"),
         time = case_when(time == "T1" ~ "PreExc",
                          time == "T2" ~ "MidExc",
                          time == "T4" ~ "PostExc"))

Sequenced_samples$study <- "Alpha/Omega"
unique(Sequenced_samples$time)

pre <- Sequenced_samples %>%
  filter(time == "PreExc")

saveRDS(pre, "data/preexercise_data/Alpha_Omega_PreExc_metadata.RDS")

range(Sequenced_samples$age)

hist(Sequenced_samples$age)
colnames(ids)
colnames(Sequenced_samples)  
unique(Sequenced_samples$condition)
length(unique(Sequenced_samples$participant))

unique(Sequenced_samples$time)
unique(Sequenced_samples$sex)