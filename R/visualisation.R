library(ggplot2)
library(dplyr)
library(tidyverse)

## Color scale ##
colors <-  c("#d7191c", "#fdae61", "#ffffbf", "#abd9e9", "#2c7bb6", "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")




all_full_metadata <- readRDS("data/all_full_metadata.RDS")
colnames(all_full_metadata)


ggplot(all_full_metadata, aes(x = age, fill = study)) +
  geom_density(alpha = 0.4) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Age Distribution Across Studies",
    x = "Age",
    y = "Density"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank()
  )

# ggplot(all_full_metadata, aes(x = age)) +
#   geom_histogram(binwidth = 5, fill = "grey70", color = "black") +
#   facet_wrap(~ study, ncol = 2) +
#   theme_minimal(base_size = 12) +
#   labs(
#     title = "Age Distribution by Study",
#     x = "Age",
#     y = "Count"
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     strip.text = element_text(face = "bold")
#   )




meta_unique <- all_full_metadata %>%
  distinct(participant, .keep_all = TRUE)%>%
  mutate(study = recode(study,
                        "ReLiEf" = "RELIEF",
                        "copd" = "COPD",
                        "ct" = "ContraTRAIN",
                        "vol" = "VOLUME"))

# create counts per study

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


max_count <- ggplot_build(
  ggplot(meta_unique, aes(x = age)) +
    geom_histogram(binwidth = 5)
)$data[[1]]$count %>% max()


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
    title = "Age distribution ofall participants",
    x = "Age",
    y = "Number of participants"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank()
  )
