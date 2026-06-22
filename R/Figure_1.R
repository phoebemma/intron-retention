# load the file with the data
source("R/Load_data_for_visualisation.R")

library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggrepel)

# plot a distribution of the participants in each study

# In the metadata, most participants contributed two samples

# subset the data to only count one participant once
meta_unique <- metadata %>%
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
all <- ggplot(meta_unique, aes(x = age, fill = sex)) +
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
    hjust = 0.6,
    vjust = 0.01,
    size = 4,
    fontface= "bold.italic"
  )+
  # facet_wrap(~ study, ncol = 2) +
  theme_minimal(base_size = 16, base_family = "Arial") +
  labs(
    title = "Age distribution of all participants",
    x = "Age",
    y = NULL
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face= "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 16, face = "bold"), 
    axis.text.y = element_text(size = 16, face= "bold"),
    axis.text.x = element_text(size = 16, face= "bold"),
    axis.title = element_text(size = 16, face= "bold")
  )


# ggsave("SVG_files/all_participants.svg", plot = all, width = 14, height = 12)



by_study <- ggplot(meta_unique, aes(x = age, fill = sex)) +
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
  facet_wrap(~ study, ncol = 3, scales = "free") +
  geom_text(
    data = counts_df, aes(x = min(meta_unique$age), y= 0.5 * max_count,
                          label = paste0("n = ", n,
                                         "\nMales = ", male,
                                         "\nFemales = ", female)),
    inherit.aes = F,
    hjust = 0.1,
    vjust = 1.5,
    size = 4,
    fontface= "bold.italic"
  )+
  theme_minimal(base_size = 14, base_family = "Arial") +
  labs(
    title = " Distribution of Participants across all studies",
    x = "Age",
    y = "Number of participants"
  ) +
  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face= "bold"),
     legend.title = element_blank(),
     legend.text = element_text(size = 12, face = "bold"), 
   # legend.position = "none",
    axis.text.y = element_text(size = 16, face= "bold"),
    axis.text.x = element_text(size = 16, face= "bold"),
    axis.title = element_text(size = 16, face= "bold")
  )


ggsave("SVG_files/participants_by_study.svg", plot = by_study, width = 14, height = 12)

# Load the bubble image
# bubble_img <- png::readPNG("Figures/Bubble_chart.PNG")
# image_grob <- grid::rasterGrob(bubble_img, interpolate = TRUE)
# 
# image_plot <- wrap_elements(image_grob)

# combine with the histograms from the participants' displa


# 
# Figure_1 <- by_study + all / image_plot 
# 
# 
# Figure_1 +
#   plot_annotation(tag_levels = "A") &
#   #plot_layout(widths = c(1.1, 1))
#   theme(
#     plot.tag = element_text(size = 14, face = "bold")
#     ,
#     plot.tag.position = c(0.08, 0.98)
#   ) 

# ggsave("Figures/Trainome_Figure_2.png", bg = colors[4], width = 15, height = 10, dpi = 400)
