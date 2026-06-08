# load the file with the data
source("R/Load_data_for_visualisation.R")

library(ggplot2)
library(patchwork)


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
all <- ggplot(meta_unique, aes(x = age, fill = sex)) +
  geom_histogram(position = position_dodge(width = 5),
                 alpha = 0.7,
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
    title = "Age and gender distribution of all participants",
    x = "Age",
    y = NULL
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.title = element_blank(),
    legend.text = element_text(size = 12, face = "bold"), 
    
    axis.title.x = element_text(size = 14, face = "bold"),
  #  axis.title.y = element_text(size = 14, face = "bold"),
    
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
  plot.background  = element_rect(fill = "white", colour = NA),
  panel.background = element_rect(fill = "white", colour = NA)
  
    
  )



# Visualise the participants by study
by_study <- ggplot(meta_unique, aes(x = age, fill = sex)) +
  geom_histogram(position = "dodge",
                 alpha = 0.7,
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
     fontface= "bold.italic"
  )+
  theme_minimal(base_size = 12) +
  labs(
    title = " Distribution of Participants across all studies",
    x = "Age",
    y = "Number of participants"
  ) +
theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "None", 
    
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    strip.text = element_text(face = "bold", size = 14),
    strip.background = element_rect(fill = "grey90", colour = NA),
    
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
    
    )





Fig1_plot <-  by_study + all
Fig1_plot +
  plot_annotation(tag_levels = "A") &
  #plot_layout(widths = c(1.1, 1))
  theme(
    plot.tag = element_text(size = 14, face = "bold")
    ,
    plot.tag.position = c(0.08, 0.98)) 
# ggsave("Figures/Figure_1.png",  width = 20, height = 12, dpi = 400, device = ragg::agg_png)