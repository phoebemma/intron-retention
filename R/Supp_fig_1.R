source("R/Load_data_for_visualisation.R")

library(gridExtra)
library(grid)


# PLot the mean splicing efficiency across all 
long_splice <- all_splice_df %>% 
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID))

long_df <- long_splice %>%
  inner_join(metadata, by = "seq_sample_id") %>%
  group_by(participant,  study, age) %>%
  summarise(mean_SE = mean(SE), .groups = "drop") %>%
  mutate(study = recode(study,
                 "ReLiEf" = "RELIEF",
                 "copd" = "COPD",
                 "ct" = "ContraTRAIN",
                 "vol" = "VOLUME"))
 
  



    
# Distribution of SE values across particpants in the data
supp_figure1 <- ggplot(long_df, aes(x = age, y = mean_SE, color = study)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "loess", se = FALSE) +
  labs(title = "Splicing efficiency across partiicpant ages and cohorts",
    x = "Age (years)",
       y = "Mean SE per sample", 
       color = "Study") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    # plot.subtitle = element_text(hjust = 0.5, size = 13),
    strip.text = element_text(face= "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12), 
    # legend.position = ,
    axis.text.y = element_text(size = 16, face= "bold"),
    axis.text.x = element_text(size = 16, face= "bold"),
    axis.title = element_text(size = 16, face= "bold"),
    panel.grid = element_blank(), 
  )


#ggsave("SVG_files/Supplementary_figure_1.svg", plot = supp_figure1, width = 14, height = 12)






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
splice_reordered <- Relief_full_splice[, c("transcript_ID",Relief_full_meta$seq_sample_id)] 


# plot the table showing completely retained introns across all samples

# extract transcript IDs that were zeros all through
zero_introns <- all_splice_df %>%
  filter(if_all(-transcript_ID, ~ . == 0)) %>%
  dplyr::select(transcript_ID) %>%
  separate(transcript_ID, into = c("transcript_ID", "intron_ID", "chr"),
           sep = "_") %>%
  inner_join(gene_annotation,
             by = c("transcript_ID" = "ensembl_transcript_id_version"))  %>%
  mutate(
    gene_label  = ifelse(
      is.na(external_gene_name) | external_gene_name == "",
      ensembl_gene_id, external_gene_name
    ),
    gene_intron = paste(gene_label, intron_ID, sep = " : ")
  )


# extract for the Relief study
zero_introns_relief <- splice_reordered%>%
  filter(if_all(-transcript_ID, ~ . == 0)) %>%
  dplyr::select(transcript_ID) %>%
  separate(transcript_ID, into = c("transcript_ID", "intron_ID", "chr"),
           sep = "_") %>%
  inner_join(gene_annotation,
             by = c("transcript_ID" = "ensembl_transcript_id_version")) %>%
  mutate(
    gene_label  = ifelse(
      is.na(external_gene_name) | external_gene_name == "",
      ensembl_gene_id, external_gene_name
    ),
    gene_intron = paste(gene_label, intron_ID, sep = " : ")
  )


retained <- zero_introns_relief %>%
  dplyr::select(gene_label, intron_ID, gene_intron, transcript_biotype)


table_gropb <- tableGrob(retained, rows = NULL)

# Plot the introns with age related splicing efficiency at baseline

relief_contrasts_df <- relief_contrasts %>%
  filter(component != "Overall retention")

baseline_age_effect <- relief_contrasts_df %>%
  filter(hypothesis == "Age effect at baseline") %>%
  filter(sig)

table(baseline_age_effect$effect)


# SUPPLEMENATARY FIGURE 1A
base_chart <- ggplot(baseline_age_effect , aes(x = estimate, y = reorder(gene_intron, estimate), color = effect)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  #facet_wrap(~ hypothesis, scales = "free") +
  scale_color_manual(values = c("Improved SE" = colors[6],
                                "Reduced SE" = colors[1]),
                     name = "Effect") +
  labs(
    x = "Effect size",
    y = NULL,
    title = "Introns with age-related splicing efficiency at baseline"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
# plot.subtitle = element_text(hjust = 0.5, size = 13),
    strip.text = element_text(face= "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12, face = "bold"), 
    # legend.position = "none",
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 12, face= "bold"),
    axis.title = element_text(size = 12, face= "bold")
  )






# the global trajectories of the ds introns as found in the secondary analysis


# combine with the supp fig1

supp_figure1 +( wrap_elements(table_gropb) / base_chart + plot_layout(heights = c(0.25, 0.75))) 
  plot_annotation(tag_levels = "A")

# ggsave("Figures/supp_Figure_1_.png", bg = colors[4], width = 25, height = 15, dpi = 400)
  
  
# supplemmentary Figure 2
  
  aged_trajectory  <- zi_predictions %>%
    filter(target %in% baseline_age_effect$target) %>% # extract top 9
    mutate(
      SE       = 1 - estimate,
      CI_low   = 1 - conf.high,
      CI_high  = 1 - conf.low,
      real_age = scaled_age * (age_max - age_min) + age_min
    ) %>%
    left_join(
      (zi_age_slopes_fdr  %>%
         dplyr::select(target, gene_intron) %>% distinct()),
      by = "target"
    ) %>%
    filter(time == "PreExc") %>%
    ggplot(aes(x = real_age, y = SE, colour =  gene_intron)) +
    geom_ribbon(aes(ymin = CI_low, ymax = CI_high),
                alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ gene_intron, scales = "free_y", ncol = 3) +
    labs(
      title    = "Age-related trajectory of the ds introns reflected in secondary analysis",
      x        = "Age (years)",
      y        = "Splicing efficiency",
      colour   = NULL, fill = NULL
    ) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle   = element_text(hjust = 0.5),
      strip.text      = element_text(face = "bold", size = 16),
      legend.text = element_text(size = 12), 
      # legend.position = ,
      axis.text.y = element_text(size = 16, face= "bold"),
      axis.text.x = element_text(size = 16, face= "bold"),
      axis.title = element_text(size = 16, face= "bold")
    )
  
  # ggsave("Figures/supp_Figure_2.png", bg = colors[4], width = 20, height = 15, dpi = 400)
  # 
  # ggsave("SVG_files/Supplementary_figure_2.svg", plot = aged_trajectory, width = 14, height = 12)
  