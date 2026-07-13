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
    hjust = 0.001,
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
    title = " Distribution of Participants across all 5 studies",
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
    axis.title = element_text(size = 16, face= "bold"),
   legend.position = "none"
  )


# ggsave("SVG_files/participants_by_study.svg", plot = by_study, width = 14, height = 12)








# Distribution of splicing efficiency values across the dataset
long_splice <- all_splice_df %>% 
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID))

dist_df <- long_splice %>%
  filter(!is.na(SE)) %>%
  mutate(
    SE_range = case_when(
      SE <= 0.25 ~ "SE from 0.00-0.25",
      SE <= 0.49 ~ "SE from 0.26-0.49",
      SE <= 0.75 ~ "SE from 0.50-0.75",
      TRUE       ~ "SE from 0.76-1.00"
    ),
    SE_range = factor(SE_range,
                      levels = c("SE from 0.00-0.25", "SE from 0.26-0.49",
                                 "SE from 0.50-0.75", "SE from 0.76-1.00"))
  ) %>%
  count(SE, SE_range) %>%
  mutate(percent = n / sum(n) * 100,
         label = ifelse(SE %in% c(0, 1),
                        paste0(round(percent, 1), "%"),
                        NA)) 


introns_distribution <-   ggplot(dist_df,
                                 aes(x = factor(SE, levels = sort(unique(SE))),
                                     y = percent,
                                     fill = case_when(
                                       SE == 1 ~ "one",
                                       SE == 0 ~ "zero",
                                       TRUE    ~ "other"
                                     ))) +
  geom_col() +
  
  # geom_text(
  #   data = subset(dist_df, SE %in% c(0, 1)),
  #   aes(label = paste0(round(percent, 1), "%")),
  #   vjust = 0.5,
  #   hjust=0.5,
  #   size = 5,
  #   fontface = "bold"
  # ) +
  # Label for SE = 0
  geom_text(
    data = subset(dist_df, SE == 0),
    aes(label = paste0(round(percent, 1), "%")),
    vjust = 1.5,
    hjust = -0.2   # adjust differently
  ) +
  
  # Label for SE = 1
  geom_text(
    data = subset(dist_df, SE == 1),
    aes(label = paste0(round(percent, 1), "%")),
    vjust = 1.5,
    hjust = 1.2    # adjust differently
  )+

  facet_wrap(~SE_range, scales = "free") +
  
  scale_fill_manual(
    values = c(
      "zero"  = colors[1],  
      "one"   = colors[6],   
      "other" = colors[11]
    ),
    guide = "none"
  ) +
  
  scale_x_discrete(breaks = seq(0, 1, by = 0.1)) +
  theme_minimal(base_size = 16, base_family = "Arial") +
  labs(
    x = "Splicing Efficiency (SE)",
    y = "Percentage",
    title = "Distribution of Splicing Efficiency Values across the pooled dataset"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_blank(), # removes facet titles
    legend.title = element_blank(),
    legend.text = element_text(size = 16, face = "bold"), 
    axis.text.y = element_text(size = 16, face= "bold"),
   # axis.text.y = element_blank(),  remove y axis title
    axis.text.x = element_text(angle = 90, size = 16, face= "bold"),
    axis.title = element_text(size = 16, face= "bold"),
   axis.title.y = element_blank(),
    panel.spacing.y = unit(2, "lines"),
    panel.grid = element_blank(), # removes grid lines 
    
    panel.border = element_rect( # add a subtle border
      colour = "grey80",   
      fill = NA,
      linewidth = 0.5
    ))


# ggsave("SVG_files/distribution_SE.svg", plot = introns_distribution, width = 16, height = 15)






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




long_relief <- splice_reordered %>% 
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID))

dist_df_relief <- long_relief %>%
  filter(!is.na(SE)) %>%
  mutate(
    SE_range = case_when(
      SE <= 0.25 ~ "SE from 0.00-0.25",
      SE <= 0.49 ~ "SE from 0.26-0.49",
      SE <= 0.75 ~ "SE from 0.50-0.75",
      TRUE       ~ "SE from 0.76-1.00"
    ),
    SE_range = factor(SE_range,
                      levels = c("SE from 0.00-0.25", "SE from 0.26-0.49",
                                 "SE from 0.50-0.75", "SE from 0.76-1.00"))
  ) %>%
  count(SE, SE_range) %>%
  mutate(percent = n / sum(n) * 100,
         label = ifelse(SE %in% c(0, 1),
                        paste0(round(percent, 1), "%"),
                        NA)) 


relief_distribution <-   ggplot(dist_df_relief,
                                 aes(x = factor(SE, levels = sort(unique(SE))),
                                     y = percent,
                                     fill = case_when(
                                       SE == 1 ~ "one",
                                       SE == 0 ~ "zero",
                                       TRUE    ~ "other"
                                     ))) +
  geom_col() +
  
  # geom_text(
  #   data = subset(dist_df_relief, SE %in% c(0, 1)),
  #   aes(label = paste0(round(percent, 1), "%")),
  #   vjust = 1.5,
  #   hjust=0.5,
  #   size = 4,
  #   fontface = "bold"
  # ) +
  # Label for SE = 0
  geom_text(
    data = subset(dist_df_relief, SE == 0),
    aes(label = paste0(round(percent, 1), "%")),
    vjust = 1.5,
    hjust = -0.2   # adjust differently
  ) +
  
  # Label for SE = 1
  geom_text(
    data = subset(dist_df_relief, SE == 1),
    aes(label = paste0(round(percent, 1), "%")),
    vjust = 1.5,
    hjust = 1.2    # adjust differently
  ) +

  
  facet_wrap(~SE_range, scales = "free") +
  
  scale_fill_manual(
    values = c(
      "zero"  = colors[1],  
      "one"   = colors[6],   
      "other" = colors[11]
    ),
    guide = "none"
  ) +
  
  scale_x_discrete(breaks = seq(0, 1, by = 0.1)) +
  theme_minimal(base_size = 16, base_family = "Arial") +
  labs(
    x = "Splicing Efficiency (SE)",
    y = "Percentage",
    title = "Distribution of Splicing Efficiency values across the RELIEF cohort"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_blank(), # removes facet titles
    legend.title = element_blank(),
    legend.text = element_text(size = 16, face = "bold"), 
    axis.text.y = element_text(size = 16, face= "bold"),
    axis.text.x = element_text(angle = 90, size = 16, face= "bold"),
    axis.title = element_text(size = 16, face= "bold"),
    panel.spacing.y = unit(2, "lines"),
    panel.grid = element_blank(), # removes grid lines 

 panel.border = element_rect( # add a subtle border
   colour = "grey80",
   fill = NA,
   linewidth = 0.5
   )
)
print(relief_distribution)

ggsave("SVG_files/distribution_Relief_SE.svg", plot = relief_distribution, width = 16, height = 15)





# Load the bubble image
bubble_img <- png::readPNG("Figures/Bubble_chart.PNG")
image_grob <- grid::rasterGrob(bubble_img, interpolate = TRUE)

image_plot <- wrap_elements(image_grob)




# explore the introns completely retained, or perfectly spliced across all samples


long_splice <- all_splice_df %>% 
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID))

dist_df <- long_splice %>%
  filter(!is.na(SE)) %>%
  mutate(
    SE_range = case_when(
      SE <= 0.25 ~ "SE from 0.00-0.25",
      SE <= 0.49 ~ "SE from 0.26-0.49",
      SE <= 0.75 ~ "SE from 0.50-0.75",
      TRUE       ~ "SE from 0.76-1.00"
    ),
    SE_range = factor(SE_range,
                      levels = c("SE from 0.00-0.25", "SE from 0.26-0.49",
                                 "SE from 0.50-0.75", "SE from 0.76-1.00"))
  ) %>%
  count(SE, SE_range) %>%
  mutate(percent = n / sum(n) * 100,
         label = ifelse(SE %in% c(0, 1),
                        paste0(round(percent, 1), "%"),
                        NA)) 

all_zeros <- long_splice %>%
  filter(SE == 0) %>%
  pivot_wider(names_from = seq_sample_id,
              values_from = SE)

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


one_introns <- all_splice_df %>%
  filter(if_all(-transcript_ID, ~ . == 1)) %>%
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




ego_ones <- enrichGO(gene =  unique(one_introns$external_gene_name),
                     keyType = "SYMBOL",
                     universe = gene_exp_df$gene_name,
                     OrgDb = org.Hs.eg.db, 
                     ont = "BP", 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05, 
                     readable = T)


## Output results from GO analysis to a table
cluster_ones <- data.frame(ego_ones)

go_ones <- dotplot(ego_ones,
                   showCategory = 6,
                   font.size = 8, title = "Enriched biological processes in genes containing introns perfectly spliced acrosss all samples") +
  theme(axis.text = element_text(size = 16), axis.text.y = element_text(size = 16), axis.title.x = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.text = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"))

# ggsave("SVG_files/gene_ontology_ones.svg", plot = go_ones, width = 16, height = 15)

# Repeat for RELIEF

one_introns_relief <- splice_reordered %>%
  filter(if_all(-transcript_ID, ~ . == 1)) %>%
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



ego_ones_relief <- enrichGO(gene =  unique(one_introns_relief$external_gene_name),
                            keyType = "SYMBOL",
                            universe = gene_exp_df$gene_name,
                            OrgDb = org.Hs.eg.db, 
                            ont = "BP", 
                            pAdjustMethod = "BH", 
                            qvalueCutoff = 0.05, 
                            readable = T)


## Output results from GO analysis to a table
cluster_ones_relief <- data.frame(ego_ones_relief)

go_ones_relief <- dotplot(ego_ones_relief,
                          showCategory = 6,
                          font.size = 16, title = "Enriched biological processes in genes containing introns perfectly spliced in Relief samples") +
  theme(axis.text = element_text(size = 16), axis.text.y = element_text(size = 16), axis.title.x = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.text = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold") )

#ggsave("SVG_files/gene_ontology_ones_RELIEF.svg", plot = go_ones_relief, width = 16, height = 15)


# 
 Figure_1 <- ( by_study + all + image_plot) /(relief_distribution + introns_distribution) /
   (go_ones_relief + go_ones)
# 
# 
Figure_1 +
  plot_annotation(tag_levels = "A") &
  #plot_layout(widths = c(1.1, 1))
  theme(
    plot.tag = element_text(size = 14, face = "bold")
    ,
    plot.tag.position = c(0.08, 0.98)
  )

 ggsave("Figures/Trainome_Figure_1_.png", bg = colors[4], width = 40, height = 25, dpi = 400)
