library(dplyr)
library(tidyverse)
library(seqwrap)
library(glmmTMB)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(ggridges)
library(UpSetR)
library(ggplot2)
library(ggplotify)
library(clusterProfiler)
library(org.Hs.eg.db)
library(magick)
library(forcats)
library(patchwork)

# Load both splicing and metadata files

all_splice_df <- readRDS("data/all_splice.RDS") %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) 



all_full_metadata <- readRDS("data/all_full_metadata.RDS")
colnames(all_full_metadata)


# REORDER THE SEQUENCE ID TO MATCH BOTH DATAFRAMMES
all_splice_reordered <- all_splice_df[, c("transcript_ID",all_full_metadata$seq_sample_id)] 


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
    y = NULL
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank()
  )



# ggsave("Figures/Dist_all_participants.png", bg = colors[4], scale = 2.5, dpi = 400)


# Visualise the participants by study
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
    y = "Number of participants"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    legend.position = "None"
  )





Fig1_plot <-  by_study + all
Fig1_plot +
  plot_annotation(tag_levels = "A") &
  #plot_layout(widths = c(1.1, 1))
theme(
  plot.tag = element_text(size = 14, face = "bold")
  ,
  plot.tag.position = c(0.08, 0.98)
) 


# ggsave("Figures/Figure_1.png", bg = colors[4], scale = 4, dpi = 400)



# Load the binary model
binom_results<- readRDS("data/binom_model.RDS")


# Load the beta-binomial model
full_model<- readRDS("data/full_model.RDS")



# extract the model summaries in the beta binomial model
informed_binom_sum <- seqwrap_summarise(binom_results)


# filter significantly differantially spliced introns
binom_model_outputs <- informed_binom_sum$summaries %>% 
  dplyr::select(-group) %>%
  dplyr::filter(term != "(Intercept)") %>%
  drop_na() %>%
  group_by(term) %>%                                  
  mutate(adj.p = p.adjust(p.value, method = "fdr")) %>% 
  ungroup() %>%
  # select significantly differentially spliced introns
  filter(adj.p <= 0.05) %>%
  mutate(term = recode(term,
                       "scaled_age" = "Aging",
                       "timePostExc" = "Resistance Training"))


# extract the results of the non-binarised model

full_model_sum <- seqwrap_summarise(full_model)

# The non-binarised model 
full_model_outputs <- full_model_sum$summaries %>% 
  dplyr::select(-group) %>%
  dplyr::filter(term != "(Intercept)") %>%
  drop_na() %>%
  group_by(term) %>%                                  
  mutate(adj.p = p.adjust(p.value, method = "fdr")) %>% 
  ungroup() %>%
  # select significantly differentially spliced introns
  filter(adj.p <= 0.05) %>%
  mutate(term = recode(term,
                       "scaled_age" = "Aging",
                       "timePostExc" = "Resistance Training"))






# BELOW IS EXPLORATION OF THE DATA OUTPUTS

# Load the gene annotation file
gene_annotation <- readRDS("data/ensembl_gene_annotation.RDS")

# Load one file from which we will extract intron length
# This is valid as only introns quantified in all samples were included in the analyses
intron_length <- readr::read_tsv("data_new/Alpha_Omega_SpliceQ_outputs/A_102.tsv") %>%
  
  distinct(across(6:ncol(.)), .keep_all = T) %>% # Removes duplicates based on columns 6 to end
  mutate(transcript_ID = paste0(transcript_ID, "_", intron_ID, "_", chr),
         intron_length = abs((sj3start - sj5end) + 1) ) %>% # Ensures positive length regardless of strand
  dplyr::select(transcript_ID, intron_length)


# based on the model outputs,
# create a column that shows if SE is improved or reduced
binom_model_outputs <- binom_model_outputs %>%
  inner_join(intron_length, by = c("target" = "transcript_ID")) %>%
  mutate(effect = case_when(estimate > 0 & adj.p <= 0.05 ~ "Improved SE", 
                            estimate < 0 & adj.p <= 0.05 ~ "Reduced SE" ,
                            estimate < 0 & adj.p > 0.05 ~ "No effect",
                            estimate > 0 & adj.p > 0.05 ~ "No effect")) %>%
  mutate(transcript_ID = str_split(target, "_",simplify= T) [,1]) %>%
  inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version"))


# plot distribution of biotypes
biotype_binom <- binom_model_outputs %>%
  filter(term != "sexmale") %>%
  distinct(target, transcript_biotype, .keep_all = T) %>%
  ggplot(aes(transcript_biotype, fill = transcript_biotype))+
  geom_bar()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")+ # no legend
  stat_count(geom = "Text", aes(label = ..count..), vjust = -0.1) +
  labs(
    title = "Biotypes of genes containing ds introns in the binomial model",
    x = "Biotype",
    y = NULL
  ) 
  



# plot distribution of intron length
intron_length_binom <- binom_model_outputs %>%
  filter(term != "sexmale") %>%
  dplyr::select(target, intron_length, effect) %>%
  ggplot(aes(x = intron_length, fill = effect)) +
  geom_histogram() +
  labs(
    title = "Distribution of intron length of ds introns the binomial model",
    x = "Intron Length",
    y = NULL
  ) +
  # scale_x_continuous(limits = c(70, 30000)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))



# Evaluating the betabinomial model

full_model <- full_model_outputs %>%
  inner_join(intron_length, by = c("target" = "transcript_ID")) %>%
  mutate(effect = case_when(estimate > 0 & adj.p <= 0.05 ~ "Improved SE", 
                            estimate < 0 & adj.p <= 0.05 ~ "Reduced SE" ,
                            estimate < 0 & adj.p > 0.05 ~ "No effect",
                            estimate > 0 & adj.p > 0.05 ~ "No effect")) %>%
  mutate(transcript_ID = str_split(target, "_",simplify= T) [,1]) %>%
  inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version"))



# plot distribution of phenotypes
biotype_full <- full_model %>%
  filter(term != "sexmale") %>%
  distinct(target, transcript_biotype, .keep_all = T) %>%
  ggplot(aes(transcript_biotype, fill = transcript_biotype))+
  geom_bar()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")+ # no legend
  stat_count(geom = "Text", aes(label = ..count..), vjust = -0.1) +
  labs(
    title = "Biotypes of genes containing ds introns in the betabinomial model",
    x = "Biotype",
    y = "Number of introns per biotype"
  ) 


# Plot the intron lengths
intron_length_full <- full_model %>%
  filter(term != "sexmale") %>%
  dplyr::select(target, intron_length, effect) %>%
  ggplot(aes(x = intron_length, fill = effect)) +
  geom_histogram() +
  labs(
    title = "Distribution of Intron Lengths of ds introns in the betabinomial model",
    x = "Intron Length",
    y = "Number of introns"
  ) +
  # scale_x_continuous(limits = c(70, 30000)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")



# full_model %>%
#   filter(term == "Aging") %>%
#   dplyr::select(target, intron_length, effect) %>%
#   ggplot(aes(x = intron_length, fill = effect)) +
#   geom_histogram() +
#   labs(
#     title = "Distribution of Intron Lengths in those affected by aging",
#     x = "Intron Length",
#     y = "Number of Introns"
#   ) +
#   # scale_x_continuous(limits = c(70, 30000)) +
#   theme_minimal() +
#   theme(plot.title = element_text(hjust = 0.5))



sup_Fig1 <-  biotype_full + biotype_binom + intron_length_full +intron_length_binom 
sup_Fig1 +
  plot_annotation(tag_levels = "A") &
  #plot_layout(widths = c(1.1, 1))
  theme(
    plot.tag = element_text(size = 14, face = "bold")
    ,
    plot.tag.position = c(0.08, 0.98)
  ) 





# simplify visualisation by generating gene_intro
# this uses the gene name followed by semicolon and intron_id number
# it makes identifying it easier in charts

df <- binom_model_outputs %>%
  dplyr::filter(term != "sexmale") %>%
  separate(target, into = c(NA, "intron_ID", NA), sep = "_") %>%
  mutate(
    gene_label = ifelse(
      is.na(external_gene_name) | external_gene_name == "",
      ensembl_gene_id,
      external_gene_name
    ), # If gene-name isnt available, use ensembl_gene_id
    gene_intron = paste(gene_label, intron_ID, sep = " : ")
  )%>%
  arrange(gene_label, estimate) %>%
  mutate(gene_intron = factor(gene_intron, levels = unique(gene_intron)))




# visualise the dofferentially spliced introns in the binomial model

ds_binom <- ggplot(df, aes(x = estimate, y =(gene_intron), color = effect)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ term, scales = "free_y") +
  scale_color_manual(values = c("Improved SE" =  colors[6], "Reduced SE" = colors[1]),
                     name = "Effect") +
  labs(
    x = "Effect size",
    y = NULL,
    title = "DS introns due to Aging and Resistance Training (RT)",
    subtitle = "Binomial model (splicing efficiency coded as 0/1)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 8),
    strip.text = element_text(face= "bold"),
    legend.position = "none"
  )



# exstract and visualise data from the beta-binomial model
full <- full_model %>%
  dplyr::filter(term != "sexmale") %>%
  separate(target, into = c(NA, "intron_ID", NA), sep = "_") %>%
  # mutate(gene_intron = paste(external_gene_name, intron_ID, sep = " : "))
  mutate(
    gene_label = ifelse(
      is.na(external_gene_name) | external_gene_name == "",
      ensembl_gene_id,
      external_gene_name
    ),
    gene_intron = paste(gene_label, intron_ID, sep = " : ")
  ) %>%
  arrange(gene_label, estimate) %>%
  mutate(gene_intron = factor(gene_intron, levels = unique(gene_intron)))



ds_beta <- ggplot(full, aes(x = estimate, y = gene_intron, color = effect)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ term, scales = "free_y") +
  scale_color_manual(values = c("Improved SE" = colors[6],
                                "Reduced SE" = colors[1]),
                     name = "Effect") +
  labs(
    x = "Effect size",
    y = NULL,
    title = "DS introns due to Aging and Resistance Training (RT)",
    subtitle = "Beta-binomial model (0,1)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 8),
    strip.text = element_text(face= "bold")
  )




# plot those affected by both aging and exercise in the beta-binomial model

shared <- full %>%
  group_by(gene_intron) %>%
  filter(n_distinct(term) == 2) %>%   # keeps only those seen in both terms
  ungroup()


shared_gene <- full %>%
  group_by(gene_label) %>%
  filter(n_distinct(term) == 2) %>%
  ungroup() %>%
  arrange(gene_label, estimate) %>%   # gene first, then introns
  mutate(gene_intron = factor(gene_intron, levels = unique(gene_intron)))

# plot introns affected by both aging and RT
shared_plot <- ggplot(shared, aes(x = estimate, y = gene_intron, estimate, color = effect)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ term, scales = "free_y") +
  scale_color_manual(values = c("Improved SE" = colors[6],
                                "Reduced SE" = colors[1]),
                     name = "Effect") +
  labs(
    x = "Effect size",
    y = NULL,
    title = "Introns Affected by Both Age and Training"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 8),
        strip.text = element_text(face= "bold"),
        legend.position = "none")



# plot genes affected by both aging and RT

shared_gene_plot <- ggplot(shared_gene, aes(x = estimate, y = gene_intron, color = effect)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ term, scales = "free_y") +
  scale_color_manual(values = c("Improved SE" = colors[6],
                                "Reduced SE" = colors[1]),
                     name = "Effect") +
  labs(
    x = "Effect size",
    y = NULL,
    title = "Genes Containing Introns Affected by Both Age and Training"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 8),
        strip.text = element_text(face= "bold"))





layout <- "
AB
AC
"

final_plot <- ds_beta + ds_binom + shared_plot +
  plot_layout(design = layout, widths = c(1.2,1))

final_plot +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(
   # legend.position = "",
    plot.tag = element_text(size = 14, face = "bold"),
    plot.tag.position = c(0.08, 0.98)
  ) &
  guides(color = guide_legend(title = "Effect"))



# ggsave("Figures/Figure_2.png", bg = colors[4], width = 13, height = 10, dpi = 400)





# quantify the terms in each 
full %>%
  group_by(term, effect) %>%
  summarise(n = n_distinct(gene_label))


df %>%
  group_by(term, effect) %>%
  summarise(n = n_distinct(gene_label))






# extract the differentially spliced introns as a table ranked by their model estimates
aging_table <- full %>%
  filter(term == "Aging" ) %>%
  arrange(desc(abs(estimate)))%>%
  dplyr::select(gene_intron, effect,  estimate, gene_label )

# saveRDS(aging_table, "tables/aging_table.rds")


RT_table <- full %>%
  filter(term == "Resistance Training") %>%
  arrange(desc(abs(estimate))) %>%
  dplyr::select(gene_intron, effect,  estimate, gene_label )

# saveRDS(RT_table, "tables/RT_table.rds")



# Functional annotation of genes affecetd by aging and RT

# First load the gene expression dataset
gene_exp_df <- readRDS("data_new/gene_counts/batch_corrected_genecounts.RDS") 

full_ds_RT <- full %>%
  filter(term == "Resistance Training")

# Functional annotation of the genes affected
ego_RT <- enrichGO(gene =  full_ds_RT$external_gene_name,
                   keyType = "SYMBOL",
                   universe = gene_exp_df$gene_name,
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = T)


## Output results from GO analysis to a table
cluster_RT <- data.frame(ego_RT)

go_RT <- dotplot(ego_RT,
                 
                 font.size = 8, title = "Enriched biological processes in genes containing introns with RT-associated SE") +
  theme(axis.text = element_text(size = 10), axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 10),
        plot.title = element_text(hjust = 0) )

print(go_RT)


# 
# RT_Reduced <- full_ds_RT %>%
#   dplyr::filter(effect == "Reduced SE") %>%
#   arrange(desc(abs(estimate))) #%>%
# # select(gene_intron, effect,  estimate, gene_label )
# 
# RT_Improved <- full_ds_RT %>%
#   dplyr::filter(effect == "Improved SE") %>%
#   arrange(desc(abs(estimate)))
# 
# 
# 
# # Functional annotation of the genes affected
# ego_RT_improved <- enrichGO(gene =  RT_Improved$external_gene_name,
#                             keyType = "SYMBOL",
#                             universe = gene_exp_df$gene_name,
#                             OrgDb = org.Hs.eg.db, 
#                             ont = "BP", 
#                             pAdjustMethod = "BH", 
#                             qvalueCutoff = 0.05, 
#                             readable = T)
# 
# 
# ## Output results from GO analysis to a table
# cluster_RT_improved <- data.frame(ego_RT_improved)
# 
# go_RT_improved <- dotplot(ego_RT_improved,
#                           
#                           font.size = 8, title = "Enriched biological processes in genes containing introns with RT-associated improved SE") +
#   theme(axis.text = element_text(size = 10), axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 10),
#         plot.title = element_text(hjust = 0) )
# 
# print(go_RT_improved)
# 
# 
# 
# 
# # Functional annotation of the genes affected
# ego_RT_reduced <- enrichGO(gene =  RT_Reduced$external_gene_name,
#                            keyType = "SYMBOL",
#                            universe = gene_exp_df$gene_name,
#                            OrgDb = org.Hs.eg.db, 
#                            ont = "BP", 
#                            pAdjustMethod = "BH", 
#                            qvalueCutoff = 0.05, 
#                            readable = T)
# 
# 
# ## Output results from GO analysis to a table
# cluster_RT_reduced <- data.frame(ego_RT_reduced)
# 
# go_RT_reduced <- dotplot(ego_RT_reduced,
#                          
#                          font.size = 8, title = "Enriched biological processes in genes containing introns with RT-associated improved SE") +
#   theme(axis.text = element_text(size = 10), axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 10),
#         plot.title = element_text(hjust = 0) )
# 
# print(go_RT_reduced)

ggsave("Figures/GO_RT.png", bg = colors[4], scale=2.5, dpi = 400)

full_ds_aging <- full %>%
  filter(term == "Aging")

# Functional annotation of the genes affected
ego_aging <- enrichGO(gene =  full_ds_aging$external_gene_name,
                      keyType = "SYMBOL",
                      universe = gene_exp_df$gene_name,
                      OrgDb = org.Hs.eg.db, 
                      ont = "BP", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05, 
                      readable = T)


## Output results from GO analysis to a table
cluster_aging <- data.frame(ego_aging)

go_aging <- dotplot(ego_aging,
                    
                    font.size = 8, title = "Enriched biological processes in genes containing introns with aging-associated SE") +
  theme(axis.text = element_text(size = 10), axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 10),
        plot.title = element_text(hjust = 0) )

print(go_aging)

ggsave("Figures/GO_aging.png", bg = colors[4],height = 20, width = 15, dpi = 400)




Fig4 <-  go_RT + go_aging 
Fig4 +
  plot_annotation(tag_levels = "A") &
  #plot_layout(widths = c(1.1, 1))
  theme(
    plot.tag = element_text(size = 14, face = "bold")
    ,
    plot.tag.position = c(0.08, 0.98)
  ) 


ggsave("Figures/Figure_4.png", bg = colors[4], scale = 4, dpi = 400)

# Plot the top differentially expressed genes
top_age_introns <- full_model_outputs %>%
  dplyr::filter(term == "Aging") %>%
  dplyr::arrange(desc(abs(estimate)))%>%
  dplyr::slice(1:5) 

age_introns_df <- all_splice_df %>%
   dplyr::filter(transcript_ID %in% top_age_introns$target) %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID) ) %>%
  inner_join(top_age_introns, by = c("transcript_ID" = "target")) %>%
  separate(transcript_ID, into = c("transcript_name", "intron_ID", NA), sep = "_", remove = F) %>%
  inner_join(gene_annotation, by= c("transcript_name" = "ensembl_transcript_id_version")) %>%
  
  mutate(
    gene_label = ifelse(
      is.na(external_gene_name) | external_gene_name == "",
      ensembl_gene_id,
      external_gene_name
    ),
    gene_intron = paste(gene_label, intron_ID, sep = " : ")) %>%
  inner_join(all_full_metadata, by = "seq_sample_id")%>% 
  group_by(scaled_age, gene_intron, transcript_ID) %>%
  summarise(mean_SE = mean(SE, na.rm = TRUE), .groups = "drop")


aging_df <- ggplot(age_introns_df, aes(x = scaled_age, y = mean_SE, group = gene_intron, color = gene_intron)) +
  geom_point(alpha = 1, size = 2) +
  geom_smooth(aes(color = gene_intron), method = "lm", se = FALSE, size = 0.5) +
  theme_minimal() +
  labs(
    y = "Splicing efficiency",
    x = "scaled age of participants",
    title = "Top 5 introns with aging-associated changes in SE"
  )+
  theme(plot.title = element_text(hjust = 0.2),
        strip.text = element_text(face= "bold"))



# explore the top five RT associated introns

RT_introns <- full_model_outputs %>%
  dplyr::filter(term == "Resistance Training") %>%
  dplyr::arrange(desc(abs(estimate)))%>%
  dplyr::slice(1:5) # %>%
# pull(target)



RT_introns_df <- all_splice_df %>%
  dplyr::filter(transcript_ID %in% RT_introns$target) %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID) ) %>%
  inner_join(RT_introns, by = c("transcript_ID" = "target")) %>%
  separate(transcript_ID, into = c("transcript_ID", "intron_ID", NA), sep = "_") %>%
  inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version")) %>%
  
  mutate(
    gene_label = ifelse(
      is.na(external_gene_name) | external_gene_name == "",
      ensembl_gene_id,
      external_gene_name
    ),
    gene_intron = paste(gene_label, intron_ID, sep = " : ")) %>%
  inner_join(all_full_metadata, by = "seq_sample_id") %>%
  group_by(time, gene_intron) %>%
  summarise(mean_SE = mean(SE, na.rm = TRUE), .groups = "drop")


RT_df <- ggplot(RT_introns_df, aes(x = time, y = mean_SE, colour = gene_intron)) +
  geom_point(size = 3) +
 # geom_jitter(alpha = 1, width = 0.5, size = 1) +
  theme_minimal()+
  labs(
    y = "Splicing efficiency",
    x = "Time",
    title = "Top 5 introns with RT-associated changes in SE"
  )+
  theme(plot.title = element_text(hjust = 0.1),
        strip.text = element_text(face= "bold"))




# exploring if the expression of genes follow the ds patterns



RT_expression <- full_model %>%
  dplyr::filter(target %in% RT_introns$target)

exp_df <- gene_exp_df %>%
  dplyr::filter(gene_name %in% RT_expression$external_gene_name) %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "gene_count",
               cols = -(gene_name) ) %>%
  inner_join(all_full_metadata, by = "seq_sample_id") %>%
  group_by(time, gene_name) %>%
  summarise(mean_count = mean(gene_count, na.rm = TRUE), .groups = "drop")







train <- ggplot(exp_df, aes(x = time, y = mean_count, colour= gene_name)) +
  geom_point(size = 3) +
  theme_minimal()+
  labs(
    y = "Gene expression value",
    x = "Time",
    title = "Genes containing top 5 ds introns due to training"
  )+
  theme(plot.title = element_text(hjust = 0.1),
        strip.text = element_text(face= "bold"))




# do same for the aging affected ones

aging_expression <- full_model %>%
  dplyr::filter(target %in% top_age_introns$target)


age_exp_df <- gene_exp_df %>%
  dplyr::filter(gene_name %in% aging_expression$external_gene_name) %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "gene_count",
               cols = -(gene_name) ) %>%
  inner_join(all_full_metadata, by = "seq_sample_id") %>%
  group_by(scaled_age, gene_name) %>%
  summarise(mean_count = mean(gene_count, na.rm = TRUE), .groups = "drop")


age <- ggplot(age_exp_df, aes(x = scaled_age, y = mean_count, colour= gene_name)) +
  geom_point(alpha = 1, size = 2) +
  geom_smooth(aes(color = gene_name), method = "lm", se = FALSE, size = 0.5) +
  theme_minimal() +
  labs(
    y = "Gene expression value",
    x = "scaled age of participants",
    title = "Genes containing top 5 ds introns due to aging"
  )+
  theme(plot.title = element_text(hjust = 0.2),
        strip.text = element_text(face= "bold"))



Fig2_plot <- (aging_df | age) /
  plot_spacer() / 
  (RT_df | train)
Fig2_plot +
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(1,0.15, 1))
  theme(
    plot.tag = element_text(size = 14, face = "bold"),
    plot.tag.position = c(0.08, 0.98)
  ) 

# ggarrange( aging_df, age,  RT_df ,train,
#            labels = c("A", "B", "C", "D"),
#            font.label = list(size = 8), 
#            align = "hv",
#            label.x = 0.05)

ggsave("Figures/Figure_3.png", bg = colors[4], scale = 2.5, dpi = 400)
