library(dplyr)
library(ggplot2)
library(broom.mixed)
library(purrr)
library(tidyr)
library(seqwrap)
library(stringr)

## Color scale ##
colors <-  c("#d7191c", "#fdae61", "#ffffbf", "#abd9e9", "#2c7bb6", "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")

# Load the gene annotation file
gene_annotation <- readRDS("data/ensembl_gene_annotation.RDS")

# Load the spliing data
all_splice_df <- readRDS("data/all_splice.RDS") %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) 


# Load the metadata 
metadata <- readRDS("data/all_full_metadata.RDS") 

# Load one file from which we will extract intron length
# This is valid as only introns quantified in all samples were included in the analyses
intron_length <- readr::read_tsv("data_new/Alpha_Omega_SpliceQ_outputs/A_102.tsv") %>%
  
  distinct(across(6:ncol(.)), .keep_all = T) %>% # Removes duplicates based on columns 6 to end
  mutate(transcript_ID = paste0(transcript_ID, "_", intron_ID, "_", chr),
         intron_length = abs((sj3start - sj5end) + 1)) %>%
  group_by(gene_ID) %>%
  mutate(number_introns = n()) %>%
  ungroup

# Load the betabinomial spline model 

Beta_model <- readRDS("data/splined_beta_binomial_model.RDS") %>%
  seqwrap_summarise()

# extract the effect of training
beta_model_outputs <- Beta_model$summaries %>% 
  dplyr::select(-group) %>%
 # dplyr::filter(term == "timePostExc") %>%
  drop_na() %>%
  mutate(adj.p = p.adjust(p.value, method = "fdr"),
         term = recode(term,
                       "timePostExc" = "Resistance Training"),
         effect = case_when(estimate > 0 & adj.p <= 0.05 ~ "Improved SE", 
                            estimate < 0 & adj.p <= 0.05 ~ "Reduced SE" ,
                            estimate < 0 & adj.p > 0.05 ~ "No effect",
                            estimate > 0 & adj.p > 0.05 ~ "No effect"),
         transcript_ID = str_split(target, "_",simplify= T) [,1],
         sig = adj.p <= 0.05,
         neg_log10_fdr = -log10(adj.p)) %>% 
  inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version")) %>%
  separate(target, into = c(NA, "intron_ID", NA), sep = "_", remove = F) %>%
  # create a gene and intron label using the gene name and intron number
  mutate(
    gene_label = ifelse(
      is.na(external_gene_name) | external_gene_name == "",
      ensembl_gene_id,
      external_gene_name
    ),
    gene_intron = paste(gene_label, intron_ID, sep = " : "),
    rank_score = -log10(adj.p) * abs(estimate)) %>%
  inner_join(intron_length, by = c("ensembl_gene_id_version" = "gene_ID")) %>%
  arrange(desc(rank_score)) #%>%
  dplyr::select(target, effect, estimate, adj.p, sig, neg_log10_fdr, transcript_biotype,
                external_gene_name, gene_label, gene_intron,
                rank_score, intron_length, number_introns)




summary <- beta_model_outputs %>%
  summarise(
    n_total = n(),
    n_sig = sum(sig),
    perc_sig = 100 * n_sig / n_total,
    x = max(estimate, na.rm = TRUE),
    y = max(neg_log10_fdr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      "Ds introns: ",
      n_sig, "/", n_total,
      " (", round(perc_sig, 1), "%)"
    )
  )



top9_introns <- beta_model_outputs %>% slice_head(n = 9)


top9_genes <- unique(top9_introns$gene_label)

RT_introns_df <- all_splice_df %>%
  dplyr::filter(transcript_ID %in% top9_introns$target) %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID) ) %>%
  inner_join(top9_introns, by = c("transcript_ID" = "target")) %>%
  inner_join(metadata, by = "seq_sample_id")%>% 
  group_by(scaled_age, gene_intron, transcript_ID) %>%
  summarise(mean_SE = mean(SE, na.rm = TRUE), .groups = "drop")

# Create SE plots, one per intron




intron_plots <- lapply(unique(RT_introns_df$gene_intron), function(intron) {
  ggplot(
    RT_introns_df %>% filter(gene_intron == intron),
    aes(x = scaled_age, y = mean_SE)
  ) +
    geom_point(size = 1.5, alpha = 0.9, colour = "black") +
    geom_smooth(method = "lm", se = FALSE, colour = "red", linewidth = 0.6) +
    labs(title = intron, x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      axis.text = element_text(size = 10),
      axis.title = element_blank()
    )
})


# stack them into a panel


mini_panel <- plot_grid(
  plotlist = intron_plots,
  ncol = 3
)


volcano_plot <- ggplot(beta_model_outputs , aes(estimate, neg_log10_fdr, colour = effect)) +
  geom_point(aes(colour = effect), alpha = 0.7, size = 2) +
  
  geom_text_repel(
    data = top6_introns,
    aes(label = gene_intron),
    size = 6,
    max.overlaps = Inf
  ) +
  
  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  
  scale_color_manual(values = c("Improved SE" =  colors[6], "Reduced SE" = colors[1]),
                     name = "Effect") +
  
  geom_text(
    data = summary,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.1,
    vjust = 1.0,
    size = 5,
    fontface = "italic"
  ) +
  
  labs(
    title = "Differentially spliced(ds) introns due to Aging",
    subtitle = "Beta-binomial model (splicing efficiency coded as 0,1)",
    x = "Effect size",
    y = expression(-log[10]("FDR value"), clip = "off")
  ) +
  
  # theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 13),
    legend.title = element_blank(),
    legend.text = element_text(size = 9, face = "bold"), 
    
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"), 
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )


# combine volcano + insets

final_volcano <- ggdraw() +
  draw_plot(volcano_plot) +
  draw_plot(
    mini_panel,
    x = 0.55,   # adjust horizontally
    y = 0.45,  # adjust vertically
    width = 0.38,
    height = 0.48
  )

final_volcano

# Betabinomial predictions from marginal effects

Beta_pred <- readRDS("data/splined_UPDATED_beta_binomial_predictions.RDS")

binomial_pred <- readRDS("data/splined_UPDATED_binomial_predictions.RDS")

binom_time_effects <- readRDS("data/splined_binomial_time_effect.RDS") %>%
  group_by(target, scaled_age) %>%
  mutate(adj.p = p.adjust(p.value, method = "BH")) %>%
  filter(adj.p <= 0.05)


# Find introns where the exercise training effect is significantly active at ANY age milestone

time_effects <- readRDS("data/splined_beta_binomial_time_effect.RDS") %>%
  group_by(target, scaled_age) %>%
  mutate(adj.p = p.adjust(p.value, method = "BH")) %>%
  filter(adj.p <= 0.05)

  
# time_effects  <- readRDS("data/splined_beta_binomial_time_effect.RDS") %>%
#   group_by(target) %>%
#   summarise(
#     min_p = min(p.value),
#     adj.p = p.adjust(min_p, method = "BH")) %>%
#   filter(adj.p <= 0.05)
length(unique(time_effects$target))

time_effects %>%
  summarise(
    mean_abs = mean(abs(estimate)),
    max_abs = max(abs(estimate))
  )



# visualise a global curve showing the non-linear relationship between aging and intron splicing efficiency

Beta_pred %>%
  group_by(time, scaled_age) %>%
  summarise(mean_IR = mean(estimate)) %>%
  ggplot(aes(scaled_age, mean_IR, color = time)) +
  geom_line() +
  theme_bw()+
  labs(
    title = "Global curve of Intron retention using betabinomial model",
    x = "Scale participants' age",
    y = "Splicing Efficiency"
  ) 

  
binomial_pred %>%
  group_by(time, scaled_age) %>%
  summarise(mean_IR = mean(estimate)) %>%
  ggplot(aes(scaled_age, mean_IR, color = time)) +
  geom_line() +
  theme_bw()+
  labs(
    title = "Global curve of Intron retention using binomial model",
    x = "Scale participants' age",
    y = "Splicing Efficiency"
  ) 

  # Repeat this for only introns that are fdr-significant
unique(Beta_pred$p.value)

Beta_pred %>%
  filter(target %in% time_effects$target) %>%
  group_by(time, scaled_age) %>%
  summarise(mean_IR = mean(estimate)) %>%
  ggplot(aes(scaled_age, mean_IR, color = time)) +
  geom_line() +
  theme_bw()+
  labs(
    title = "Global curve of Intron retention using betabinomial model",
    x = "Scale participants' age",
    y = "Splicing Efficiency"
  ) 


binomial_pred %>%
  filter(target %in% binom_time_effects$target) %>%
  group_by(time, scaled_age) %>%
  summarise(mean_IR = mean(estimate)) %>%
  ggplot(aes(scaled_age, mean_IR, color = time)) +
  geom_line() +
  theme_bw()+
  labs(
    title = "Global curve of Intron retention using binomial model",
    x = "Scale participants' age",
    y = "Splicing Efficiency"
  ) 



old_age_responders <- time_effects %>%
  # 1. Correct p-values globally within each specific age group
  group_by(scaled_age) %>%
  mutate(adj.p = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  # 2. Filter for robust significance (FDR < 0.05) strictly at advanced milestones
  filter(adj.p < 0.05 & scaled_age >= 0.75) # %>%
  # pull(target) %>%
  # unique()

# Find targets significant in youth to subtract them (ensuring age-specificity)
young_age_responders <- time_effects %>%
  group_by(scaled_age) %>%
  mutate(adj.p = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  filter(adj.p < 0.05 & scaled_age == 0.25) #%>%
  # pull(target) %>%
  # unique()

length(unique(old_age_responders$target))

length(unique(young_age_responders$target))


Beta_pred %>%
  filter(target %in% young_age_responders$target) %>%
  group_by(time, scaled_age) %>%
  summarise(mean_IR = mean(estimate)) %>%
  ggplot(aes(scaled_age, mean_IR, color = time)) +
  geom_line() +
  theme_bw()+
  labs(
    title = "Global curve of Intron retention using betabinomial model",
    x = "Scale participants' age",
    y = "Splicing Efficiency"
  ) 

# detect peaks in betabinomial 

slopes_beta <- readRDS("data/splined_UPDATED_beta_binomial_age_slopes.RDS")
