library(dplyr)
library(tidyverse)
library(seqwrap)
library(glmmTMB)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(scales)
library(biomaRt)
library(patchwork)
library(ggrepel)


# source("R/Trainome_functions.R")
# load the individual metadata 
#COPD metadata
copd_metadata <- readRDS("data/copd_metadata.RDS") %>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)


# Volume_data
Vol_metadata <- readRDS("data/volume_metadata.RDS")%>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age) 


# Contratrain_data
ct_metadata <- readRDS("data/contratrain_metadata.RDS") %>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)


# Alpha and Omega data

A_Omega_metadata <- readRDS("data/Alpha_Omega_metadata.RDS")%>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)



Relief_full_meta <- readRDS("data/Relief_metadata.RDS")%>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age)


# Merge the metadata into one
metadata <- rbind(copd_metadata, Vol_metadata)%>%
  rbind(ct_metadata) %>%
  rbind(A_Omega_metadata) %>%
  rbind(Relief_full_meta) %>%
  mutate(across(c("age"), round, 0)) %>%
  mutate(sex = factor(sex, levels = c("female", "male")),
         time = factor(time, levels = c("PreExc", "PostExc")),
         scaled_age = round(rescale(age), digits = 2),
         participant = paste0(study, "_", participant))







# Load the individual splicing data
copd_splice_df <- readRDS("data/copd_splicing_data.RDS") %>%
  drop_na()

vol_splice_df <- readRDS("data/volume_splicing_data.RDS") %>%
  drop_na()
ct_splice_df <- readRDS("data/contratrain_splicing_data.RDS") %>%
  drop_na()

AOD_splice_df <- readRDS("data/Alpha_Omega_splicing_data.RDS") %>%
  drop_na()
Relief_full_splice <- readRDS("data/Relief_splicing_data.RDS") %>%
  drop_na()





all_splice_df <-copd_splice_df  %>%
  inner_join(vol_splice_df, by = "transcript_ID") %>%
  inner_join(ct_splice_df, by = "transcript_ID") %>%
  inner_join(AOD_splice_df, by = "transcript_ID") %>%
  inner_join(Relief_full_splice, by = "transcript_ID") 




# select only columnames in the splicing data that match sequence ids in the metadata
intersect <- intersect(colnames(all_splice_df), metadata$seq_sample_id)

all_splice_df <-all_splice_df %>%
  subset(select = c("transcript_ID", intersect)) %>%
  drop_na()

# IDs in metadata not found in splice data
# missing_from_splice <- setdiff(metadata$seq_sample_id, colnames(all_splice_df))
# 
# # IDs in splice data not found in metadata
# missing_from_meta <- setdiff(colnames(all_splice_df)[-1], metadata$seq_sample_id)
# 
# head(missing_from_splice)
# head(missing_from_meta)

# saveRDS(all_splice_df, "data/Trainome_all_splice_df.RDS")
# 
# saveRDS(metadata, "data/Trainome_metadata.RDS")

# REORDER THE SEQUENCE ID TO MATCH BOTH DATAFRAMMES
all_splice_reordered <- all_splice_df[, c("transcript_ID",metadata$seq_sample_id)] 

# Check if everything matches except the transcript_id
match(colnames(all_splice_reordered), metadata$seq_sample_id)


# visualise the data

## Color scale ##
colors <-  c("#d7191c", "#fdae61", "#ffffbf", "#abd9e9", "#2c7bb6", "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")


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
    hjust = 0.01,
    vjust = 0.5,
    size = 4,
    fontface= "bold.italic"
  )+
  # facet_wrap(~ study, ncol = 2) +
  theme_minimal(base_size = 12, base_family = "Arial") +
  labs(
    title = "Age distribution of all participants",
    x = "Age",
    y = NULL
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face= "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12, face = "bold"), 
    axis.text.y = element_text(size = 12, face= "bold"),
    axis.text.x = element_text(size = 12, face= "bold"),
    axis.title = element_text(size = 12, face= "bold")
  )






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
  facet_wrap(~ study, ncol = 2, scales = "free") +
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
  theme_minimal(base_size = 12, base_family = "Arial") +
  labs(
    title = " Distribution of Participants across all studies",
    x = "Age",
    y = "Number of participants"
  ) +
  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face= "bold"),
    # legend.title = element_blank(),
    # legend.text = element_text(size = 12, face = "bold"), 
    legend.position = "none",
    axis.text.y = element_text(size = 12, face= "bold"),
    axis.text.x = element_text(size = 12, face= "bold"),
    axis.title = element_text(size = 12, face= "bold")
  )


# Load the bubble image
# bubble_img <- png::readPNG("Figures/Bubble_chart.PNG")
# image_grob <- grid::rasterGrob(bubble_img, interpolate = TRUE)
# 
# image_plot <- wrap_elements(image_grob)

# combine with the histograms from the participants' displa



Figure_1 <- by_study + all  

  
Figure_1 +
  plot_annotation(tag_levels = "A") &
  #plot_layout(widths = c(1.1, 1))
  theme(
    plot.tag = element_text(size = 14, face = "bold")
    ,
    plot.tag.position = c(0.08, 0.98)
  ) 




ggsave("Figures/Trainome_Figure_2.png", bg = colors[4], width = 15, height = 10, dpi = 400)
# Build binomial model
# This model investigates the question, "given an intron,
# what is the probability of perfect splicing as a function of age and resistance exercise training"



# derive a matrix that indicates 0 if SE is not 1
one_inflated_mat <- all_splice_reordered

one_inflated_mat[-1] <- lapply(
  one_inflated_mat[-1],
  function(x) as.integer(x == 1)
)





binom_container <- seqwrap_compose(
  data       = one_inflated_mat,
  metadata   = metadata,
  samplename = "seq_sample_id",
  modelfun   = glmmTMB::glmmTMB,
  arguments  =  alist(
    formula = y ~ splines::ns(scaled_age, df = 4) + time + sex +
      (1 | study) + (1 | participant),
    family  = binomial(link = "logit")
  )
)


binom_results <- seqwrap(
  binom_container,
  return_models = TRUE,
  cores = 10,
  # subset = 1:200
)

saveRDS(binom_results, "data/splined_binom_models.RDS")




# The second model
# Models the degree of intron retention among introns that are not perfectly spliced.
# Flip SE so perfect splicing (1) becomes structural zeros for the zi component.

zi_mat <- all_splice_reordered

zi_mat[-1] <- 1 - all_splice_reordered[-1]

zi_container <- seqwrap_compose(
  data       = zi_mat,
  metadata   = metadata,
  samplename = "seq_sample_id",
  modelfun   = glmmTMB::glmmTMB,
  arguments  = alist(
    formula   = y ~ splines::ns(scaled_age, df = 4) + time + sex +
                    (1 | study) + (1 | participant),
    ziformula = ~ splines::ns(scaled_age, df = 4) + time,
    family    = glmmTMB::beta_family(link = "logit")
  )
)

zi_results <- seqwrap(
  zi_container,
  return_models = TRUE,
  cores = 10
)



saveRDS(zi_results, "data/splined_zi_results.RDS")
# zi_results <- readRDS("data/splined_zi_results.RDS")
 
zi_container_int <- seqwrap_compose(
  data       = zi_mat,
  metadata   = metadata,
  samplename = "seq_sample_id",
  modelfun   = glmmTMB::glmmTMB,
  arguments  = alist(
    formula   = y ~ splines::ns(scaled_age, df = 4) * time + sex +
      (1 | study) + (1 | participant),
    ziformula = ~ splines::ns(scaled_age, df = 4) * time,
    family    = glmmTMB::beta_family(link = "logit")
  )
)

zi_results_int <- seqwrap(
  zi_container_int,
  return_models = TRUE,
  cores = 10
)

saveRDS(zi_results_int, "data/splined_zi_interaction_results.RDS")
# a Betabinomial model that investigates interaction


















# Load the gene annotation file
gene_annotation <- readRDS("data/ensembl_gene_annotation.RDS")

# Load one file from which we will extract intron length
# This is valid as only introns quantified in all samples were included in the analyses
intron_length <- readr::read_tsv("data_new/Alpha_Omega_SpliceQ_outputs/A_102.tsv") %>%
  
  distinct(across(6:ncol(.)), .keep_all = T) %>% # Removes duplicates based on columns 6 to end
  mutate(transcript_ID = paste0(transcript_ID, "_", intron_ID, "_", chr),
         intron_length = abs((sj3start - sj5end) + 1) ) %>% # Ensures positive length regardless of strand
  dplyr::select(transcript_ID, intron_length)

#  load the gene expression dataset
gene_exp_df <- readRDS("data_new/gene_counts/batch_corrected_genecounts.RDS") 

 
 # extract the model summaries in the beta binomial model
 informed_binom_sum <- seqwrap_summarise(binom_results)
 
 
 # filter the summary and create new columns
 binom_model_outputs <- informed_binom_sum$summaries %>% 
   dplyr::select(-group) %>%
   inner_join(intron_length, by = c("target" = "transcript_ID")) %>%
   filter(term != "(Intercept)", term != "sexmale") %>%
   drop_na() %>%
   group_by(term) %>%
   mutate(
     adj.p = p.adjust(p.value, method = "fdr"),
     term = recode(term,
                   "scaled_age" = "Aging",
                   "timePostExc" = "Resistance Training"),
     effect = case_when(estimate > 0 & adj.p <= 0.05 ~ "Improved SE", 
                        estimate < 0 & adj.p <= 0.05 ~ "Reduced SE" ,
                        estimate < 0 & adj.p > 0.05 ~ "No effect",
                        estimate > 0 & adj.p > 0.05 ~ "No effect"),
     transcript_ID = str_split(target, "_",simplify= T) [,1]) %>%
   ungroup() %>%
   mutate(
     sig = adj.p <= 0.05,
     neg_log10_fdr = -log10(adj.p)
   ) %>%
   inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version")) %>%
   separate(target, into = c(NA, "intron_ID", NA), sep = "_", remove = F) %>%
   mutate(
     gene_label = ifelse(
       is.na(external_gene_name) | external_gene_name == "",
       ensembl_gene_id,
       external_gene_name
     ), # If gene-name isnt available, use ensembl_gene_id
     gene_intron = paste(gene_label, intron_ID, sep = " : "))
 # %>%
 #   arrange(gene_label, estimate) %>%
 #   mutate(gene_intron = factor(gene_intron, levels = unique(gene_intron)))
 
 
 
 # extract the results of the non-binarised model
 
 full_model_sum <- seqwrap_summarise(full_model)
 
 # The non-binarised model 
 full_model_outputs <- full_model_sum$summaries %>% 
   dplyr::select(-group) %>%
   inner_join(intron_length, by = c("target" = "transcript_ID")) %>%
   filter(term != "(Intercept)", term != "sexmale") %>%
   drop_na() %>%
   group_by(term) %>%
   mutate(
     adj.p = p.adjust(p.value, method = "fdr"),
     term = recode(term,
                   "scaled_age" = "Aging",
                   "timePostExc" = "Resistance Training"),
     effect = case_when(estimate > 0 & adj.p <= 0.05 ~ "Improved SE", 
                        estimate < 0 & adj.p <= 0.05 ~ "Reduced SE" ,
                        estimate < 0 & adj.p > 0.05 ~ "No effect",
                        estimate > 0 & adj.p > 0.05 ~ "No effect"),
     transcript_ID = str_split(target, "_",simplify= T) [,1]) %>%
   ungroup() %>%
   mutate(
     sig = adj.p <= 0.05,
     neg_log10_fdr = -log10(adj.p)
   ) %>%
   inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version")) %>%
   separate(target, into = c(NA, "intron_ID", NA), sep = "_", remove = F) %>%
   mutate(
     gene_label = ifelse(
       is.na(external_gene_name) | external_gene_name == "",
       ensembl_gene_id,
       external_gene_name
     ), # If gene-name isnt available, use ensembl_gene_id
     # simplify visualisation by generating gene_intro
     # this uses the gene name followed by semicolon and intron_id number
     # it makes identifying it easier in charts
     gene_intron = paste(gene_label, intron_ID, sep = " : "))
 
 
 
 
 
 
 # BELOW IS EXPLORATION OF THE DATA OUTPUTS
 
 
 # plot distribution of biotypes
 biotype_binom <- binom_model_outputs %>%
   filter(adj.p <= 0.05) %>%
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
   filter(adj.p <= 0.05) %>%
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
   filter(adj.p <= 0.05)
 
 
 
 # plot distribution of phenotypes
 biotype_full <- full_model %>%
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
 
 
 
 
 
 sup_Fig1 <-  biotype_full + biotype_binom + intron_length_full +intron_length_binom 
 sup_Fig1 +
   plot_annotation(tag_levels = "A") &
   #plot_layout(widths = c(1.1, 1))
   theme(
     plot.tag = element_text(size = 14, face = "bold")
     ,
     plot.tag.position = c(0.08, 0.98)
   ) 
 
 
 
 
 

 
 
 
 
 
 # Volcano plots
 
 # visualise the dofferentially spliced introns in the binomial model
 ggplot(binom_model_outputs, aes(estimate, neg_log10_fdr, colour = effect)) +
  geom_point(aes(colour = effect), alpha = 0.7, size = 2) +
  
  # geom_text_repel(
  #   data = beta_top10_labels,
  #   aes(label = gene_intron),
  #   size = 4,
  #   max.overlaps = Inf
  # ) +
  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  
  scale_color_manual(values = c("Improved SE" =  colors[6], "Reduced SE" = colors[1]),
                     name = "Effect") +
  
  facet_wrap(~term, scales = "free") +
 # coord_cartesian(xlim = c(-1, 1))+
  
  labs(
    title = "Differentially spliced introns due to Aging and Resistance Training",
    subtitle = "Binomial model (splicing efficiency coded as 0/1)",
    x = "Effect size",
    y = expression(-log[10]("FDR value"), clip = "off")
  ) +
  
  # theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 13),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, face = "bold"), 
    
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"), 
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
    )


 
 
 
 
 # Volcano plot full model
 
 # visualise the dofferentially spliced introns in the binomial model
 ggplot(full_model_outputs, aes(estimate, neg_log10_fdr, colour = effect)) +
   geom_point(aes(colour = effect), alpha = 0.7, size = 2) +
   
   
   geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
   
   scale_color_manual(values = c("Improved SE" =  colors[6], "Reduced SE" = colors[1]),
                      name = "Effect") +
   
   facet_wrap(~term, scales = "free") +
   # coord_cartesian(xlim = c(-1, 1))+
   
   labs(
     title = "Differentially spliced introns due to Aging and Resistance Training",
     subtitle = "Beta-binomial model (splicing efficiency coded as 0,1)",
     x = "Effect size",
     y = expression(-log[10]("FDR value"), clip = "off")
   ) +
   
   # theme_minimal(base_size = 16) +
   theme(
     plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
     plot.subtitle = element_text(hjust = 0.5, size = 13),
     legend.title = element_blank(),
     legend.text = element_text(size = 14, face = "bold"), 
     
     axis.title.x = element_text(size = 14, face = "bold"),
     axis.title.y = element_text(size = 14, face = "bold"),
     strip.text = element_text(size = 14, face = "bold"), 
     axis.text.x = element_text(size = 14, face = "bold"),
     axis.text.y = element_text(size = 14, face = "bold"),
     plot.background  = element_rect(fill = "white", colour = NA),
     panel.background = element_rect(fill = "white", colour = NA)
   )
 
 
 
 
# Visualize the effect plots
# binomial model
 ds_binom <- ggplot((binom_model_outputs %>%
                       filter(adj.p <= 0.05) %>%
                       arrange(gene_label, estimate)), aes(x = estimate, y =(gene_intron), color = effect)) +
   geom_point(size = 4) +
   geom_vline(xintercept = 0, linetype = "dashed") +
   facet_wrap(~ term, scales = "free") +
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
     plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
     plot.subtitle = element_text(hjust = 0.5, size = 13),
     # legend.title = element_blank(),
     # legend.text = element_text(size = 14, face = "bold"), 
     # 
     axis.title.x = element_text(size = 14, face = "bold"),
     axis.title.y = element_text(size = 14, face = "bold"),
     strip.text = element_text(size = 14, face = "bold"), 
     axis.text.x = element_text(size = 10, face = "bold"),
     axis.text.y = element_text(size = 10, face = "bold"),
     plot.background  = element_rect(fill = "white", colour = NA),
     panel.background = element_rect(fill = "white", colour = NA),
     legend.position = "none"
   )
 
 
 
 # betabinomial model
 ds_beta <- ggplot((full_model %>%
                     arrange(gene_label, estimate)) , aes(x = estimate, y = gene_intron, color = effect)) +
   geom_point(size = 3) +
   geom_vline(xintercept = 0, linetype = "dashed") +
   facet_wrap(~ term, scales = "free") +
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
     plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
     plot.subtitle = element_text(hjust = 0.5, size = 13),
     # legend.title = element_blank(),
     # legend.text = element_text(size = 14, face = "bold"), 
      
     axis.title.x = element_text(size = 14, face = "bold"),
     axis.title.y = element_text(size = 10, face = "bold"),
     strip.text = element_text(size = 14, face = "bold"), 
     axis.text.x = element_text(size = 12, face = "bold"),
     axis.text.y = element_text(size = 8, face = "bold"),
     plot.background  = element_rect(fill = "white", colour = NA),
     panel.background = element_rect(fill = "white", colour = NA),
     legend.position = "none"
   )
 
 
 
 
 
 # plot those affected by both aging and exercise in the beta-binomial model
 
 shared <- full_model %>%
   group_by(gene_intron) %>%
   filter(n_distinct(term) == 2) %>%   # keeps only those seen in both terms
   ungroup()
 
 
 shared_gene <- full_model %>%
   group_by(gene_label) %>%
   filter(n_distinct(term) == 2) %>%
   ungroup() %>%
   arrange(gene_label, estimate) %>%   # gene first, then introns
   mutate(gene_intron = factor(gene_intron, levels = unique(gene_intron)))
 
 # plot introns affected by both aging and RT
 shared_plot <- ggplot(shared, aes(x = estimate, y = gene_intron, estimate, color = effect)) +
   geom_point(size = 4) +
   geom_vline(xintercept = 0, linetype = "dashed") +
   facet_wrap(~ term, scales = "free") +
   scale_color_manual(values = c("Improved SE" = colors[6],
                                 "Reduced SE" = colors[1]),
                      name = "Effect") +
   labs(
     x = "Effect size",
     y = NULL,
     title = "Introns Affected by Both Age and Training"
   ) +
   theme_minimal() +
   theme(
     plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
     plot.subtitle = element_text(hjust = 0.5, size = 13),
     # legend.title = element_blank(),
     # legend.text = element_text(size = 14, face = "bold"), 
     # 
     axis.title.x = element_text(size = 14, face = "bold"),
     axis.title.y = element_text(size = 14, face = "bold"),
     strip.text = element_text(size = 14, face = "bold"), 
     axis.text.x = element_text(size = 10, face = "bold"),
     axis.text.y = element_text(size = 10, face = "bold"),
     plot.background  = element_rect(fill = "white", colour = NA),
     panel.background = element_rect(fill = "white", colour = NA),
     legend.position = "none"
   )
 
 
 # plot genes affected by both aging and RT
 
 # shared_gene_plot <- ggplot(shared_gene, aes(x = estimate, y = gene_intron, color = effect)) +
 #   geom_point(size = 3) +
 #   geom_vline(xintercept = 0, linetype = "dashed") +
 #   facet_wrap(~ term, scales = "free_y") +
 #   scale_color_manual(values = c("Improved SE" = colors[6],
 #                                 "Reduced SE" = colors[1]),
 #                      name = "Effect") +
 #   labs(
 #     x = "Effect size",
 #     y = NULL,
 #     title = "Genes Containing Introns Affected by Both Age and Training"
 #   ) +
 #   theme_minimal() +
 #   theme(axis.text.y = element_text(size = 8),
 #         plot.title = element_text(hjust = 0.5),
 #         plot.subtitle = element_text(hjust = 0.5, size = 8),
 #         strip.text = element_text(face= "bold"))
 # 
 # 
 
 
 
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
 # aging_table <- full_model %>%
 #   filter(term == "Aging" ) %>%
 #   arrange(desc(abs(estimate)))%>%
 #   dplyr::select(gene_intron, effect,  estimate, gene_label )
 
 aging_table <- full_model %>%
   filter(term == "Aging") %>%
   mutate(rank_score = -log10(adj.p) * abs(estimate)) %>%
   arrange(desc(rank_score)) #%>%
   #dplyr::select(gene_intron, effect, estimate, adj.p, gene_label, rank_score)
 
#  saveRDS(aging_table, "tables/Trainome_aging_table.rds")
 
 
 RT_table <- full_model %>%
   filter(term == "Resistance Training") %>%
   mutate(rank_score = -log10(adj.p) * abs(estimate)) %>%
   arrange(desc(rank_score)) #%>%
#   dplyr::select(gene_intron, effect, estimate, adj.p, gene_label, rank_score)
#  saveRDS(RT_table, "tables/Trainome_RT_table.rds")
 
 
 
 # Functional annotation of genes affecetd by aging and RT
 
 # First load the gene expression dataset
 gene_exp_df <- readRDS("data_new/gene_counts/batch_corrected_genecounts.RDS") 
 
 # select the top 6 most age-affected introns
 top6_introns <- aging_table %>% slice_head(n = 6)
 
 top6_genes <- unique(top6_introns$gene_label)
 
 
 age_exp_df <- gene_exp_df %>%
   filter(gene_name %in% top6_genes) %>%
   pivot_longer(
     cols = -gene_name,
     names_to = "seq_sample_id",
     values_to = "gene_count"
   ) %>%
   inner_join(metadata, by = "seq_sample_id") %>%
   inner_join(top6_introns, by = c("gene_name" = "external_gene_name")) %>%
   group_by(scaled_age, gene_name, estimate) %>%
   summarise(mean_count = mean(gene_count, na.rm = TRUE), .groups = "drop")
 
 
 
 
 ggplot( age_exp_df, aes(x = mean_count, y = estimate)) +
   geom_point(alpha = 0.7, size = 2) +
   geom_smooth(method = "lm", se = FALSE) +
   scale_x_log10() +
   theme_minimal() +
   labs(x = "Mean count (log10)", y = "Effective size")
 
 
 
 age_introns_df <- all_splice_df %>%
   dplyr::filter(transcript_ID %in% top6_introns$target) %>%
   pivot_longer(names_to = "seq_sample_id",
                values_to = "SE",
                cols = -(transcript_ID) ) %>%
   inner_join(top6_introns, by = c("transcript_ID" = "target")) %>%
   inner_join(metadata, by = "seq_sample_id")%>% 
   group_by(scaled_age, gene_intron, transcript_ID) %>%
   summarise(mean_SE = mean(SE, na.rm = TRUE), .groups = "drop")
 
 # Create SE plots, one per intron

 

 
 intron_plots <- lapply(unique(age_introns_df$gene_intron), function(intron) {
   ggplot(
     age_introns_df %>% filter(gene_intron == intron),
     aes(x = scaled_age, y = mean_SE)
   ) +
     geom_point(size = 1.5, alpha = 0.9, colour = "black") +
     geom_smooth(method = "lm", se = FALSE, colour = "red", size = 0.6) +
     labs(title = intron, x = NULL, y = NULL) +
     theme_minimal() +
     theme(
       plot.title = element_text(size = 7, face = "bold"),
       axis.text = element_text(size = 6),
       axis.title = element_blank()
     )
 })
 
 
 # stack them into a panel

 
 mini_panel <- plot_grid(
   plotlist = intron_plots,
   ncol = 3
 )
 
 
 volcano_plot <- ggplot((full_model_outputs %>%
                           filter(term == "Aging")), aes(estimate, neg_log10_fdr, colour = effect)) +
   geom_point(aes(colour = effect), alpha = 0.7, size = 2) +
   
   geom_text_repel(
     data = top6_introns,
     aes(label = gene_intron),
     size = 4,
     max.overlaps = Inf
   ) +
   
   
   geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
   
   scale_color_manual(values = c("Improved SE" =  colors[6], "Reduced SE" = colors[1]),
                      name = "Effect") +
   
   # facet_wrap(~term, scales = "free") +
   # coord_cartesian(xlim = c(-1, 1))+
   
   labs(
     title = "Differentially spliced introns due to Aging",
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
     y = 0.35,  # adjust vertically
     width = 0.38,
     height = 0.48
   )
 
 final_volcano
 
 
 
 
 
 # select the top 6 most age-affected introns
 top6_RT_introns <- RT_table %>% slice_head(n = 6)
 
 top6_RT_genes <- unique(top6_RT_introns$gene_label)
 
 
 
 RT_introns_df <- all_splice_df %>%
   dplyr::filter(transcript_ID %in% top6_RT_introns$target) %>%
   pivot_longer(names_to = "seq_sample_id",
                values_to = "SE",
                cols = -(transcript_ID) ) %>%
   inner_join(top6_RT_introns, by = c("transcript_ID" = "target")) %>%
   inner_join(metadata, by = "seq_sample_id")%>% 
   group_by( time, gene_intron, transcript_ID) %>%
   summarise(mean_SE = mean(SE, na.rm = TRUE), .groups = "drop")
 
 # Create SE plots, one per intron
 
 
 
 
 RT_intron_plots <- lapply(unique(RT_introns_df$gene_intron), function(intron) {
   ggplot(
     RT_introns_df %>% filter(gene_intron == intron),
     aes(x = time, y = mean_SE, colour = time, e)
   ) +
     geom_point(size = 1.5, alpha = 0.9) +
     geom_smooth(method = "lm", se = FALSE, size = 0.6) +
     labs(title = intron, x = NULL, y = NULL) +
     theme_minimal() +
     theme(
       plot.title = element_text(size = 7, face = "bold"),
       axis.text = element_text(size = 6),
       axis.title = element_blank()
     )
 })
 
 
 # stack them into a panel
 
 
 RT_mini_panel <- plot_grid(
   plotlist = RT_intron_plots,
   ncol = 3
 )
 
 
 RT_volcano_plot <- ggplot((full_model_outputs %>%
                           filter(term == "Resistance Training")), aes(estimate, neg_log10_fdr, colour = effect)) +
   geom_point(aes(colour = effect), alpha = 0.7, size = 2) +
   
   geom_text_repel(
     data = top6_RT_introns,
     aes(label = gene_intron),
     size = 4,
     max.overlaps = Inf
   ) +
   
   
   geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
   
   scale_color_manual(values = c("Improved SE" =  colors[6], "Reduced SE" = colors[1]),
                      name = "Effect") +
   
   # facet_wrap(~term, scales = "free") +
   # coord_cartesian(xlim = c(-1, 1))+
   
   labs(
     title = "Differentially spliced introns due to Resistance Training",
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
 
 RT_final_volcano <- ggdraw() +
   draw_plot(RT_volcano_plot) +
   draw_plot(
     RT_mini_panel,
     x = 0.55,   # adjust horizontally
     y = 0.55,  # adjust vertically
     width = 0.38,
     height = 0.4
   )
 
 RT_final_volcano
 
 
 
 
 
 # Repeat for Reistance training
 
 
 full_ds_RT <- full_model %>%
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
 
 
 
 
# ggsave("Figures/GO_RT.png", bg = colors[4], scale=2.5, dpi = 400)
 
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
 
# ggsave("Figures/GO_aging.png", bg = colors[4],height = 20, width = 15, dpi = 400)
 
 
 
 
 Fig4 <-  go_RT + go_aging 
 Fig4 +
   plot_annotation(tag_levels = "A") &
   #plot_layout(widths = c(1.1, 1))
   theme(
     plot.tag = element_text(size = 14, face = "bold")
     ,
     plot.tag.position = c(0.08, 0.98)
   ) 
 
 
# ggsave("Figures/Figure_4.png", bg = colors[4], scale = 4, dpi = 400)
 
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
 
