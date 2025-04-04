library(dplyr)
library(tidyverse)

# Load baseline model that uses scaled age
scaled_age <- readRDS("data_new/models/new_edition/scaled_age_baseline.RDS") %>%
  filter(p.val<= 0.05  ) %>%
  subset(coef != "(Intercept)") %>%
  drop_na()
  
 
#  preds <- c(coef(m)[1] + coef(m)[2] * seq(from = 0, to = 1, by = 0.1))
length(unique(scaled_age$target))

ggplot(scaled_age, aes(x = estimate, y = odds_ratio)) +
  geom_point()+
  facet_wrap(~coef)+
geom_ribbon(aes(ymin = exp(cil), ymax = exp(ciu)),fill = "grey", alpha = 0.5) 


grouped_age <- readRDS("data_new/models/new_edition/grouped_age_baseline.RDS") %>%
filter(p.val<= 0.05  ) %>%
  subset(coef != "(Intercept)") %>%
  drop_na()
ggplot(grouped_age, aes(x = estimate, y = odds_ratio)) +
  geom_point()+
  facet_wrap(~coef)+
  geom_ribbon(aes(ymin = exp(cil), ymax = exp(ciu)),fill = "grey", alpha = 0.5) 



# RT effect
RT_scaled_age <- readRDS("data_new/models/new_edition/RT_full_scaledage_model.RDS")  %>%
  filter(p.val<= 0.05  ) %>%
  subset(coef != "(Intercept)") %>%
  drop_na()

ggplot(RT_scaled_age, aes(x = estimate, y = odds_ratio)) +
  geom_point()+
  facet_wrap(~coef)+
  geom_ribbon(aes(ymin = exp(cil), ymax = exp(ciu)),fill = "grey", alpha = 0.5) 



RT_grouped_age <- readRDS("data_new/models/new_edition/RT_full_groupedage_model.RDS")%>%
  filter(p.val<= 0.05  ) %>%
  subset(coef != "(Intercept)") %>%
  drop_na()

ggplot(RT_grouped_age, aes(x = estimate, y = odds_ratio)) +
  geom_point()+
  facet_wrap(~coef)+
  geom_ribbon(aes(ymin = exp(cil), ymax = exp(ciu)),fill = "grey", alpha = 0.5) +
  ylim(0,7)


batch_corrected_counts <- readRDS("data_new/gene_counts/batch_corrected_genecounts.RDS") %>%
  mutate(across(where(is.numeric), round, digits = 2))

# Load the gene annotation file
gene_annotation <- readRDS("data_new/ensembl_gene_annotation.RDS")

pre_scaled_sep <- scaled_age %>%
  # select only those ds by age
  filter (coef == "scaled_age") %>%
  separate("target", c("transcript_ID", "intron_ID", "chr"), sep = "_") %>%
  inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version")) %>%
  dplyr::select(estimate, transcript_ID,intron_ID,  transcript_biotype, ensembl_gene_id,ensembl_gene_id_version, external_gene_name, transcript_length)%>%
  mutate(intron_ID = paste0(transcript_ID, "_", intron_ID),
         group = if_else(estimate > 0, "improved", "reduced"))


pre_counts <- RT <- batch_corrected_counts[batch_corrected_counts$gene_id %in% pre_scaled_sep$ensembl_gene_id_version,]


pre_df  <- pre_counts %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "count",
               cols = -(gene_id))%>%
  group_by(gene_id)

merged <- merge(pre_df, pre_scaled_sep , by.x = "gene_id" , by.y= "ensembl_gene_id_version") %>%
  group_by(gene_id, group) %>%
  summarize(mean_count = mean(count))

ggplot(merged, aes(y = mean_count, x = group ))+
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Gene Expression Across Conditions",
       x = "Condition",
       y = "Mean Expression") +
  theme_minimal()


x <- pre_scaled_sep %>%
                filter(group == "reduced")

y <- pre_scaled_sep%>%
                filter(group == "improved")

ggplot(pre_scaled_sep, aes(group))+
  geom_bar()

# View RT and the scaled age

RT_scaled_sep <- RT_scaled_age %>%
  separate("target", c("transcript_ID", "intron_ID", "chr"), sep = "_") %>%
  inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version")) %>%
  dplyr::select(coef, estimate, transcript_ID,intron_ID,  transcript_biotype, ensembl_gene_id,ensembl_gene_id_version, external_gene_name, transcript_length)%>%
  mutate(intron_ID = paste0(transcript_ID, "_", intron_ID),
         group = if_else(estimate > 0, "improved", "reduced"))

ggplot(RT_scaled_sep, aes(group))+
  geom_bar()+
  facet_wrap(~ coef)


# View RT and the age groups
RT_grouped_sep <- RT_grouped_age %>%
  separate("target", c("transcript_ID", "intron_ID", "chr"), sep = "_") %>%
  inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version")) %>%
  dplyr::select(coef, estimate, transcript_ID,intron_ID,  transcript_biotype, ensembl_gene_id,ensembl_gene_id_version, external_gene_name, transcript_length)%>%
  mutate(intron_ID = paste0(transcript_ID, "_", intron_ID),
         group = if_else(estimate > 0, "improved", "reduced"))

ggplot(RT_grouped_sep, aes(group))+
  geom_bar()+
  facet_wrap(~ coef)



