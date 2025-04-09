library(dplyr)
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(ggridges)
library(UpSetR)
library(ggplot2)
library(ggplotify)
library(clusterProfiler)
library(org.Hs.eg.db)

# Load the postexercise model summary

RT_model_summary <- readRDS("data_new/simpler_RT_model_summary.RDS") %>%
  # select the model results of those affected by RT alone or interaction with aging
  filter(coef == "timePostExc" | coef == "scaled_age:timePostExc")%>%
  # add the odds ratio for each Estimate
  mutate(odds_ratio = exp(Estimate),
         type = "model coefficient")

# load the predictions 

RT_predictions <- readRDS("data_new/predictions_simpler_RT_model.RDS")%>%
  # Extract the transcript_id from the target
  mutate(transcript_ID = str_split(target, "_",simplify= T) [,1])


# Load the gene annotation file
gene_annotation <- readRDS("data_new/ensembl_gene_annotation.RDS")


# merge dataset
RT_merged <- RT_predictions %>%
  inner_join(RT_model_summary, by = "target") %>%
  drop_na() %>%
  mutate(effect = case_when(Estimate > 0 & Pr...z.. <= 0.05 ~ "Improved SE",
                            Estimate < 0 & Pr...z.. <= 0.05 ~ "Reduced SE" ,
                            Estimate < 0 & Pr...z.. > 0.05 ~ "No effect",
                            Estimate > 0 & Pr...z.. > 0.05 ~ "No effect")) %>%
  inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version"))

saveRDS(RT_merged, "data_new/RT_model_outputs/RT_model_df.RDS")


# select only those affected by RT
RT_effect <- RT_merged %>%
  filter(effect != "No effect" )



# Extract the number of the various effect groups
summary_df <- RT_effect %>%
  group_by(effect) %>%
  summarize(num_targets = n_distinct(target))


ggplot(RT_effect, aes(x = scaled_age, y = fit,  group = target, colour = time)) +
  geom_line(aes(alpha = 0.5, colour = "grey"), show.legend = F) + 
  theme_minimal()+
  labs(title = "Relationship between age and splicing efficiency", x = "Scaled Age", y = "Splicing efficiency") +
  facet_wrap(~effect) +
  scale_alpha_identity()+
  scale_color_manual(values = c("grey"), guide = "none")+
  geom_text(data = summary_df, aes(x = Inf, y = Inf, label = paste("Number of introns:", num_targets)), 
            hjust = 1.1, vjust = 1.1, size = 3, color = "black", inherit.aes = FALSE) +
  theme(plot.title = element_text(hjust = 0.5))



# 
# RT_merged %>%
#   ggplot(aes(x = scaled_age, y = round(fit, 2), colour = effect))+
#   geom_point() +
#   labs(title = "Interaction Effects Plot",
#        x = "Predictor",
#        y = "Response") +
#   theme_minimal()



RT_alone <- RT_effect %>%
  filter(coef == "timePostExc") 
# RT_alone %>%
#   group_by(target) %>%
#   ggplot(aes(y= Estimate, x = odds_ratio, colour = time, group = target))+
#   geom_line() +
#   labs(title = "RT Effects Plot",
#        x = "scaled age",
#        y = "splicing efficiency") +
# #  theme_minimal() +
#   facet_wrap(~effect)


# 
# ggplot(RT_alone , aes(x = Estimate, y = odds_ratio, colour = time)) +
#   geom_point(size = 3, alpha = 0.5) +
#   theme_minimal() +
#   labs(title = "Relationship between Fit and Scaled Age", x = "Scaled Age", y = "Fit") +
#   scale_color_manual(values = c("PreExc" = "gray20", "PostExc" = "gray70"))   # Custom color mapping
#  geom_smooth(method = "lm", col = "red")  # Adding a regression line

ggplot(RT_alone, aes(x = scaled_age, y = fit, group = target)) +
  geom_line(alpha = 0.1) +
  theme_minimal() +
  labs(title = "SE of introns affected by RT alone", x = "Scaled Age", y = "Splicing efficiency")+
  facet_wrap(~time)#+
 # scale_color_manual(values = c("PreExc" = "black", "PostExc" = "gray70"))   # Custom color mapping
#  geom_smooth(method = "lm", col = "red")  # Adding a regression line




# bin <- RT_alone[1:20,]
# 
# write.csv(bin, "x.csv")

RT_with_age  <- RT_effect %>%
  filter(coef == "scaled_age:timePostExc")

ggplot(RT_with_age, aes(x = scaled_age, y = fit, group = target)) +
  geom_line(alpha = 0.1) +
  theme_minimal() +
  labs(title = "SE of introns affected by RT in interaction with aging", x = "Scaled Age", y = "Splicing efficiency")+
  facet_wrap(~time)

# RT_with_age %>%
#  # filter(effect == "effect" ) %>%
# #  group_by(target) %>%
#   ggplot(aes(x = scaled_age, y = fit, colour = time, group=paste(target,time)))+
#   geom_line() +
#   labs(title = "Interaction Effects RT with age",
#        x = "scaled age",
#        y = "splicing efficiency") +
#   theme_minimal()
# 
# 
# RT_scaled_age <- RT_merged %>%
#   filter(coef == "scaled_age" & fcthreshold == "s")
# 
# RT_scaled_age %>%
#   ggplot(aes(x = scaled_age, y = fit, colour = time, group = time))+
#   geom_point() +
#   labs(title = "Interaction Effects Plot",
#        x = "scaled age",
#        y = "splicing efficiency") +
#   theme_minimal()
# 
# 
# ggplot(RT_merged, aes(x = estimate, y = odds_ratio, colour = time)) +
#   geom_point()+
#   facet_wrap(~coef) 
