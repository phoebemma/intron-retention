
library(dplyr)
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(effects)



# Load the postexercise model summary

RT_model_summary <- readRDS("data_new/model_summary_RT_model.RDS") %>%
  filter(coef != "(Intercept)")

# load the predictions 

RT_predictions <- readRDS("data_new/predictions_RT_model.RDS")



RT_merged <- RT_predictions %>%
  inner_join(RT_model_summary, by = "target") %>%
  drop_na() %>%
  mutate(effect = ifelse(p.val <= 0.05,  "effect", "no_effect"))

RT_merged %>%
  ggplot(aes(x = scaled_age, y = round(fit, 2), colour = effect))+
  geom_point() +
  labs(title = "Interaction Effects Plot",
       x = "Predictor",
       y = "Response") +
  theme_minimal()



RT_alone <- RT_merged %>%
  filter(coef == "timePostExc" & fcthreshold == "s") 
RT_alone %>%
  filter(effect == "effect") %>%
  group_by(target) %>%
  ggplot(aes(x = scaled_age, y = round(fit, 2), colour = time))+
  geom_point() +
  labs(title = "RT Effects Plot",
       x = "scaled age",
       y = "splicing efficiency") +
  theme_minimal()




RT_with_age  <- RT_merged %>%
  filter(coef == "scaled_age:timePostExc" & fcthreshold == "s")

RT_with_age %>%
 # filter(effect == "effect" ) %>%
#  group_by(target) %>%
  ggplot(aes(x = scaled_age, y = round(fit, 2), colour = time))+
  geom_point() +
  labs(title = "Interaction Effects RT with age",
       x = "scaled age",
       y = "splicing efficiency") +
  theme_minimal()


RT_scaled_age <- RT_merged %>%
  filter(coef == "scaled_age" & fcthreshold == "s")

RT_scaled_age %>%
  ggplot(aes(x = scaled_age, y = fit, colour = time, group = time))+
  geom_point() +
  labs(title = "Interaction Effects Plot",
       x = "scaled age",
       y = "splicing efficiency") +
  theme_minimal()


ggplot(RT_merged, aes(x = estimate, y = odds_ratio, colour = time)) +
  geom_point()+
  facet_wrap(~coef) 
