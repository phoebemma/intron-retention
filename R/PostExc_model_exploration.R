
library(dplyr)
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(hexbin)
library(ggridges)

# Load the postexercise model summary

RT_model_summary <- readRDS("data_new/simpler_RT_model_summary.RDS") %>%
  filter(coef != "(Intercept)")

# load the predictions 

RT_predictions <- readRDS("data_new/predictions_simpler_RT_model.RDS")
colnames(RT_model_summary)

colnames(RT_merged)
RT_merged <- RT_predictions %>%
  inner_join(RT_model_summary, by = "target") %>%
  drop_na() %>%
  mutate(effect = case_when(Estimate > 0 & Pr...z.. <= 0.05 ~ "improved_SE",
                            Estimate < 0 & Pr...z.. <= 0.05 ~ "reduced_SE" ,
                            Estimate < 0 & Pr...z.. > 0.05 ~ "no effect",
                            Estimate > 0 & Pr...z.. > 0.05 ~ "no effect"))
# 
# RT_merged %>%
#   ggplot(aes(x = scaled_age, y = round(fit, 2), colour = effect))+
#   geom_point() +
#   labs(title = "Interaction Effects Plot",
#        x = "Predictor",
#        y = "Response") +
#   theme_minimal()



RT_alone <- RT_merged %>%
  filter(coef == "timePostExc") 
RT_alone %>%
  group_by(target) %>%
  ggplot(aes(y= Estimate, x = odds_ratio, colour = time))+
  geom_point() +
  labs(title = "RT Effects Plot",
       x = "scaled age",
       y = "splicing efficiency") +
#  theme_minimal() +
  facet_wrap(~effect)



ggplot(RT_alone , aes(x = Estimate, y = odds_ratio, colour = time)) +
  geom_point(size = 3, alpha = 0.5) +
  theme_minimal() +
  labs(title = "Relationship between Fit and Scaled Age", x = "Scaled Age", y = "Fit") +
  scale_color_manual(values = c("PreExc" = "blue", "PostExc" = "red"))   # Custom color mapping
#  geom_smooth(method = "lm", col = "red")  # Adding a regression line

ggplot(RT_alone, aes(x = scaled_age, y = round(fit, 2))) +
  geom_smooth() +
  theme_minimal() +
  labs(title = "Smooth Plot of Fit vs. Scaled Age", x = "Scaled Age", y = "Fit")+
  facet_wrap(~effect)

# bin <- RT_alone[1:20,]
# 
# write.csv(bin, "x.csv")

RT_with_age  <- RT_merged %>%
  filter(coef == "scaled_age:timePostExc")

ggplot(RT_with_age, aes(x = scaled_age, y = round(fit, 2), colour = time)) +
  geom_smooth() +
  theme_minimal() +
  labs(title = "Smooth Plot of Fit vs. Scaled Age", x = "Scaled Age", y = "Fit")+
  facet_wrap(~effect)

RT_with_age %>%
 # filter(effect == "effect" ) %>%
#  group_by(target) %>%
  ggplot(aes(x = scaled_age, y = fit, colour = time, group=paste(target,time)))+
  geom_line() +
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
