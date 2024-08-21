#Exploratory analyses of the data

library(dplyr)
library(tidyverse)
library(gridExtra)
library(ggpubr)


#Explore the splicing data

splice_df <- readRDS("data/model/full_splice_data.RDS") 

#load the metadata
metadf <- readRDS("data/model/full_metadata.RDS")

   
long_df <- splice_df %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID) )%>%
  inner_join(metadf, by = "seq_sample_id") %>%
  #create a new column that states whether intron is retained or not
  mutate(status = if_else(SE == 1, "none_retained", "retained"),
         reverse_status = case_when(SE == 0.0 ~ "completely_retained", 
                                    SE > 0.0 $ SE < 1.0 ~ "partially_retained", 
                                    SE == 1.0 ~ "completely_spliced"))

hist(long_df$SE)

unique(long_df$reverse_status)

long_prExc_df <- long_df %>%
  subset(time == "PreExc")
long_postExc <- long_df %>%
  subset(time == "PostExc")

old_pre_df <- long_prExc_df %>%
  subset(age_group == "Old") %>%
  group_by(reverse_status)%>%
  count()%>% #count the occurence
  ungroup()%>%
  mutate(perc = n/sum(n))%>% #extract the percentage
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))


young_pre_df <- long_prExc_df %>%
  subset(age_group == "Young") %>%
  group_by(reverse_status)%>%
  count()%>% #count the occurence
  ungroup()%>%
  mutate(perc = n/sum(n))%>% #extract the percentage
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

old_post_df <- long_postExc %>%
  subset(age_group == "Old") %>%
  group_by(reverse_status)%>%
  count()%>% #count the occurence
  ungroup()%>%
  mutate(perc = n/sum(n))%>% #extract the percentage
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

young_post_df <- long_postExc %>%
  subset(age_group == "Young") %>%
  group_by(reverse_status)%>%
  count()%>% #count the occurence
  ungroup()%>%
  mutate(perc = n/sum(n))%>% #extract the percentage
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

old_pre <- ggplot(old_pre_df, aes(x = "", y = perc, fill = reverse_status)) +
  geom_col()+
  geom_text(aes(label = labels, x = 1.6),
            position = position_stack(vjust = 0.8)) +
  coord_polar(theta = "y")+
  theme_void()+
  ggtitle("percentage of  introns retained in old participants at baseline")

young_pre <- ggplot(young_pre_df, aes(x = "", y = perc, fill = reverse_status)) +
  geom_col()+
  geom_text(aes(label = labels, x = 1.6),
            position = position_stack(vjust = 0.8)) +
  coord_polar(theta = "y")+
  theme_void()+
  ggtitle("percentage of  introns retained in young participants at baseline")

old_post <-  ggplot(old_post_df, aes(x = "", y = perc, fill = reverse_status)) +
  geom_col()+
  geom_text(aes(label = labels, x = 1.6),
            position = position_stack(vjust = 0.8)) +
  coord_polar(theta = "y")+
  theme_void()+
  ggtitle("percentage of  introns retained in old participants at postexercise")


young_post <- ggplot(young_post_df, aes(x = "", y = perc, fill = reverse_status)) +
  geom_col()+
  geom_text(aes(label = labels, x = 1.6),
            position = position_stack(vjust = 0.8)) +
  coord_polar(theta = "y")+
  theme_void()+
  ggtitle("percentage of  introns retained in young participants at postexercise")


ggarrange(old_pre,old_post, young_pre,young_post, #you can also specify rremove("x.text) to remove texts on the x axis 
         labels = c("A", "B", "C", "D"), #Add labels like done in publications
         nrow = 2, ncol = 2) 
