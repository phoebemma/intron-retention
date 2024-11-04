library(AOData)
library(dplyr)
library(cowplot)

# Load the participant details
ids <- idkeys %>%
  select(participant, treat, age,  sex )

Sequenced_samples <- seq_samples %>%
  select(participant, condition, time) %>%
  inner_join(ids, by = "participant") %>%
  mutate(age_group = ifelse(age <=40,  "Young","Old"),
         time = case_when(time == "T1" ~ "PreExc",
                          time == "T2" ~ "MidExc",
                          time == "T4" ~ "PostExc"))
Sequenced_samples["participant"] <- as.character(Sequenced_samples$participant)

Sequenced_samples$study <- "Alpha/Omega"
Sequenced_samples$seq_sample_id <- "Alpha/Omega"

Sequenced_samples["sex"][Sequenced_samples["sex"] == "m" ] <- "male"

Sequenced_samples["sex"][Sequenced_samples["sex"] == "f" ] <- "female"
unique(Sequenced_samples$time)

pre <- Sequenced_samples %>%
  filter(time == "PreExc") %>%
  select(study, participant, sex, time, seq_sample_id, age, age_group)
  # Add a dummy sequence ID since its not available yet


saveRDS(pre, "data/preexercise_data/Alpha_Omega_PreExc_metadata.RDS")

range(Sequenced_samples$age)

hist(Sequenced_samples$age)
colnames(ids)
colnames(Sequenced_samples)  
unique(Sequenced_samples$condition)
length(unique(Sequenced_samples$participant))

unique(Sequenced_samples$time)
unique(Sequenced_samples$sex)


# The code below runs only because I have run the preexercise_model.R in the global environment

all_pre_metadata_with_alpha <- rbind(copd_metadata, volume_metadata)%>%
  rbind(Contratrain_metadata) %>%
  
  mutate(age_group = case_when(study == "copd" ~ "Old", 
                               study == "vol" ~ "Young",
                               study == "ct" ~ "Young")) %>%
  rbind(SRP102542_metadata) %>%
  rbind(SRP280348_metadata) %>% 
  rbind(pre) %>%
  
  #Copd and volume and SRP102542 have age data in decimal
  mutate(across(c("age"), round, 0)) %>%
  # Create youps where group 1 = those below 30
  #group 2 is those above 30 but below 51
  # group 3 those above 50 but below 71
  # group 4 is those above 70
  
  mutate(group = case_when(age <=25 ~ "<=25" ,
                           age > 25 & age <= 50 ~ ">25 & <=50", 
                           age > 50 & age <= 70 ~ ">50 & <=70",
                           age > 70 ~ ">70")) %>%
  mutate(sex = factor(sex, levels = c("female", "male")),
         age_group = factor(age_group, levels= c("Young", "Old")),
         group = factor(group, levels = c("<=25", ">25 & <=50", ">50 & <=70", ">70"))) 



 A <- ggplot(all_pre_metadata_with_alpha, aes(age, fill = age_group )) +
  geom_bar()+
  ggtitle("Distribution of baseline data with Alpha/Omega dataset")+
  theme(plot.title = element_text(hjust = 0.5))


B <- ggplot(all_pre_metadata_with_alpha, aes(group, fill = group)) +
  geom_bar()+
 # ggtitle("Distribution of baseline data with Alpha/Omega dataset")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("Age range of participants")+
  stat_count(geom = "Text", aes(label = ..count..), vjust = 1.5)


C <- ggplot(all_pre_metadata_with_alpha, aes(age_group, fill = age_group)) +
  geom_bar()+
 # ggtitle("Distribution of baseline data with Alpha/Omega dataset")+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = 1.5)+
  xlab("Age group of participants")



D <- ggplot(all_pre_metadata_with_alpha, aes(sex, fill = age_group)) +
  geom_bar()+
  #ggtitle("Distribution of baseline data with Alpha/Omega dataset")+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = 1.5)



a <- ggplot(all_pre_metadata, aes(age, fill = age_group )) +
  geom_bar()+
  ggtitle("Distribution of baseline data")+
  theme(plot.title = element_text(hjust = 0.5))


b <- ggplot(all_pre_metadata, aes(group, fill = group)) +
  geom_bar()+
 # ggtitle("Distribution of baseline data")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("Age range of participants")+
  stat_count(geom = "Text", aes(label = ..count..), vjust = 1.5)


c <- ggplot(all_pre_metadata, aes(age_group, fill = age_group)) +
  geom_bar()+
 # ggtitle("Distribution of baseline data")+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = 1.5)+
  xlab("Age group of participants")



d <- ggplot(all_pre_metadata, aes(sex, fill = age_group)) +
  geom_bar()+
  # ggtitle("Distribution of baseline data")+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = 1.5)

plot_grid(a, A, b, B, c, C, d, D, ncol = 2)
