#Load the file with libraries
source("R/libraries.R")

#Load the trainome functions file
source("R/Trainome_functions.R")


#Load the full splicing_data
all_splice_df <- readRDS("./data/Full_splice_data.RDS")

#Load metadata

all_metadata <- readr::read_csv("./data/all_metadata.csv")%>%
  #drop the "MidExc rows 
  subset(time != "MidExc") %>%
  mutate(time = factor(time, levels = c("PreExc",  "PostExc")), 
         age = factor(age, levels = c("young", "old")),
         sex = factor (sex, levels = c("male", "female")))

#check the distribution of old versus young in the metadata
plot(all_metadata$age)

#distribution time
plot(all_metadata$time)

#plot distrribution study
plot(table(all_metadata$study))


#get the sample_ids that exist in both splicing data and metadata
splice_intersect <- (intersect(colnames(all_splice_df),
                               all_metadata$sample_id))

#subset the splicing data to only include the intersects
splice_data <- all_splice_df %>%
  subset( select = c("transcript_ID", splice_intersect))
colnames(splice_data)



#select only metadata that intersect
meta_df <- all_metadata %>%
  dplyr::filter(sample_id %in% c(splice_intersect)) 


plot(meta_df$age, main = "Distribution of age data")

plot(meta_df$time, main = "Distribution of sampling time")


#convert the splicing data to a longer form and merge with metadat
full_splice_df <- all_splice_df %>%
  pivot_longer(names_to = "sample_id",
               values_to = "SE",
               cols = -(transcript_ID) )

full_splice_df <-  merge(meta_df, full_splice_df, by= "sample_id")

#save the data as RDS
#saveRDS(full_splice_df, "./data/all_splice_data_with_metadata_long_form.RDS")

#Subset to the Pre-exercise data

Pre_exc_df <- full_splice_df %>%
  subset(time == "PreExc")


#Subset the post exercise data

Post_ex_df <- full_splice_df %>%
  subset(time == "PostExc")


int_both <- (intersect(Pre_exc_df$participant,
                                           Post_ex_df$participant))

#subset the splicing data to only include the intersects

Pre_Exc_in <- Pre_exc_df[Pre_exc_df$participant %in% int_both,]

Post_exc_int <- Post_ex_df[Post_ex_df$participant %in% int_both,]


splice_df_full <- rbind(Pre_Exc_in, Post_exc_int)




plot(splice_df_full$age, main = "Distribution of datapoints across age groups")


#Plot the full dataframe

splice_df_full %>%  dplyr::mutate(n = n(),
                                          bins = cut(SE, 100, labels = FALSE)) %>%
  group_by(age, time, transcript_ID, bins) %>%
  dplyr::summarise(freq = n() / mean(n)) %>%
  #dplyr::filter(bins < 50)%>%
  ggplot(aes(bins, freq, fill = time)) +
  geom_col(position = "dodge",  size = 8, width = 1.5) + 
  facet_wrap(~age)+
  ggtitle("Splicing efficiency old versus young")+
  theme(axis.text = element_text(size = 15), text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))





#Zoom into the data above 0.5 SE
splice_df_full %>%  dplyr::mutate(n = n(),
                                  bins = cut(SE, 100, labels = FALSE)) %>%
  group_by(age, time,transcript_ID, bins) %>%
  dplyr::summarise(freq = n() / mean(n)) %>%
  dplyr::filter(bins > 50)%>%
  ggplot(aes(bins, freq, fill = time)) +
  geom_col(position = "dodge",  size = 8, width = 1.5) + 
  facet_wrap(~age)+
  ggtitle("Splicing efficiency old versus young above 0.5")+
  theme(axis.text = element_text(size = 15), text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))


#Zoom into the data above 0.9 SE
splice_df_full %>%  dplyr::mutate(n = n(),
                                  bins = cut(SE, 100, labels = FALSE)) %>%
  group_by(age, time,transcript_ID, bins) %>%
  dplyr::summarise(freq = n() / mean(n)) %>%
  dplyr::filter(bins > 80)%>%
  ggplot(aes(bins, freq, fill = time)) +
  geom_col(position = "dodge",  size = 8, width = 1.5) + 
  facet_wrap(~age)+
  ggtitle("Splicing efficiency old versus young above 0.8")+
  theme(axis.text = element_text(size = 15), text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))



#Zoom into the data below 0.2 SE
splice_df_full %>%  dplyr::mutate(n = n(),
                                  bins = cut(SE, 100, labels = FALSE)) %>%
  group_by(age, time, transcript_ID,bins) %>%
  dplyr::summarise(freq = n() / mean(n)) %>%
  dplyr::filter(bins < 20)%>%
  ggplot(aes(bins, freq, fill = time)) +
  geom_col(position = "dodge",  size = 8, width = 1.5) + 
  facet_wrap(~age)+
  ggtitle("Splicing efficiency old versus young below 0.2")+
  theme(axis.text = element_text(size = 15), text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))





#Plot the full dataframe
#Line plot

# splice_df_full %>%  dplyr::mutate(n = n(),
#                                   bins = cut(SE, 100, labels = FALSE)) %>%
#   group_by(age, time, transcript_ID, bins) %>%
#   dplyr::summarise(freq = n() / mean(n)) %>%
#   #dplyr::filter(bins < 50)%>%
#   ggplot(aes(bins, freq, fill = time)) +
#   #geom_col(position = "dodge",  size = 8, width = 1.5) + 
#   geom_line()+
#   geom_point()+
#   facet_wrap(~age)+
#   ggtitle("Splicing efficiency old versus young")+
#   theme(axis.text = element_text(size = 15), text = element_text(size = 15),
#         plot.title = element_text(hjust = 0.5))
# 
# 



### Repeat the intersect analyses for introns. That is, introns recorded in both pre
#and post exercise


introns_both <- (intersect(Pre_exc_df$transcript_ID,
                       Post_ex_df$transcript_ID))

Pre_Exc_introns <- Pre_exc_df[Pre_exc_df$transcript_ID %in% introns_both,]

Post_exc_introns <- Post_ex_df[Post_ex_df$transcript_ID %in% introns_both,]

#How many unique introns do we have i the pre exercise data
length(unique(Pre_Exc_introns$transcript_ID))


#How many in the post exercise data
length(unique(Post_exc_introns$transcript_ID))


#How many unique participants in pre exercise data

length(unique(Pre_Exc_introns$participant))


#How many unique participants in the pre exercise filtered based on intercepting participants with old
length(unique(Pre_Exc_in$participant))


#How many unique participants in post exercise data
length(unique(Post_exc_introns$participant))

length(unique(Post_exc_int$participant))


# Plot the postexercise data

#Plot the full dataframe

Post_exc_introns %>%  dplyr::mutate(n = n(),
                                  bins = cut(SE, 100, labels = FALSE)) %>%
  group_by(age,  transcript_ID, bins) %>%
  dplyr::summarise(freq = n() / mean(n)) %>%
  #dplyr::filter(bins < 50)%>%
  ggplot(aes(bins, freq, fill = age)) +
  geom_col(position = "dodge",  size = 8, width = 1.5) + 
  #facet_wrap(~age)+
  ggtitle("Splicing efficiency old versus young post exercise")+
  theme(axis.text = element_text(size = 15), text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))


Pre_Exc_introns %>%  dplyr::mutate(n = n(),
                                bins = cut(SE, 100, labels = FALSE)) %>%
  group_by(age,  transcript_ID, bins) %>%
  dplyr::summarise(freq = n() / mean(n)) %>%
  #dplyr::filter(bins < 50)%>%
  ggplot(aes(bins, freq, fill = age)) +
  geom_col(position = "dodge",  size = 8, width = 1.5) + 
  #facet_wrap(~age)+
  ggtitle("Splicing efficiency old versus young pre exercise")+
  theme(axis.text = element_text(size = 15), text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))
