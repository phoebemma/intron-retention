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
