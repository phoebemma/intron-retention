#Load the file with libraries
source("R/libraries.R")

#Load the trainome functions file
source("R/Trainome_functions.R")


#Load the long form data
#This contains the splicing data merged with the metadata
full_splice_df <- readRDS("./data/all_splice_data_with_metadata_long_form.RDS")


#Subset to the Pre-exercise data

#Pre_exc_df <- full_splice_df %>%
#  subset(time == "PreExc")


#Subset the post exercise data

# Post_ex_df <- full_splice_df %>%
#   subset(time == "PostExc")



#Subset only the introns wih SE 0

fully_retained_ints <- full_splice_df %>%
  subset(SE == 0)


table(fully_retained_ints$age)
table(fully_retained_ints$time)


#Plot the data
plot(fully_retained_ints$time)
