source("./R/libraries.R")
source("./R/Trainome_functions.R")

#est <- read_Splice_Q("./data/Contratrain_SpliceQ_outputs/1-subj1sample1.tsv")

#The last two underscores in the column "transcript_ID" are the intron_Id and chromosome respectively

#Load contratrain data
contratrain_data <- extract_splice_q("./data/Contratrain_SpliceQ_outputs/")

idx <- sapply(contratrain_data, class)== "numeric"
contratrain_data[, idx] <- lapply(contratrain_data[, idx], round, 2)
#saveRDS(contratrain_data, "./data/contratrain_splicing_data.RDS")



#Load Volume data
volume_data <- extract_splice_q("./data/Volume_SpliceQ_outputs/")
idx <- sapply(volume_data, class)== "numeric"
volume_data[, idx] <- lapply(volume_data[, idx], round, 2)

#remove everything before the . in sample_id. 
colnames(volume_data) <- gsub(".*?\\.", "", colnames(volume_data) )
#saveRDS(volume_data, "./data/volume_splicing_data.RDS")




#Load COPD data
copd_data <- extract_splice_q("./data/COPD_spliceQ_outputs/")
#Rename the ccolumn names by removing everything after the second underscore
 colnames(copd_data) <- gsub("_.*", "", colnames(copd_data) )
 #COPD sequence data contains duplicates. The code below should remove duplicate columns
copd_data <-   copd_data[, !duplicated(colnames(copd_data))]

#the gsub functio above also removed the underscore in "transcript_ID, replace it
colnames(copd_data)[1] <- "transcript_ID"
#remove everything before the . in sample_id. 
colnames(copd_data) <- gsub(".*?\\.", "", colnames(copd_data) )
#Add an X before the sample_ID. This makes it match to the naming in metadata
colnames(copd_data)[-1] <- paste0("X", colnames(copd_data)[-1] )

idx <- sapply(copd_data, class)== "numeric"
copd_data[, idx] <- lapply(copd_data[, idx], round, 2)
#saveRDS(copd_data, "./data/copd_splicing_data.RDS")






#Load the publicly available data


#This contains pre and postexercise (RT) data
#see article https://pubmed.ncbi.nlm.nih.gov/28273480/
SRP102542_data <- extract_splice_q("./data/SRP102542_GSE97084_SpliceQ_outputs/")

idx <- sapply(SRP102542_data, class)== "numeric"
SRP102542_data[, idx] <- lapply(SRP102542_data[, idx], round, 2)
#saveRDS(SRP102542_data, "./data/SRP102542_splicing_data.RDS")



#The following contains only baseline data in young participants
#See here https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58608 
#Participants were young see https://faseb.onlinelibrary.wiley.com/doi/epdf/10.1096/fj.14-255000
SRP043368_data <- extract_splice_q("./data/SRP043368_GSE58608_SpliceQ_outputs/")
idx <- sapply(SRP043368_data, class)== "numeric"
SRP043368_data[, idx] <- lapply(SRP043368_data[, idx], round, 2)

#saveRDS(SRP043368_data, "./data/SRP043368_splicing_data.RDS")




#Contains the baseline and post-14 weeks training of old adults. To be selected only those with placebo and 
#RT
#Also contains only baseline of young individuals

#see article https://pubmed.ncbi.nlm.nih.gov/33071237/
SRP280348_data <- extract_splice_q("./data/SRP280348_GSE157585_SpliceQ_outputs/")
idx <- sapply(SRP280348_data, class)== "numeric"
SRP280348_data[, idx] <- lapply(SRP280348_data[, idx], round, 2)

#saveRDS(SRP280348_data, "./data/SRP280348_splicing_data.RDS")



#Load  the metadata
metadata <- readr::read_csv("./data/all_metadata.csv")  

#Subset to only the pre_exercise data
        #drop the "MidExc rows 
Pre_Exc_metadata <- metadata %>%
        subset(time == "PreExc") %>%
        mutate(age = factor(age, levels = c("young", "old")),
               sex = factor (sex, levels = c("male", "female"))) 

unique(Pre_Exc_metadata$study)
#saveRDS(Pre_Exc_metadata, "./data/PreExc_metadata.RDS")



# Extract sample_ids that are common between splicing and metadata
splice_intersect <- (intersect(colnames(volume_data), Pre_Exc_metadata$sample_id))


#subsett the Volume splicing data to only include the intersects
PreExc_vol_data <- volume_data %>%
        subset( select = c("transcript_ID", splice_intersect))
colnames(PreExc_vol_data)
#reduce the efficiency values to 2 decimal places


#saveRDS(PreExc_vol_data, "./data/PreExc_volume_splicing_data.RDS")





#update splice intersect for the contratrain data
splice_intersect <- (intersect(colnames(contratrain_data), Pre_Exc_metadata$sample_id))


#subsett the Volume splicing data to only include the intersects
PreExc_ct_data <- contratrain_data %>%
        subset( select = c("transcript_ID", splice_intersect))
colnames(PreExc_ct_data)
#saveRDS(PreExc_ct_data, "./data/PreExc_contratrain_splicing_data.RDS")




#update splice intersect for the copd data
splice_intersect <- (intersect(colnames(copd_data), Pre_Exc_metadata$sample_id))


#subsett the Volume splicing data to only include the intersects
PreExc_copd_data <- copd_data %>%
        subset( select = c("transcript_ID", splice_intersect))
colnames(PreExc_copd_data)
#saveRDS(PreExc_copd_data, "./data/PreExc_copd_splicing_data.RDS")





#merge the three dataframes based on their transcript_ID
PreExc_all_data <- PreExc_copd_data %>% 
        left_join(PreExc_vol_data, by= "transcript_ID") %>%
        left_join(PreExc_ct_data, by = "transcript_ID") %>%
        drop_na()
#saveRDS(PreExc_all_data, "./data/all_PreExc_splicing_data.RDS")


#Extracting the full splicing data
Full_splice_data <- copd_data %>%
        left_join(volume_data, by = "transcript_ID") %>%
        left_join(contratrain_data, by = "transcript_ID") %>%
        drop_na()

unique(metadata$time)


#saveRDS(Full_splice_data, "./data/Full_splice_data.RDS")
