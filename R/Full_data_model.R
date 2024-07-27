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

#get the sample_ids that exist in both splicing data and metadata
splice_intersect <- (intersect(colnames(all_splice_df),
                               all_metadata$sample_id))

#subset the splicing data to only include the intersects
splice_data <- all_splice_df %>%
  subset( select = c("transcript_ID", splice_intersect))

print(colnames(splice_data))



#select only metadata that intersect
meta_df <- all_metadata %>%
  dplyr::filter(sample_id %in% c(splice_intersect)) 

#invert the data to make it suitable for zero inflated analyses
splice_df_inverted <- splice_data %>%
  mutate(across(X102PostExcVLL14:X98.subj40sample7 , function(x)1-x))

#Relace all the 1s in the dataframe to 0.99

splice_df_inverted[splice_df_inverted == 1 ] <- 0.995


#check for missing values
any(is.na(splice_df_inverted))
any(is.na(meta_df))


set.seed(12345)
#Model building

args <- list(formula = y ~  time + age + age:time + sex  + (1|participant),
             ziformula = ~1,
             family = glmmTMB::beta_family()
)



### These data are problematic for the beta-family distribution
temp_df <- splice_df_inverted %>%
  dplyr::filter(transcript_ID == "ENST00000233190.11_3_2") %>%
  pivot_longer(-transcript_ID, names_to = "sample_id") %>%
  
  mutate(zero_one = if_else(value == 0, 1, 0)) %>%
  
  
  left_join(meta_df) %>%
  print()
  
  
  
temp_df %>%  
  ggplot(aes(value)) + geom_histogram() + facet_grid(study ~ sex)
  

## Attempt at modelling 

m1 <- glmmTMB(zero_one ~ time * age + (1|participant), 
              data = temp_df, 
           
              family = betabinomial)

summary(m1)

sim <- DHARMa::simulateResiduals(m1, n = 100)

DHARMa::testDispersion(sim, plot = TRUE)
DHARMa::testUniformity(sim, plot = TRUE)



fits <- seq_wrapper(    fitting_fun = glmmTMB::glmmTMB,
                        arguments = args,
                        data = splice_df_inverted,
                        metadata = meta_df,
                        samplename = "sample_id",
                        summary_fun = sum_fun,
                        eval_fun = eval_mod,
                        additional_vars = NULL,
                        exported = list(),
                        subset = 200:201,
                        cores = ncores)

#check the summary for the first model
summary(fits$model_fits[[7]])


names(fits$model_fits)


sim <- DHARMa::simulateResiduals(fits$model_fits[[1]], n = 1000)

disp <- DHARMa::testDispersion(sim, plot = TRUE)
unif <- DHARMa::testUniformity(sim, plot = TRUE)




bind_rows(fits$model_evaluations) %>%
  mutate(target = names(fits$model_evaluations)) %>%
  #print()
  
  #dplyr::filter(pval.unif < 0.05)%>%
  
  
  ggplot(aes(pval.unif)) + geom_histogram()
  print()




bind_rows(fits$model_summarises) %>%
  mutate(target = rep(names(fits$model_summarises), each = 5)) %>%
  dplyr::filter(coef == "ageold") %>%
  ggplot(aes(Pr...z..)) + geom_histogram(bins=80)+
  ggtitle("Splicing efficiency; impact of age") +
  theme(axis.text = element_text(size = 15), text = element_text(size = 15))

