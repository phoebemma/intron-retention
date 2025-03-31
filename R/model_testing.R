library(marginaleffects)


all_splice_df <- readRDS("data_new/processed_data/all_splice_data.RDS")
all_full_metadata <- readRDS("data_new/processed_data/all_full_metadata.RDS")

long_df <- all_splice_df %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "SE",
               cols = -(transcript_ID) )%>%
  inner_join(all_full_metadata, by = "seq_sample_id")

ggplot(long_df, aes(z = SE, y = age, x = as.numeric(time)))+
  geom_contour_filled() +
  labs(title = "Contour Plot of Interaction Effects", x = "time", y = "age")


# Load the Trainome functions
# Contains functions needed for model building
source("R/Trainome_functions.R")


# Load metadata
all_pre_metadata <- readRDS("data_new/Pre_Exercise/all_prexercise_metadata.RDS")

# Standardize the age by scaling them 0 to 1
all_pre_metadata$scaled_age <- round(rescale(all_pre_metadata$age), digits = 2)



# Load Splicing data
all_pre_splice <- readRDS("data_new/Pre_Exercise/all_pre_Exc_splicing_data.RDS")



# select only samples present in the metadata
all_pre_metadata <- all_pre_metadata %>%
  filter((seq_sample_id %in% colnames(all_pre_splice[,-1]))) 


# reorder the sample ids to match how they occur in the metadata
all_pre_splice_reordered <- all_pre_splice[ , c("transcript_ID",all_pre_metadata$seq_sample_id)]


# Check if everything matches except the transcript_id
match(colnames(all_pre_splice_reordered), all_pre_metadata$seq_sample_id)



# convert the 1.0 to 0.999. This is because beta-model accepts only values between 0 and one
all_pre_splice_reordered[all_pre_splice_reordered == 1 ] <- 0.999

# This argument would estimate the intercept, and the slope separately
# with uncorrelated random intercept and random slope within each study
arg_1<- list(formula = y ~  scaled_age + sex + (1|study) + (scaled_age+0|study) +(1|participant), 
             family = glmmTMB::beta_family())


model_1 <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                   arguments = arg_1,
                   data = all_pre_splice_reordered,
                   metadata = all_pre_metadata,
                   samplename = "seq_sample_id",
                   summary_fun = sum_fun2,
                   eval_fun = eval_mod,
                   exported = list(),
                   save_models = FALSE,
                   return_models = FALSE,
                    subset = 1:100,
                   cores = ncores-2)


model_1$summaries

model_1$summaries$ENST00000023939.8_6_20



#Remove all that have output NULL
excl_1<- names(which(model_1$summaries == "NULL"))
geneids_1 <- names(which(model_1$summaries != "NULL"))



# Collect all model summaries
mod_sum_1 <- bind_rows(within(model_1$summaries, rm(excl_1))) %>%
  mutate(target = rep(geneids_1, each = 3)) %>%
  #subset(coef != "(Intercept)")  %>%
  mutate(adj.p = p.adjust( p.val, method = "fdr"),
         log2fc = estimate/log(2),
         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns"),
         odds_ratio = exp(estimate))


# Bind all model evaluations
mod_eval_1 <- bind_rows(within(model_1$evaluations, rm(excl_1)))%>%
  mutate(target = geneids_1)


# Merge the evaluations and summaries
model_cont_1 <- mod_sum_1 %>%
  inner_join(mod_eval_1, by = "target") # %>%
# filter(adj.p <= 0.05)



# filter based on p values
filt_pre_group <- model_cont_1 %>%
  filter(p.val<= 0.05  ) %>%
  subset(coef == "scaled_age") %>%
  drop_na()

ggplot(filt_pre_group, aes(x = estimate, y = odds_ratio)) +
  geom_point()+
  geom_ribbon(aes(ymin = exp(cil), ymax = exp(ciu)),fill = "grey", alpha = 0.5) 



 # Test the exercise model

# Load the full metadata 
all_full_metadata <- readRDS("data_new/processed_data/all_full_metadata.RDS")%>%
  mutate(group = case_when(age <=20 ~ "20 and below" ,
                           age > 20 & age < 30 ~ "21 to 29",
                           age >= 30 & age < 40 ~ "30 to 39", 
                           age >= 40 & age < 50 ~ "40 to 49",
                           age >= 50 & age < 60 ~ "50 to 59",
                           age >= 60 & age < 70 ~ "60 to 69",
                           age >= 70 & age < 80 ~ "70 to 79",
                           age >= 80 ~ "80 and above")) %>%
  mutate(group = factor(group, levels = c("20 and below", "21 to 29", "30 to 39",
                                          "40 to 49", "50 to 59",  "60 to 69",
                                          "70 to 79", "80 and above" )))

# Standardize the age by scaling them 0 to 1

all_full_metadata$scaled_age <- round(rescale(all_full_metadata$age), digits = 2)
# Load full splice data
all_splice_df <- readRDS("data_new/processed_data/all_splice_data.RDS")



# Evaluate the impact of RT on the full dataset


# Select only the splicing data whose metadata is available
all_full_metadata <- all_full_metadata %>%
  filter((seq_sample_id %in% colnames(all_splice_df[,-1]))) 


# REORDER THE SEQUENCE ID TO MATCH BOTH DATAFRAMMES
all_splice_reordered <- all_splice_df[, c("transcript_ID",all_full_metadata$seq_sample_id)]


# Check if everything matches except the transcript_id
match(colnames(all_splice_reordered), all_full_metadata$seq_sample_id)

# convert the 1.0 to 0.999. This is becasue beta-model accepts only values between 0 and one
all_splice_reordered[all_splice_reordered == 1 ] <- 0.999


# Since only few of the ds by age introns were contained in the full data, might be nice to look at the 
# Impact of age and exercise in one go


args_full <-list(formula = y ~  scaled_age*time + sex*time + (1|study) + (scaled_age+0|study) +(1|participant), 
                 family = glmmTMB::beta_family())



full_RT_model <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                         arguments = args_full,
                         data = all_splice_reordered,
                         metadata = all_full_metadata,
                         samplename = "seq_sample_id",
                         summary_fun = sum_fun2,
                         eval_fun = eval_mod,
                         exported = list(),
                         save_models = FALSE,
                         return_models = FALSE,
                          subset = 1:100,
                         cores = ncores-2)


full_RT_model$summaries


missing_full <- names(which(full_RT_model$summaries == "NULL"))
avail_full <- names(which(full_RT_model$summaries != "NULL"))



mod_sum <- bind_rows(within(full_RT_model$summaries, rm(missing_full))) %>%
  mutate(target = rep(avail_full, each = 6)) %>%
  subset(coef != "(Intercept)")  %>%
  mutate(adj.p = p.adjust( p.val, method = "fdr"),
         log2fc = estimate/log(2),
         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns"),
         odds_ratio = exp(estimate))


# 
# mod_eval <- bind_rows(within(full_RT_model$evaluations, rm(missing_full)))%>%
#   mutate(target = avail_full)
# 
# 
# 
# model_cont <- mod_sum %>%
#   inner_join(mod_eval, by = "target") # %>%
# filter(Pr...z..<= 0.05)

filt_pre_group <- mod_sum %>%
  filter(p.val<= 0.05  ) %>%
  subset(coef != "(Intercept)") %>%
  drop_na()

ggplot(filt_pre_group, aes(x = estimate, y = odds_ratio)) +
  geom_point()+
  geom_ribbon(aes(ymin = exp(cil), ymax = exp(ciu)),fill = "grey", alpha = 0.5) +
  facet_wrap(~coef)






args_group <-list(formula = y ~  group*time + (1|study) + (1|participant), 
                  family = glmmTMB::beta_family())



group_RT_model <- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                          arguments = args_group,
                          data = all_splice_reordered,
                          metadata = all_full_metadata,
                          samplename = "seq_sample_id",
                          summary_fun = sum_fun,
                          eval_fun = eval_mod,
                          exported = list(),
                          save_models = FALSE,
                          return_models = FALSE,
                          # subset = 1:10,
                          cores = ncores-2)


group_RT_model$summaries


missing_group <- names(which(group_RT_model$summaries == "NULL"))
avail_group <- names(which(group_RT_model$summaries != "NULL"))



mod_sum_group <- bind_rows(within(group_RT_model$summaries, rm(missing_group))) %>%
  mutate(target = rep(avail_group, each = 16)) %>%
  subset(coef != "(Intercept)")  %>%
  mutate(.by = coef,
         adj.p = p.adjust(Pr...z.., method = "fdr"),
         log2fc = Estimate/log(2),
         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) 



mod_eval_group <- bind_rows(within(group_RT_model$evaluations, rm(missing_group)))%>%
  mutate(target = avail_group)



model_cont_group <- mod_sum_group %>%
  inner_join(mod_eval_group, by = "target") 