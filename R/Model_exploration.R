#This is where the models built are explored

library(dplyr)
library(ggplot2)
library(biomaRt)
library(tmod)
#library(plotrix)

#Load the  model

#These are models from the primary data

Full_group_model <- readRDS("data/re_models/primary_full_group_unfiltered_model.RDS")%>%
  subset(coef != "(Intercept)") %>% 
  drop_na()



Full_count_model <- readRDS("data/re_models/primary_full_count_model.RDS")%>%
  subset(coef != "(Intercept)") %>% 
  drop_na()




# Filter the results for those different in old versus young
age_coef <- Full_group_model %>%
  subset(coef == "age_groupOld")%>%
  mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
         log2fc = Estimate/log(2),
         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) %>%
  filter(fcthreshold == "s" & adj.p <= 0.05 )
saveRDS(age_coef, "data/re_models/primary_model_extracts/age_group_grouped_model.RDS")



# Filter for time, ie impact of RT

time_coef <- Full_group_model %>%
  subset(coef == "timePostExc")%>%
  mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
         log2fc = Estimate/log(2),
         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) %>%
  filter(fcthreshold == "s" & Pr...z..<= 0.05 )





sex_coef <- Full_group_model %>%
  subset(coef == "sexmale")%>%
  mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
         log2fc = Estimate/log(2),
         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) %>%
  filter(fcthreshold == "s" & adj.p <= 0.05 )




# Interaction age and RT

int_coef <- Full_group_model %>%
  subset(coef == "age_groupOld:timePostExc")%>%
  mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
         log2fc = Estimate/log(2),
         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) %>%
  filter(fcthreshold == "s" & adj.p <= 0.05 )

#The target contains the transcript_id, intron_id and chromosome in thats order.

# see "read_Splice_Q" in Trainome_functions.R
colnames(pre_model)

unique(pre_model$coef)

pre_df <- pre_model %>%
  separate("target", c("transcript_ID", "intron_ID", "chr"), sep = "_") %>%
  subset(coef == "timePostExc")

#unique(pre_df$coef)
ggplot(pre_df, aes(chr)) +
  geom_bar()+
  labs(x = "Chromosomes")+
  ggtitle("chromosomes containing DR introns in baseline model")+
theme(plot.title = element_text(hjust = 0.5))

length(unique(pre_df$transcript_ID))


#Use the ensemble database to get the annotation of the transcripts
listMarts()
ensembl=useMart("ENSEMBL_MART_ENSEMBL")
datasets <- listDatasets(ensembl)
ensembl=useDataset("hsapiens_gene_ensembl",mart=ensembl)
attributes <- listAttributes(ensembl)

# #extract GENE biotypes and  names
annotation<- getBM(attributes = c("transcript_biotype", "external_transcript_name",
                                  "ensembl_gene_id", "ensembl_transcript_id_version"), values =pre_df$transcript_ID, mart = ensembl )
annotation <- inner_join(pre_df, annotation, by= c("transcript_ID" = "ensembl_transcript_id_version"))

#Explore the biotypes of the transcripts
x <- annotation %>%
  group_by(transcript_biotype)%>%
  count()%>% #count the occurence
  ungroup()%>%
  mutate(perc = n/sum(n))%>% #extract the percentage
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))
  

ggplot(x, aes(x = "", y = perc, fill = transcript_biotype)) +
  geom_col()+
  geom_text(aes(label = labels, x = 1.6),
            position = position_stack(vjust = 0.8)) +
  coord_polar(theta = "y")+
  theme_void()#+
  #ggtitle("Distribution of biotypes of transcripts containing DR introns by age")
  
# pie3D(x$perc, labels = x$transcript_biotype,# pie3D(x$percTRUE, labels = x$transcript_biotype,
#       #explode = 0.2,
#       labelcex = 0.65,
#       radius = 1, theta = pi/1)


