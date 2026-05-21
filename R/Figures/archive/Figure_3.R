library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(cowplot)
library(gridExtra)
library(grid)
library(gt)
library(biomaRt)
library(org.Hs.eg.db)
library(clusterProfiler)

# This explores the RT model outputs


# Figure 4
# Visualization of the impact of RT on Se
RT_effects <- readRDS("data_new/processed_data/annotation_RT_effects.RDS")

RT_effects %>%
  ggplot(aes(transcript_biotype))+
  geom_bar()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  stat_count(geom = "Text", aes(label = ..count..), vjust = -0.5) +
  ggtitle("Annotation of genes containing differentially spliced introns in full data") +
  ylab("Number of genes")




# Extract those ds by exercise
filt_exc <- RT_effects %>%
  filter(coef == "timePostExc")
hist(filt_exc$Estimate)


filt_int <- RT_effects %>%
  filter(coef == "scaled_age:timePostExc")
hist(filt_int$Estimate)

filt_int_neg <- filt_int %>%
  filter(Estimate < 0)

resistance <-  filt_exc %>%
  ggplot(aes(x=Estimate))+
  geom_histogram( colour= "lightblue") +
  ylab("Number of ds introns per Estimate")+
  xlab("Change in SE after RT")
ggtitle("Model Estimate of differentially spliced introns by RT")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))+
  theme_cowplot()


# Load the full model to extract the effects of aging alone
model <- readRDS("data_new/models/full_data_RT_scaled_age_int_model.RDS") %>%
  # filter to only those with  p values at or below 0.05
  filter(Pr...z..<= 0.05)

p <- model %>%
  filter(coef ==  "scaled_age")


o <- p%>%
  #filter(fcthreshold == "s") %>%
  separate("target", c("transcript_ID", "intron_ID", "chr"), sep = "_") %>%
  inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version")) %>%
  dplyr::select(Estimate, coef,  transcript_ID,intron_ID, transcript_biotype, ensembl_gene_id, external_gene_name, transcript_length)%>%
  # merge the intron_id back to the transcript ID so we track specific introns 
  mutate(intron_ID = paste0(transcript_ID, "_", intron_ID))

saveRDS(o, "data_new/models/annotated_aging_RT_model.RDS")

aging <-  o %>%
  ggplot(aes(x=Estimate))+
  geom_histogram( colour= "lightblue") +
  ylab("Number of ds introns per Estimate")+
  xlab("Change in SE due to aging")
ggtitle("Model Estimate of differentially spliced introns by aging in full dataset")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))+
  theme_cowplot()




# Annotation
ego_df_postEXC <- enrichGO(gene = filt_exc$ensembl_gene_id,
                           keyType = "ENSEMBL",
                           OrgDb = org.Hs.eg.db, 
                           ont = "BP", 
                           pAdjustMethod = "BH", 
                           qvalueCutoff = 0.05, 
                           readable = T)

## Output results from GO analysis to a table
cluster_summary_postExc <- data.frame(ego_df_postEXC)

a <- dotplot(ego_df_postEXC,
             
             font.size = 8, title = "Enriched biological processes due to RT") +
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))





ego_df_postEXC_2 <- enrichGO(gene = filt_int$ensembl_gene_id,
                             keyType = "ENSEMBL",
                             OrgDb = org.Hs.eg.db, 
                             ont = "BP", 
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.05, 
                             readable = T)

## Output results from GO analysis to a table
cluster_summary_postExc_2 <- data.frame(ego_df_postEXC_2)

b <- dotplot(ego_df_postEXC_2,
             
             font.size = 8, title = "Enriched biological processes based on interaction between age and RT") +
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))





ggarrange(a, b,
          labels = c("A", "B"))

ggsave("Images_tables/Figure4_annotation_RT.png", bg = "white" ,  scale = 2)








y <- model %>%
  filter(coef == "scaled_age:timePostExc")
x <- model %>%
  filter(coef == "timePostExc")
intersect(x$target, y$target)

c <- intersect(a$target, x$target)



aging_in_RT <- p[p$target %in% x$target,]

RT_in_aging <- x[x$target %in% p$target,]


q <-   ggplot()+
  geom_histogram(data = aging_in_RT, aes(x = Estimate),  fill= "green") +
  geom_histogram(data = RT_in_aging, aes(x = Estimate),  fill= "darkblue")+
  ylab("Number of ds introns per Estimate")+
  xlab("Change in SE by aging and by RT")+
  #  ggtitle("Introns differentially spliced by both aging and RT")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))+
  theme_cowplot()



r <-  ggplot()+
  geom_point(data = aging_in_RT, aes(x = target,  y = Estimate ), colour = "green") +
  geom_point(data = RT_in_aging, aes(x = target, y = Estimate), colour = "darkblue")+
  theme(axis.text.y = element_text(size = 8))+
  xlab("DS introns by aging and RT")+
  ylab("Change in SE by aging and by RT")+
  #   ggtitle("Introns differentially spliced by both aging and RT")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))

#   
# aging_and_RT <- RT_in_aging %>%
#   inner_join( aging_in_RT, by = "target") %>%
#   dplyr::select(target, Estimate.x, Estimate.y)
# aging_and_RT %>%
#   ggplot(aes(x = Estimate.x, y = Estimate.y)) +
#   geom_point()
# 
# cor.test(aging_and_RT$Estimate.x, aging_and_RT$Estimate.y)

# EXpplore those affected by aging but captured in interaction between aging and RT

int_in_aging <- y[y$target %in% a$target,]

aging_in_int <- a[a$target %in% y$target,]

s <-  ggplot()+
  geom_point(data = int_in_aging, aes(x = target,  y = Estimate ), colour = "red") +
  geom_point(data = aging_in_int, aes(x = target, y = Estimate), colour = "purple")+
  theme(axis.text.y = element_text(size = 8))+
  xlab("DS introns by aging and interaction with RT")+
  ylab("Change in SE by aging and interaction with RT")+
  #  ggtitle("Introns differentially spliced by aging and its interaction with RT")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))




t <-  ggplot()+
  geom_histogram(data = int_in_aging, aes(x = Estimate),  fill= "red") +
  geom_histogram(data = aging_in_int, aes(x = Estimate),  fill= "purple")+
  ylab("Number of ds introns per Estimate")+
  xlab("Change in SE by aging and interaction with RT")+
  #  ggtitle("Introns differentially spliced by aging and its interaction with RT")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))+
  theme_cowplot()



ggarrange(resistance+ rremove("xlab"), aging + rremove("xlab"), q + rremove("xlab"),r+ rremove("xlab"),t+ rremove("xlab"),s+ rremove("xlab"),
          # labels = c("A", "B", "C", "D", "E", "F"),
          align = "v", ncol = 2, nrow = 3)



ggsave("Images_tables/Figure5_ds_introns_aging_RT.png", bg = "white" ,  scale = 2.2, dpi = 400)






# EXplore tyhe intersect annotation

anno <- RT_in_aging %>%
  separate("target", c("transcript_ID", "intron_ID", "chr"), sep = "_") %>%
  inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version")) %>%
  dplyr::select(Estimate, coef,  transcript_ID,intron_ID, transcript_biotype, ensembl_gene_id, external_gene_name, transcript_length)%>%
  # merge the intron_id back to the transcript ID so we track specific introns 
  mutate(intron_ID = paste0(transcript_ID, "_", intron_ID))






# load the batch-corrected gene counts data
gene_counts <- readRDS("data_new/gene_counts/batch_corrected_genecounts.RDS") %>%
  # remove the version number to match gene_id
  mutate(gene_id = gsub("\\..*", "",  gene_id))





# Effect on age group
RT_age_group <- readRDS("data_new/models/agegroup_RT_model.RDS") %>%
  # filter to only those with adjusted p values at or below 0.05
  filter(Pr...z..<= 0.05) %>%
  separate("target", c("transcript_ID", "intron_ID", "chr"), sep = "_") %>%
  inner_join(gene_annotation, by= c("transcript_ID" = "ensembl_transcript_id_version")) %>%
  dplyr::select(coef, Estimate, transcript_ID, intron_ID,  transcript_biotype, ensembl_gene_id, external_gene_name, transcript_length) %>%
  # filter only the interaction effects
  filter(coef == "group21 to 29:timePostExc"| coef == "group30 to 39:timePostExc" |
           coef == "group40 to 49:timePostExc" | coef == "group50 to 59:timePostExc" |
           coef == "group60 to 69:timePostExc" | coef == "group70 to 79:timePostExc" |
           coef == "group80 and above:timePostExc") %>%
  # # merge the intron_id back to the transcript ID so we track specific introns 
  mutate(intron_ID = paste0(transcript_ID, "_", intron_ID),
         coef = recode(coef, "group21 to 29:timePostExc" = "21 to 29", 
                       "group30 to 39:timePostExc" = "30 to 39",
                       "group40 to 49:timePostExc" = "40 to 49",
                       "group50 to 59:timePostExc" = "50 to 59",
                       "group60 to 69:timePostExc" =  "60 to 69",
                       "group70 to 79:timePostExc" = "70 to 79",
                       "group80 and above:timePostExc" = "80 and above" )) 


# Visualize the coefficients
RT_effect <- ggplot( RT_age_group , aes(coef, fill = coef )) +
  geom_bar()+
  ggtitle("Effect of RT by age group")+
  xlab(" Age group of participants")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 10))



# Explore age group effect
# Twenties
twenties_RT <-RT_age_group %>%
  filter(coef == "21 to 29")
tw <- twenties_RT %>%
  ggplot()+
  geom_histogram(aes(x = Estimate)) +
  ylab("Number of ds introns per Estimate")+
  # xlab("Change in SE by aging and by RT")+
  ggtitle("Distribution of ds introns after RT in participants aged 21 to 29")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))+
  theme_cowplot()

 
 
 # gene ontology (biological process)  analysis of the genes containing introns with improved SE with aging
 ego_twenties_Rt <- enrichGO(gene = twenties_RT$ensembl_gene_id,
                          keyType = "ENSEMBL",
                          OrgDb = org.Hs.eg.db, 
                          universe = gene_counts$gene_id,
                          ont = "BP", 
                          pAdjustMethod = "BH", 
                          qvalueCutoff = 0.05, 
                          readable = T)
 
 
 
 twenties_RT_plot<-  dotplot(ego_twenties_Rt,
                          
                          font.size = 8, title = "Enriched biological processes following RT in participants aged 21 to 29") +
   theme(axis.text = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12))
 
 
 
 

# Thirties 
 thirties_RT <-RT_age_group %>%
   filter(coef == "30 to 39")
th <- thirties_RT %>%
   ggplot()+
   geom_histogram(aes(x = Estimate)) +
   ylab("Number of ds introns per Estimate")+
   # xlab("Change in SE by aging and by RT")+
   ggtitle("Distribution of ds introns after RT in participants aged 30 to 39")+
   theme(plot.title = element_text(hjust = 0.5),
         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))+
   theme_cowplot()
 
 
 
 # gene ontology (biological process)  analysis of the genes containing introns with improved SE with aging
 ego_thirties_Rt <- enrichGO(gene = thirties_RT$ensembl_gene_id,
                             keyType = "ENSEMBL",
                             OrgDb = org.Hs.eg.db, 
                             universe = gene_counts$gene_id,
                             ont = "BP", 
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.05, 
                             readable = T)
 
 
 
 thirties_RT_plot<-  dotplot(ego_thirties_Rt,
                             
                             font.size = 8, title = "Enriched biological processes following RT in participants aged 30 to 39") +
   theme(axis.text = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12))
 
 
 
 
 # Fourties 
 fourties_RT <-RT_age_group %>%
   filter(coef == "40 to 49")
 fo <- fourties_RT %>%
   ggplot()+
   geom_histogram(aes(x = Estimate)) +
   ylab("Number of ds introns per Estimate")+
   # xlab("Change in SE by aging and by RT")+
   ggtitle("Distribution of ds introns after RT in participants aged 40 to 49")+
   theme(plot.title = element_text(hjust = 0.5),
         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))+
   theme_cowplot()
 
 
 
 # gene ontology (biological process)  analysis of the genes containing introns with improved SE with aging
 ego_fourties_Rt <- enrichGO(gene = fourties$ensembl_gene_id,
                             keyType = "ENSEMBL",
                             OrgDb = org.Hs.eg.db, 
                             universe = gene_counts$gene_id,
                             ont = "BP", 
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.05, 
                             readable = T)
 
 
 
 fourties_RT_plot<-  dotplot(ego_fourties_Rt,
                             
                             font.size = 8, title = "Enriched biological processes following RT in participants aged 40 to 49") +
   theme(axis.text = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12))
 
 
 
 # fifties 
 fifties_RT <-RT_age_group %>%
   filter(coef == "50 to 59")
fi <-  fifties_RT %>%
   ggplot()+
   geom_histogram(aes(x = Estimate)) +
   ylab("Number of ds introns per Estimate")+
   # xlab("Change in SE by aging and by RT")+
   ggtitle("Distribution of ds introns after RT in participants aged 50 to 59")+
   theme(plot.title = element_text(hjust = 0.5),
         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))+
   theme_cowplot()
 
 
 
 # gene ontology (biological process)  analysis of the genes containing introns with improved SE with aging
 ego_fifties_Rt <- enrichGO(gene = fifties_RT$ensembl_gene_id,
                             keyType = "ENSEMBL",
                             OrgDb = org.Hs.eg.db, 
                             universe = gene_counts$gene_id,
                             ont = "BP", 
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.05, 
                             readable = T)
 
 
 
 fifties_RT_plot<-  dotplot(ego_fifties_Rt ,
                             
                             font.size = 8, title = "Enriched biological processes following RT in participants aged 50 to 59") +
   theme(axis.text = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12))
 
 
 
 # sixties 
 sixties_RT <-RT_age_group %>%
   filter(coef == "60 to 69")
si <- sixties_RT %>%
   ggplot()+
   geom_histogram(aes(x = Estimate)) +
   ylab("Number of ds introns per Estimate")+
   # xlab("Change in SE by aging and by RT")+
   ggtitle("Distribution of ds introns after RT in participants aged 60 to 69")+
   theme(plot.title = element_text(hjust = 0.5),
         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))+
   theme_cowplot()
 
 
 
 # gene ontology (biological process)  analysis of the genes containing introns with improved SE with aging
 ego_sixties_Rt <- enrichGO(gene = sixties_RT$ensembl_gene_id,
                            keyType = "ENSEMBL",
                            OrgDb = org.Hs.eg.db, 
                            universe = gene_counts$gene_id,
                            ont = "BP", 
                            pAdjustMethod = "BH", 
                            qvalueCutoff = 0.05, 
                            readable = T)
 
 
 
sixties_RT_plot<-  dotplot(ego_sixties_Rt ,
                            
                            font.size = 8, title = "Enriched biological processes following RT in participants aged 60 to 69") +
   theme(axis.text = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12))
 
 

# seventies 
seventies_RT <-RT_age_group %>%
  filter(coef == "70 to 79")
se <- seventies_RT %>%
  ggplot()+
  geom_histogram(aes(x = Estimate)) +
  ylab("Number of ds introns per Estimate")+
  # xlab("Change in SE by aging and by RT")+
  ggtitle("Distribution of ds introns after RT in participants aged 70 to 79")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))+
  theme_cowplot()



# gene ontology (biological process)  analysis of the genes containing introns with improved SE with aging
ego_seventies_Rt <- enrichGO(gene = seventies_RT$ensembl_gene_id,
                           keyType = "ENSEMBL",
                           OrgDb = org.Hs.eg.db, 
                           universe = gene_counts$gene_id,
                           ont = "BP", 
                           pAdjustMethod = "BH", 
                           qvalueCutoff = 0.05, 
                           readable = T)



seventies_RT_plot<-  dotplot(ego_seventies_Rt ,
                           
                           font.size = 8, title = "Enriched biological processes following RT in participants aged 70 to 79") +
  theme(axis.text = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12))


 
# eighties 
eighties_RT <-RT_age_group %>%
  filter(coef == "80 and above")
ei <- eighties_RT %>%
  ggplot()+
  geom_histogram(aes(x = Estimate)) +
  ylab("Number of ds introns per Estimate")+
  # xlab("Change in SE by aging and by RT")+
  ggtitle("Distribution of ds introns after RT in participants aged 80 and above")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10))+
  theme_cowplot()



# gene ontology (biological process)  analysis of the genes containing introns with improved SE with aging
ego_eighties_Rt <- enrichGO(gene = eighties_RT$ensembl_gene_id,
                             keyType = "ENSEMBL",
                             OrgDb = org.Hs.eg.db, 
                             universe = gene_counts$gene_id,
                             ont = "BP", 
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.05, 
                             readable = T)



eighties_RT_plot<-  dotplot(ego_eighties_Rt ,
                             
                             font.size = 8, title = "Enriched biological processes following RT in participants aged 80 and above") +
  theme(axis.text = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12))

 
ggarrange(RT_effect,  tw + rremove("xy.title"), th + rremove("xy.title"), fo + rremove("xy.title"), fi + rremove("xy.title"),
          si + rremove("xy.title"), se + rremove("y.title"),
          ei + rremove("y.title"),
          ncol = 2,
          nrow = 4,
          labels = c("A", "B", "C",
                     "D", "E", "F",
                     "G", "H"))
ggsave("Images_tables/Figure3_grouped_RT_effect.png", bg = "white" ,  scale = 2, dpi = 400)




ggarrange(twenties_RT_plot, thirties_RT_plot,
          fourties_RT_plot, fifties_RT_plot, 
          sixties_RT_plot, seventies_RT_plot,
          eighties_RT_plot, ncol = 2, nrow = 4,  align = "v", legend = "none")

ggsave("Images_tables/Figure3_grouped_RT_effectGO.png", bg = "white" ,  scale = 2, dpi = 400)
 