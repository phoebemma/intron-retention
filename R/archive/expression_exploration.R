library(dplyr)
library(tidyverse)
library(ggplot2)


low_SE <- readRDS("data_new/Pre_Exercise/annotated_low_SE_introns.RDS")

high_Se <- readRDS("data_new/Pre_Exercise/annotated_perfect_SE_introns.RDS")
length(unique(high_Se$transcript_ID))
# Load baseline metadata
pre_metadata <- readRDS("data_new/Pre_Exercise/all_prexercise_metadata.RDS")

normalized_counts <- readRDS("data_new/gene_counts/all_normalised_RSEM_counts.RDS")



# Select the pre-exercise splicing data
pre_intersect <- intersect(colnames(normalized_counts), pre_metadata$seq_sample_id)

pre_counts  <-normalized_counts %>%
  subset(select = c("transcript_id", pre_intersect))


# get the lowly SE genes in the normalised counts

low_SE_counts <- pre_counts[pre_counts$transcript_id %in% low_SE$transcript_ID,]

high_SE_counts <- pre_counts[pre_counts$transcript_id %in% high_Se$transcript_ID,]


RT_affected_genes <- readRDS("data_new/processed_data/annotation_RT_effects.RDS")

RT <- pre_counts[pre_counts$transcript_id %in% RT_affected_genes$transcript_ID,]

low_transcripts  <- low_SE_counts %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "count",
               cols = -(transcript_id))%>%
  group_by(transcript_id)
merged <- merge(low_transcripts, low_SE %>%
                  select(transcript_ID, "SE" = mean), by.x = "transcript_id" , by.y= "transcript_ID")
merged$SE <- round(merged$SE,2)

ggplot(merged, aes(x = count, y = SE ))+
  geom_point()




high_transcripts <- high_SE_counts %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "count",
               cols = -(transcript_id))
merged_high <- merge(high_transcripts,high_Se %>%
                       select(transcript_ID, "SE" = mode), by.x = "transcript_id" , by.y= "transcript_ID" )
ggplot(merged_high, aes(x = count, y = SE ))+
  geom_point()



# Visualize all the captured introns by their expression

#Load all baseline introns

pre_int <- readRDS("data_new/Pre_Exercise/all_pre_Exc_splicing_data.RDS") %>%
  separate("transcript_ID", c("transcript_ID", "intron_ID", "chr"), sep = "_") %>%
  select(-c(chr, intron_ID)) %>%
  pivot_longer(names_to = "seq_sample_id",
                              values_to = "SE",
                              cols = -(transcript_ID))

long_pre_counts <- pre_counts %>%
  pivot_longer(names_to = "seq_sample_id",
               values_to = "transcript_count",
               cols = -(transcript_id))


merged_full_pre <- long_pre_counts %>%
  inner_join(pre_int, by = c("transcript_id" = "transcript_ID", "seq_sample_id"))

ggplot(merged_full_pre, aes(x = log10(transcript_count), y = log10(SE) ))+
  geom_point()
