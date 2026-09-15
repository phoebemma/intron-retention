

Relief_full_meta <- readRDS("data/Relief_metadata.RDS")%>%
  dplyr::select(study, participant, sex, time, seq_sample_id, age) %>%
  mutate(sex = factor(sex, levels = c("female", "male")),
         time = factor(time, levels = c("PreExc", "PostExc")),
         age_group = case_when(age < 40 ~ "Young",
                               age > 40 ~ "Old"),
         age_group = factor(age_group, levels = c("Young", "Old")))


Relief_full_splice <- readRDS("data/Relief_splicing_data.RDS")  %>%
  drop_na()


young <- Relief_full_meta %>%
  filter(age_group == "Young")
young_splice <- intersect(colnames(Relief_full_splice), young$seq_sample_id)


young_splice_df <-Relief_full_splice %>%
  subset(select = c("transcript_ID", young_splice)) 

# which ones were zeros all through in the young
young_zero_introns <- young_splice_df %>%
  filter(if_all(-transcript_ID, ~ . == 0)) %>%
  dplyr::select(transcript_ID) %>%
  separate(transcript_ID, into = c("transcript_ID", "intron_ID", "chr"),
           sep = "_") %>%
  inner_join(gene_annotation,
             by = c("transcript_ID" = "ensembl_transcript_id_version")) %>%
  mutate(
    gene_label  = ifelse(
      is.na(external_gene_name) | external_gene_name == "",
      ensembl_gene_id, external_gene_name
    ),
    gene_intron = paste(gene_label, intron_ID, sep = " : ")
  )



young_perfect_introns <- young_splice_df %>%
  filter(if_all(-transcript_ID, ~ . == 1)) %>%
  dplyr::select(transcript_ID) %>%
  separate(transcript_ID, into = c("transcript_ID", "intron_ID", "chr"),
           sep = "_") %>%
  inner_join(gene_annotation,
             by = c("transcript_ID" = "ensembl_transcript_id_version")) %>%
  mutate(
    gene_label  = ifelse(
      is.na(external_gene_name) | external_gene_name == "",
      ensembl_gene_id, external_gene_name
    ),
    gene_intron = paste(gene_label, intron_ID, sep = " : ")
  )






old <- Relief_full_meta %>%
  filter(age_group == "Old")
old_splice <- intersect(colnames(Relief_full_splice), old$seq_sample_id)


old_splice_df <-Relief_full_splice %>%
  subset(select = c("transcript_ID", old_splice)) 

# which ones were zeros all through in the young
old_zero_introns <- old_splice_df %>%
  filter(if_all(-transcript_ID, ~ . == 0)) %>%
  dplyr::select(transcript_ID) %>%
  separate(transcript_ID, into = c("transcript_ID", "intron_ID", "chr"),
           sep = "_") %>%
  inner_join(gene_annotation,
             by = c("transcript_ID" = "ensembl_transcript_id_version")) %>%
  mutate(
    gene_label  = ifelse(
      is.na(external_gene_name) | external_gene_name == "",
      ensembl_gene_id, external_gene_name
    ),
    gene_intron = paste(gene_label, intron_ID, sep = " : ")
  )



old_perfect_introns <- old_splice_df %>%
  filter(if_all(-transcript_ID, ~ . == 1)) %>%
  dplyr::select(transcript_ID) %>%
  separate(transcript_ID, into = c("transcript_ID", "intron_ID", "chr"),
           sep = "_") %>%
  inner_join(gene_annotation,
             by = c("transcript_ID" = "ensembl_transcript_id_version")) %>%
  mutate(
    gene_label  = ifelse(
      is.na(external_gene_name) | external_gene_name == "",
      ensembl_gene_id, external_gene_name
    ),
    gene_intron = paste(gene_label, intron_ID, sep = " : ")
  )





ego_ones <- enrichGO(gene =  unique(old_perfect_introns$external_gene_name),
                     keyType = "SYMBOL",
                     universe = gene_exp_df$gene_name,
                     OrgDb = org.Hs.eg.db, 
                     ont = "BP", 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05, 
                     readable = T)


## Output results from GO analysis to a table
cluster_ones <- data.frame(ego_ones)

go_ones <- dotplot(ego_ones,
                   showCategory = 6,
                   font.size = 8, title = "Enriched biological processes in genes containing introns perfectly spliced acrosss all samples") +
  theme(axis.text = element_text(size = 16), axis.text.y = element_text(size = 16), axis.title.x = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.text = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"))

