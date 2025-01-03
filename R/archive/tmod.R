library(tmod)
library(msigdbr)
library(dplyr)

# Obtain gene sets from the msigdbr database
msig <- msigdbr()
msig <- makeTmodFromDataFrame(df = msig, feature_col = "gene_symbol",
                              module_col = "gs_id", title_col = "gs_name",
                              extra_module_cols = c("gs_cat", "gs_subcat", "gs_url",
                                                    "gs_exact_source", "gs_description"))
# This would perform gene ontlogy analyses of the ds sliced genes

ds_age <- readRDS("data/re_models/primary_model_extracts/annotated_improved_ds_introns_by_RT.RDS") 

gene_list <- ds_age %>%
  mutate(ciu = Estimate + qnorm(0.975) * Std..Error,
         cil = Estimate - qnorm(0.975) * Std..Error,
         msd = ifelse(Estimate > 0, cil, -ciu)) %>%
  # rank in descending order
  arrange(-msd) %>%
  
  #extract the ranked list of the protein-coding genes
  pull(external_gene_name)

x <- tmodCERNOtest(gene_list)


length(unique(all_pre_metadata$participant))

x <- all_pre_metadata %>%
  filter(age_group == "Young")
length(unique(x$participant))
