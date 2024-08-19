library(tidyverse)
library(gridExtra)
library(patchwork)
library(ggplot2)
library(ggfortify)

#source of the codes
#https://tavareshugo.github.io/data-carpentry-rnaseq/03_rnaseq_pca.html
#Load preexercise splicing data

pre_sp <- readRDS("data/model/full_splice_data.RDS")
pca_matrix <- pre_sp %>%
  column_to_rownames("transcript_ID") %>% 
  as.matrix() %>%
  t()


# Perform the PCA
sample_pca <- prcomp(pca_matrix)


pca_matrix[1:10, 1:5]
as_tibble(pca_matrix)
# Convert matrix to tibble - add colnames to a new column called "gene"
as_tibble(pca_matrix, rownames = "sample")


pc_eigenvalues <- sample_pca$sdev^2

# create a "tibble" manually with 
# a variable indicating the PC number
# and a variable with the variances
pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                         variance = pc_eigenvalues) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_eigenvalues


pc_eigenvalues %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

# The PC scores are stored in the "x" value of the prcomp object
pc_scores <- sample_pca$x


pc_scores <- pc_scores %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")

# print the result
pc_scores



pc_scores %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()


#eXPLORING THE RELATIONSHIP bewteen the splicing efficency and the pcs

pc_loadings <- sample_pca$rotation

pc_loadings <- pc_loadings %>% 
  as_tibble(rownames = "transcript_ID")

# print the result
pc_loadings



top_genes <- pc_loadings %>% 
  # select only the PCs we are interested in
  select(transcript_ID, PC1, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows
  slice(1:10) %>% 
  # pull the gene column as a vector
  pull(transcript_ID) %>% 
  # ensure only unique genes are retained
  unique()

top_genes


top_loadings <- pc_loadings %>% 
  filter(transcript_ID %in% top_genes)




# Put them together
autoplot(sample_pca)

#Load the meatadata
met_df <- readRDS("data/model/full_metadata.RDS") 

unique(met_df$study)
autoplot(sample_pca, data = met_df, colour = "time", shape ="study" )
