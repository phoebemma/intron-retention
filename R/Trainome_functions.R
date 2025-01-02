

#Function for summarising seq_model summary 
sum_fun <- function(x){
  
  cond_effects <- data.frame(cbind(data.frame(coef = rownames(coef(summary(x))$cond))),
                             coef(summary(x))$cond, 
                             
                             row.names = NULL)
  
  return(cond_effects)
  
}


#Function for evaluating seq_model output
eval_mod <- function(x) {
  
  sim <- DHARMa::simulateResiduals(x, n = 1000)
  
  disp <- DHARMa::testDispersion(sim, plot = FALSE)
  unif <- DHARMa::testUniformity(sim, plot = FALSE)
  zinfl <- DHARMa::testZeroInflation(sim, plot = FALSE)
  
  results <- data.frame(pval.disp = disp$p.value, 
                        pval.unif = unif$p.value, 
                        pval.zinfl = zinfl$p.value)
  
  return(results)
}




#check number of cores. Used in seq_model
ncores <- parallel::detectCores()


model_eval <-  function(x){ 
  bind_rows(x$evaluations) %>%
    mutate(target = names(x$evaluations))
  
}


#Function to extract model summary. Takes s input the model file and the number of coefficients
model_sum <- function(x, y){
  bind_rows(x$summaries) %>%
    mutate(target = rep(names(x$summaries), each = y))%>%
    subset(!coef == "(Intercept)") %>%
    mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
           log2fc = Estimate/log(2),
           
           fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns"))
}




#function to read Kallisto files
read_kallisto_output <- function(file){
  df <- readr::read_table(file)
  df <- separate(df, target_id, c("transcript_ID", "gene_ID", "Havana_gene_ID",
                                  "Havana_transcript_ID", "transcript_name",
                                  "gene_name", "sequence_length", 
                                  "transcript_biotype"), sep = "\\|")
  df <- df %>% select(transcript_ID, gene_ID, transcript_name, gene_name,
                      transcript_biotype, length, est_counts, tpm)
  #extract all transcripts with est_count below 1
  # df <- filter(df, est_counts > 1)
  return(df)
}





## a function to extract all the lncRNAs based on EBI's definition in the link
#https://www.ensembl.org/info/genome/genebuild/biotypes.html
read_biotype_lncRNAs <- function(file){
  df <- read_kallisto_output(file)
  df <- df %>% filter(transcript_biotype %in% c("processed_transcript", "lncRNA",
                                                "lincRNA", "3prime_overlapping_ncrna",
                                                "antisense", "non_coding", "sense_intronic",
                                                "sense_overlapping", "TEC", "known_ncrna",
                                                "bidirectional_promoter_lncrna",
                                                "macro_lncRNA"))
  return(df)
}






#function to extract those that are lncRNAs among the lncRNA biotype
read_lncRNAs <- function(file){
  df <- read_kallisto_output(file)
  df <- filter(df, transcript_biotype == "lncRNA")
  return(df)
}





#A function to load load non-coding RNAs and protein coding RNAs
read_coding_and_lncRNAs <- function(file){
  df <- read_kallisto_output(file)
  df <- df %>% filter(transcript_biotype %in% c("processed_transcript","protein_coding", "lncRNA",
                                                "lincRNA", "3prime_overlapping_ncrna",
                                                "antisense", "non_coding", "sense_intronic",
                                                "sense_overlapping", "TEC", "known_ncrna",
                                                "bidirectional_promoter_lncrna",
                                                "macro_lncRNA"))
  return(df)
}




## Using data table with large number of rows instead...
# function to remove all columns except TPM and transcript name
# the function also combines all files into a data fram

####The subfunction which extracts the desired biotype is interchangable

####extract all transcripts
extract_all_transcripts <- function(folder){
  
  files <- list.files(folder, pattern =".tsv") 
  
  transcripts <- list()
  
  for(i in 1:length(files)){
    
    transcripts[[i]] <- read_kallisto_output(paste0(folder,"/", files[i])) %>% 
      mutate(file_id = gsub("^[^-]*-", "", gsub(".tsv", "", files[i]))) %>% ## This removes the file number
      dplyr::select(transcript_ID, file_id, tpm)
  }
  
  comd.df <- data.table::rbindlist(transcripts)
  
  comb.df <- data.table::dcast(comd.df, transcript_ID ~ file_id, value.var = "tpm")
  
  return(data.frame(comb.df))
}





#extracts proteing coding and lncRNAs 

extract_coding_lncRNA <- function(folder){
  
  files <- list.files(folder, pattern =".tsv") 
  
  imports <- list()
  
  for(i in 1:length(files)){
    
    imports[[i]] <- read_coding_and_lncRNAs(paste0(folder,"/", files[i])) %>% 
      mutate(file_id = gsub("^[^-]*-", "", gsub(".tsv", "", files[i]))) %>% ## This removes the file number
      dplyr::select(transcript_ID, file_id, tpm)
    
    
    
  }
  
  comd.df <- data.table::rbindlist(imports)
  
  comb.df <- data.table::dcast(comd.df, transcript_ID ~ file_id, value.var = "tpm")
  
  return(data.frame(comb.df))
  
  
}





#Extracts only lncRNAs

extract_lncRNAs <- function(folder){
  
  files <- list.files(folder, pattern =".tsv") 
  
  imports1 <- list()
  
  for(i in 1:length(files)){
    
    imports1[[i]] <- read_biotype_lncRNAs(paste0(folder,"/", files[i])) %>% 
      mutate(file_id = gsub("^[^-]*-", "", gsub(".tsv", "", files[i]))) %>% ## This removes the file number
      dplyr::select(transcript_ID, file_id, tpm)
    
    
    
  }
  
  comd.df <- data.table::rbindlist(imports1)
  
  comb.df <- data.table::dcast(comd.df, transcript_ID ~ file_id, value.var = "tpm")
  
  return(data.frame(comb.df))
  
}






#function to read Rsem .gene.results files
read_Rsem_genes <- function(file){
  df <- readr::read_delim(file)
  df <- df %>% select(gene_id, length, effective_length, expected_count)
  return(df)
}



#function to read Rsem isoform.results file
read_Rsem_isoforms <- function(file){
  df <- readr::read_delim(file)
  df <- df %>% select(transcript_id, effective_length, length, expected_count,)
  return(df)
}




#Function to extract the rsem isoforms in all the  files in a folder
extract_rsem_isoform_counts <- function(folder){
  
  files <- list.files(folder, pattern =".isoforms.results") 
  
  isoforms <- list()
  
  for(i in 1:length(files)){
    
    isoforms[[i]] <- read_Rsem_isoforms(paste0(folder,"/", files[i])) %>% 
      mutate(file_id =gsub(".isoforms.results", "", files[i])) %>% ## This removes the file number
      dplyr::select(transcript_id, file_id, expected_count, effective_length, length )
  }
  
  comd.df <- data.table::rbindlist(isoforms)
  
  comb.df <- data.table::dcast(comd.df, transcript_id ~file_id, value.var = "expected_count")
  
  return(data.frame(comb.df))
  
}






#function to extract the rsem-generated genes in all files in a folder
extract_rsem_gene_counts <- function(folder){
  
  files <- list.files(folder, pattern =".genes.results") 
  
  genes <- list()
  
  for(i in 1:length(files)){
    
    genes[[i]] <- read_Rsem_genes(paste0(folder,"/", files[i])) %>% 
      mutate(file_id = gsub(".genes.results", "", files[i])) %>% ## This removes the file number
      dplyr::select(gene_id, file_id, length, effective_length, expected_count)
    
    
    
  }
  
  comd.df <- data.table::rbindlist(genes)
  
  
  comb.df <- data.table::dcast(comd.df, gene_id ~ file_id, value.var = "expected_count")
  
  return(data.frame(comb.df))
}






#extract the length of the genes. This extracts the column called "effective length"
extract_rsem_gene_lengths <- function(folder){
  
  files <- list.files(folder, pattern =".genes.results") 
  
  genes <- list()
  
  for(i in 1:length(files)){
    
    genes[[i]] <- read_Rsem_genes(paste0(folder,"/", files[i])) %>% 
      mutate(file_id = gsub(".genes.results", "", files[i])) %>% ## This removes the file number
      dplyr::select(gene_id, file_id, length, effective_length, expected_count)
    
    
    
  }
  
  comd.df <- data.table::rbindlist(genes)
  
  
  comb.df <- data.table::dcast(comd.df, gene_id ~ file_id, value.var = "effective_length")
  
  return(data.frame(comb.df))
}





extract_rsem_isoform_lengths <- function(folder){
  
  files <- list.files(folder, pattern =".isoforms.results") 
  
  isoforms <- list()
  
  for(i in 1:length(files)){
    
    isoforms[[i]] <- read_Rsem_isoforms(paste0(folder,"/", files[i])) %>% 
      mutate(file_id =gsub(".isoforms.results", "", files[i])) %>% ## This removes the file number
      dplyr::select(transcript_id, file_id, expected_count, effective_length, length )
  }
  
  comd.df <- data.table::rbindlist(isoforms)
  
  comb.df <- data.table::dcast(comd.df, transcript_id ~file_id, value.var = "effective_length")
  
  return(data.frame(comb.df))
  
}






#function to read Splice_q results file
read_Splice_Q <- function(file){
  df <- readr::read_tsv(file)
  df <- df %>% dplyr::select(transcript_ID, intron_ID, score, chr)
  df$transcript_ID <- paste0(df$transcript_ID,"_",
                             df$intron_ID,"_",
                             df$chr) 
  df %>% dplyr::select(transcript_ID, score)
  return(df)
}



#Function to extract all the spliceQ data in a folder
extract_splice_q <- function(folder){
  
  files <- list.files(folder, pattern =".tsv") 
  
  splice_list <- list()
  
  for(i in 1:length(files)){
    
    splice_list[[i]] <- read_Splice_Q(paste0(folder, files[i])) %>% 
      mutate(file_id = gsub(".tsv", "", files[i])) %>% ## This removes the file number
      dplyr::select(transcript_ID, file_id, score)
  }
  
  comd.df <- data.table::rbindlist(splice_list)
  
  comb.df <- data.table::dcast(comd.df, transcript_ID  ~ file_id, value.var = "score")
  
  return(data.frame(comb.df))
}