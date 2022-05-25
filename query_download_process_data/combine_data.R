library(tximport)
library(biomaRt)
library(tidyverse)

mart <- useEnsembl(biomart = "ensembl",
                   dataset = "celegans_gene_ensembl",
                   version = "100")

Tx.worm <- getBM(attributes=c('ensembl_transcript_id',
                              'ensembl_gene_id'),
                 mart = mart)

Tx.worm <- as_tibble(Tx.worm)
#we need to rename the two columns we just retrieved from biomart
Tx.worm <- dplyr::rename(Tx.worm, target_id = ensembl_transcript_id, 
                         gene_name = ensembl_gene_id)
bioprojects <- c()
raw_data <- c()
for (dir in list.files("../query_download_process_data/ready_for_analysis/")) {
    path <- file.path("../query_download_process_data/ready_for_analysis", dir, 
                      list.files(paste(
                        "../query_download_process_data/ready_for_analysis/", 
                        dir, sep='')), 'abundance.tsv')
    
    Txi_gene <- tximport(path, 
                         type = "kallisto", 
                         tx2gene = Tx.worm, 
                         txOut = FALSE,
                         ignoreTxVersion = FALSE)
    
    genes <- rownames(Txi_gene[["counts"]])
    samples = lapply(strsplit(path, "/"), '[', 5)
    colnames(Txi_gene[["counts"]]) <- samples
    bioprojects <- c(bioprojects, rep(dir, ncol(Txi_gene[["counts"]])))
    raw_data <- cbind(raw_data, Txi_gene[['counts']])
}

save(raw_data, file = "../common_datastore/combined_studies_raw_counts.Rdata") 
save(bioprojects, file = "../common_datastore/bioprojects.Rdata") 
