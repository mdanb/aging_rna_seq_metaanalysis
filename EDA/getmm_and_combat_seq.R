library(sva)
library(biomaRt) # an alternative for annotation
library(edgeR)
library(tidyverse)
library(ensembldb) # helps deal with ensembl

load(file = "../common_datastore/bioprojects.Rdata")
load(file = "../common_datastore/combined_studies_raw_counts.Rdata")

# GeTMM normalization
getmm_normalize <- function(gene_exp_data, lengths) {
  rpk <- gene_exp_data / lengths[rownames(gene_exp_data), 1]
  group <- c(rep("A", ncol(gene_exp_data)))
  rpk.norm = DGEList(counts=rpk, group=group)
  rpk.norm <- calcNormFactors(rpk.norm)
  norm.counts.rpk_edger <- cpm(rpk.norm)
  norm.counts.rpk_edger.log2 <- log2(norm.counts.rpk_edger)
  norm.counts.rpk_edger.log2[!is.finite(norm.counts.rpk_edger.log2)] <- 0
  return(t(norm.counts.rpk_edger.log2))
}

getGeneLength <- function(id, org) {
  id.type <- "ensembl"
  inp.id <- id
  
  id.type <- paste0(id.type, ifelse(id.type=="entrez", "gene", "_gene_id"))
  
  # setting mart
  message("Connecting to BioMart ...")
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
  ds <- listDatasets(ensembl)[,"dataset"]
  ds <- grep(paste0("^", org), ds, value=TRUE)
  if(length(ds) == 0)
    stop(paste("Mart not found for:", org))
  else if(length(ds) > 1)
  {
    message("Found several marts")
    sapply(ds, function(d)
      message(paste(which(ds==d), d, sep=": ")))
    n <- readline(paste0("Choose mart (1-", length(ds),") : "))
    ds <- ds[as.integer(n)]
  }
  
  ensembl <- useDataset(ds, mart=ensembl)
  
  message( paste0( "Downloading sequence",
                   ifelse(length(id) > 1, "s", ""), " ..."))
  if(length(id) > 100) message("This may take a few minutes ...")
  
  # download sequence
  # (1) get exon coordinates
  attrs <- c(id.type, "ensembl_exon_id",
             "chromosome_name", "exon_chrom_start", "exon_chrom_end")
  coords <- getBM(filters=id.type, attributes=attrs, values=id, mart=ensembl)
  id <- unique(coords[,id.type])
  coords <- GRangesList(sapply(id,
                               function(i)
                               {
                                 i.coords <- coords[coords[,1]== i, 3:5]
                                 g <- GRanges(i.coords[,1], 
                                              IRanges(i.coords[,2],
                                                      i.coords[,3]))
                                 return(g)
                               }), compress=FALSE)
  coords <- reduce(coords)
  len <- sum(width(coords))
  
  # (2) get genes and sequences
  sel <- c(id.type, "start_position", "end_position")
  gene.pos <- getBM(attributes = sel, filters=id.type, values=id,
                    mart=ensembl)
  gene.seqs <- getSequence(id=id,
                           type=id.type, seqType="gene_exon_intron", 
                           mart=ensembl)
  res <- data.frame(len)
  return(res)
}

if (!file.exists("lengths.Rdata")) {
  worm.anno <- useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                       dataset = "celegans_gene_ensembl")
  
  Tx.worm <- getBM(attributes=c('ensembl_transcript_id',
                                'ensembl_gene_id'),
                   mart = worm.anno)
  
  Tx.worm <- as_tibble(Tx.worm)
  # we need to rename the two columns we just retrieved from biomart
  Tx.worm <- dplyr::rename(Tx.worm, target_id = ensembl_transcript_id, 
                           gene_name = ensembl_gene_id)
  
  lengths = getGeneLength(c(Tx.worm["gene_name"]), "cel")
  save(lengths, file="lengths.Rdata")
}

load(file = "lengths.Rdata")

all_samples <- colnames(combined_studies)
sra_to_bioproject <- cbind(all_samples, bioprojects)
colnames(sra_to_bioproject) <- c('Sample', 'Bioproject')
write.csv(sra_to_bioproject, row.names = FALSE,
          file = "../common_datastore/sra_to_bioproject.csv")

getmm_gene_expression_all_data <- getmm_normalize(combined_studies, lengths)
if (!file.exists("../common_datastore/getmm_gene_expression_all_data.csv")) {
  write.csv(getmm_gene_expression_all_data, 
            file="../common_datastore/getmm_gene_expression_all_data.csv")
}

bioprojects_no_outliers <- bioprojects[!(bioprojects %in% c("PRJNA575569"))]
combined_gene_expression_no_outliers <- as.matrix(dplyr::select(
                                                as.data.frame(combined_studies), 
                                                -(SRR10222650:SRR10222661)))

getmm_gene_expression_no_outliers <- getmm_normalize(
                                  combined_gene_expression_no_outliers, lengths)

if (!file.exists("../common_datastore/getmm_gene_expression_no_outliers.csv")) {
  write.csv(getmm_gene_expression_no_outliers, 
            file="../common_datastore/getmm_gene_expression_no_outliers.csv")
}

bioprojects_no_singles <- bioprojects[!(bioprojects %in% 
                                          c("PRJNA261420", "PRJNA282784"))]
combined_gene_expression_no_singles <- as.matrix(dplyr::select(
                                                as.data.frame(combined_studies), 
                                                -(c(SRR2005820, SRR1578745))))

bioprojects_no_outliers_and_singles <- bioprojects_no_singles[
                                !(bioprojects_no_singles %in% c("PRJNA575569"))]
combined_gene_expression_no_outliers_and_singles <- as.matrix(dplyr::select(
                            as.data.frame(combined_gene_expression_no_singles), 
                            -(SRR10222650:SRR10222661)))

if (!file.exists(paste("../common_datastore/", 
                  "combat_seq_gene_expression_no_outliers_and_singles.Rdata",
                  sep=""))) {
  combat_seq_no_outliers_and_singles_gene_expression <- ComBat_seq(
                               combined_gene_expression_no_outliers_and_singles, 
                               batch=as.numeric(factor(
                               bioprojects_no_outliers_and_singles)), 
                               group=NULL)
  save(combat_seq_no_outliers_and_singles_gene_expression, 
       file=paste("../common_datastore/",
                  "combat_seq_gene_expression_no_outliers_and_singles.Rdata",
                  sep=""))
}

load(file = paste("../common_datastore/",
                  "combat_seq_gene_expression_no_outliers_and_singles.Rdata",
                  sep=""))
getmm_combat_seq_no_outliers_and_singles_gene_expression <- 
    getmm_normalize(combat_seq_no_outliers_and_singles_gene_expression, lengths)
if (!file.exists(paste("../common_datastore/",
                 "getmm_combat_seq_no_outliers_and_singles_gene_expression.csv",
                 sep=""))) {
  write.csv(getmm_combat_seq_no_outliers_and_singles_gene_expression, 
            file=paste("../common_datastore/",
                "getmm_combat_seq_no_outliers_and_singles_gene_expression.csv",
                sep=""))
}

# All data age batch correction ----
age <- read.csv("../common_datastore/age.csv")
samples <- colnames(combined_gene_expression_no_outliers)
temp <- as_tibble(t(combined_gene_expression_no_outliers))
x <- temp %>%
  add_column(samples, .before=colnames(temp)[1]) %>%
  rename(samples = 'Sample')
y <- as_tibble(age)
tib <- left_join(x, y)
age <- tib['age']$age

if (!file.exists(paste("../common_datastore/",
              "combat_seq_age_corrected_gene_expression_no_outliers.Rdata",
              sep=""))) {
  combat_seq_age_adjusted_no_outliers <- ComBat_seq(
                                          combined_gene_expression_no_outliers, 
                                          batch=as.numeric(factor(age)), 
                                          group=NULL)
  save(combat_seq_age_adjusted_no_outliers, 
  file=paste("../common_datastore/",
            "combat_seq_age_corrected_gene_expression_no_outliers.Rdata", sep=
              ""))
}

combat_seq_age_corrected_getmm_gene_expression_no_outliers <- 
  getmm_normalize(combat_seq_age_adjusted_no_outliers, lengths)
write.csv(combat_seq_age_corrected_getmm_gene_expression_no_outliers, 
          file = paste("../common_datastore/", 
          "combat_seq_age_corrected_getmm_gene_expression_no_outliers.csv",
          sep=""))

# All data experiment AND age double batch correction ----
age <- read.csv("../common_datastore/age.csv")
samples <- colnames(combined_gene_expression_no_outliers_and_singles)
temp <- as_tibble(t(combined_gene_expression_no_outliers_and_singles))
x <- temp %>%
  add_column(samples, .before=colnames(temp)[1]) %>%
  rename(samples = 'Sample')
y <- as_tibble(age)
tib <- left_join(x, y)
age <- tib['age']$age

combat_seq_no_outliers_and_singles_gene_expression <- ComBat_seq(
  combined_gene_expression_no_outliers_and_singles, 
  batch=as.numeric(factor(
    bioprojects_no_outliers_and_singles)), 
  group=NULL)

combat_seq_experiment_and_age_adjusted_no_outliers <- ComBat_seq(
  combat_seq_no_outliers_and_singles_gene_expression, 
  batch=as.numeric(factor(age)), 
  group=NULL)

save(combat_seq_experiment_and_age_adjusted_no_outliers, 
     file=paste("../common_datastore/",
                "combat_seq_experiment_and_age_corrected_gene_expression_no_outliers.Rdata", sep=
                  ""))

combat_seq_experiment_and_age_corrected_getmm_gene_expression_no_outliers <- 
  getmm_normalize(combat_seq_experiment_and_age_adjusted_no_outliers, lengths)
write.csv(combat_seq_experiment_and_age_corrected_getmm_gene_expression_no_outliers, 
          file = paste("../common_datastore/", 
                       "combat_seq_experiment_and_age_corrected_getmm_gene_expression_no_outliers.csv",
                       sep=""))

# Day 1 and older experiment batch corrected ----
combined_studies_day_1_and_older <- 
  combined_gene_expression_no_outliers_and_singles[, age >= 1]
temp <- t(combined_studies_day_1_and_older)
Sample <- rownames(temp)
combined_studies_day_1_and_older <- as_tibble((t(combined_studies_day_1_and_older))) %>%
                                    add_column(Sample, .before="WBGene00000001") %>%
                                    left_join(sra_to_bioproject, copy=TRUE)
batch_var <- as.list(combined_studies_day_1_and_older[, 22115])
combined_studies_day_1_and_older <- as.matrix(combined_studies_day_1_and_older[, 2:22114])
rownames(combined_studies_day_1_and_older) <- Sample
combined_studies_day_1_and_older <- t(combined_studies_day_1_and_older)
if (!file.exists("../common_datastore/combat_seq_day_1_and_older_experiment_corrected_no_outliers.Rdata")) {
  combat_seq_day_1_and_older_experiment_corrected_no_outliers <- 
    ComBat_seq(combined_studies_day_1_and_older, batch=as.numeric(factor(
      batch_var$Bioproject)), group=NULL)
  save(combat_seq_day_1_and_older_experiment_corrected_no_outliers, 
       file=paste("../common_datastore/",
                  "combat_seq_day_1_and_older_experiment_corrected_no_outliers.Rdata", sep=""))
}

combat_seq_getmm_day_1_and_older_experiment_corrected_no_outliers <- 
  getmm_normalize(combat_seq_day_1_and_older_experiment_corrected_no_outliers, lengths)
write.csv(combat_seq_getmm_day_1_and_older_experiment_corrected_no_outliers, 
          file = paste("../common_datastore/",
                       "combat_seq_experiment_corrected_day_1_and_older_getmm_gene_expression_no_outliers.csv", 
                       sep=""))

# L4 and younger experiment batch corrected ----
combined_studies_L4_and_younger <- 
  combined_gene_expression_no_outliers_and_singles[, age < 1]
temp <- t(combined_studies_L4_and_younger)
Sample <- rownames(temp)
combined_studies_L4_and_younger <- as_tibble((t(combined_studies_L4_and_younger))) %>%
  add_column(Sample, .before="WBGene00000001") %>%
  left_join(sra_to_bioproject, copy=TRUE)
batch_var <- as.list(combined_studies_L4_and_younger[, 22115])
combined_studies_L4_and_younger <- as.matrix(combined_studies_L4_and_younger[, 2:22114])
rownames(combined_studies_L4_and_younger) <- Sample
combined_studies_L4_and_younger <- t(combined_studies_L4_and_younger)
if (!file.exists("../common_datastore/combat_seq_L4_and_younger_experiment_corrected_no_outliers.Rdata")) {
  combat_seq_L4_and_younger_experiment_corrected_no_outliers <- 
    ComBat_seq(combined_studies_L4_and_younger, batch=as.numeric(factor(
      batch_var$Bioproject)), group=NULL)
  save(combat_seq_L4_and_younger_experiment_corrected_no_outliers, 
       file=paste("../common_datastore/",
                  "combat_seq_L4_and_younger_experiment_corrected_no_outliers.Rdata", sep=""))
}

combat_seq_getmm_L4_and_younger_experiment_corrected_no_outliers <- 
  getmm_normalize(combat_seq_L4_and_younger_experiment_corrected_no_outliers, lengths)
write.csv(combat_seq_getmm_L4_and_younger_experiment_corrected_no_outliers, 
          file = paste("../common_datastore/",
                       "combat_seq_experiment_corrected_L4_and_younger_getmm_gene_expression_no_outliers.csv", 
                       sep=""))

# For GO gene filtering
write.csv(combined_gene_expression_no_outliers_and_singles, file = "../common_datastore/raw_gene_expression_no_outliers_and_singles_for_GO_filtering.csv")

# After GO filtering in Python
GO_filtered_raw_gene_expression_data <- read.csv("../common_datastore/GO_filtered_raw_gene_expression_data.csv")
filtered_genes <- GO_filtered_raw_gene_expression_data[, 1]
all_genes <- row.names(lengths)
lengths_filtered <- as_tibble(lengths) %>%
                    add_column(all_genes, .before=colnames(lengths)[1]) %>%
                    rename(all_genes = 'Gene') %>%
                    dplyr::filter(Gene %in% filtered_genes)
lengths_filtered <- as.matrix(lengths_filtered[, 2])
row.names(lengths_filtered) <- filtered_genes
data <- as.matrix(GO_filtered_raw_gene_expression_data[, 2:dim(GO_filtered_raw_gene_expression_data)[2]])
combat_seq_GO_filtered <- ComBat_seq(data, batch=as.numeric(factor(bioprojects_no_outliers_and_singles)), group=NULL)
row.names(combat_seq_GO_filtered) = filtered_genes
combat_seq_getmm_GO_filtered <- getmm_normalize(combat_seq_GO_filtered, lengths_filtered)
write.csv(combat_seq_getmm_GO_filtered, file = "../common_datastore/combat_seq_getmm_GO_filtered_gene_expression_no_singles_and_outliers.csv")

# With singles for deep learning batch correction approaches ----
getmm_combined_studies_no_outliers <- as.matrix(dplyr::select(
                                        as.data.frame(combined_studies), 
                                        -(SRR10222650:SRR10222661)))
getmm_only <- getmm_normalize(getmm_combined_studies_no_outliers, lengths)
bioprojects_no_outliers <- bioprojects[!(bioprojects %in% c("PRJNA575569"))]

write_to_files <- function(df) {
  write.table(df[, 3:dim(df)[2]], paste0('../SAUCIE/getmm_corrected/gene_expression_', 
                                         unique(df$Bioproject), ".csv"), 
              row.names = FALSE, col.names = FALSE, sep = ',')
}

Sample <- row.names(getmm_only)

getmm_corrected <- as_tibble(getmm_only) %>%
                    add_column(Sample, .before="WBGene00000001")

getmm_corrected <- as_tibble(sra_to_bioproject) %>%
                    right_join(getmm_corrected)
                                    
getmm_corrected <- getmm_corrected %>%
                    group_by(Bioproject) %>%
                    do(write_to_files(.))

# Count filtered
cpm <- cpm(combined_gene_expression_no_outliers_and_singles)
keepers <- rowSums(cpm>1)>=236
count_filtered_gene_expression_no_outliers_and_singles <- 
                     combined_gene_expression_no_outliers_and_singles[keepers, ]
combat_seq_count_filtered_gene_expression_no_outliers_and_singles <- 
  ComBat_seq(count_filtered_gene_expression_no_outliers_and_singles, 
             batch=as.numeric(factor(bioprojects_no_outliers_and_singles)), 
             group=NULL)
getmm_combat_seq_count_filtered_gene_expression_no_outliers_and_singles <- 
  getmm_normalize(combat_seq_count_filtered_gene_expression_no_outliers_and_singles, lengths)
if (!file.exists("../common_datastore/getmm_combat_seq_count_filtered_gene_expression_no_outliers_and_singles.csv")) {
  write.csv(getmm_combat_seq_count_filtered_gene_expression_no_outliers_and_singles, file="../common_datastore/getmm_combat_seq_count_filtered_gene_expression_no_outliers_and_singles.csv")
}