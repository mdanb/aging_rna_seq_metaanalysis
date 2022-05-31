# Load packages -----
library(tidyverse)
library(limma)
library(edgeR)
library(msigdbr)
library(gage)
library(fgsea)
library(data.table)
library(biomaRt)
library(tximport)

load('../common_datastore/combat_seq_gene_expression_no_outliers_and_singles.Rdata')

dir.create('results')
dir.create('../common_datastore/single_experiments/')
# Make a DGElist from your counts ----
genes <- rownames(combat_seq_no_outliers_and_singles_gene_expression)
rownames(combat_seq_no_outliers_and_singles_gene_expression) <- genes

myDGEList <- DGEList(combat_seq_no_outliers_and_singles_gene_expression)
cpm <- cpm(myDGEList)
keepers <- rowSums(cpm>1) >= 236
myDGEList.filtered <- myDGEList[keepers,]

# Normalize your data ----
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")

labels <- as_tibble(read.csv(file = '../common_datastore/labels.csv'))
age <- as_tibble(read.csv(file = '../common_datastore/age.csv'))
labels_with_age <- inner_join(labels, age)
labels_with_age_combined_data <- labels_with_age[labels_with_age$Sample %in% 
                                                   colnames(combat_seq_no_outliers_and_singles_gene_expression), ]
temp_for_reordering <- as_tibble(colnames(combat_seq_no_outliers_and_singles_gene_expression))
colnames(temp_for_reordering) <- 'Sample'
labels_with_age_combined_data <- left_join(temp_for_reordering, labels_with_age_combined_data)

group <- factor(labels_with_age_combined_data$longevity)
design <- model.matrix(~0 + group)
colnames(design) <- c('long', 'normal', 'short')

variance_stabilized_adjusted <- voom(myDGEList.filtered.norm, design, plot = TRUE)
fit <- lmFit(variance_stabilized_adjusted, design)
contrast.matrix <- makeContrasts(long_vs_short = long - short,
                                 long_vs_normal = long - normal,
                                 short_vs_normal = short - normal,
                                 levels=design)

fits <- contrasts.fit(fit, contrast.matrix)
# get bayesian stats for your linear model fit
ebFit <- eBayes(fits)

results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=1.5)
summary(results)

diffGenes_long_v_short <- variance_stabilized_adjusted$E[results[,1] != 0,]
diffGenes_long_v_normal <- variance_stabilized_adjusted$E[results[,2] != 0,]
diffGenes_short_v_normal <- variance_stabilized_adjusted$E[results[,3] != 0,]

# convert your DEGs to a dataframe using as_tibble
diffGenes.df.long.v.short <- as_tibble(diffGenes_long_v_short, rownames = "geneID")["geneID"]
diffGenes.df.long.v.normal <- as_tibble(diffGenes_long_v_normal, rownames = "geneID")["geneID"]
diffGenes.df.short.v.normal <- as_tibble(diffGenes_short_v_normal, rownames = "geneID")["geneID"]

mart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)
diffGenes.df.long.v.short.gene.names <- getBM(mart = mart, filters = "wbps_gene_id", attributes = c("entrezgene_name"), values = diffGenes.df.long.v.short)
diffGenes.df.long.v.normal.gene.names <- getBM(mart = mart, filters = "wbps_gene_id", attributes = c("entrezgene_name"), values = diffGenes.df.long.v.normal)
diffGenes.df.short.v.normal.gene.names <- getBM(mart = mart, filters = "wbps_gene_id", attributes = c("entrezgene_name"), values = diffGenes.df.short.v.normal)

write.csv(diffGenes.df.long.v.short.gene.names, 'results/diff_genes_long_vs_short.csv', row.names = FALSE)
write.csv(diffGenes.df.long.v.normal.gene.names, 'results/diff_genes_long_vs_normal.csv', row.names = FALSE)
write.csv(diffGenes.df.short.v.normal.gene.names, 'results/diff_genes_short_vs_normal.csv', row.names = FALSE)

# GSEA ----
myTopHits.long.v.short <- topTable(ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")
myTopHits.long.v.normal <- topTable(ebFit, adjust ="BH", coef=2, number=40000, sort.by="logFC")
myTopHits.short.v.normal <- topTable(ebFit, adjust ="BH", coef=3, number=40000, sort.by="logFC")

myTopHits.df.long.v.short <- myTopHits.long.v.short %>%
  as_tibble(rownames = "geneID")
myTopHits.df.long.v.normal <- myTopHits.long.v.normal %>%
  as_tibble(rownames = "geneID")
myTopHits.df.short.v.normal <- myTopHits.short.v.normal %>%
  as_tibble(rownames = "geneID")

hs_gsea_c2 <- msigdbr(species = "Caenorhabditis elegans", # change depending on species your data came from
                      category = "C2") %>% # choose your msigdb collection of interest
  dplyr::select(gs_name, ensembl_gene)

hs_gsea_c3 <- msigdbr(species = "Caenorhabditis elegans", # change depending on species your data came from
                      category = "C3") %>% # choose your msigdb collection of interest
  dplyr::select(gs_name, ensembl_gene)

c2_pathway_to_genes <- split(hs_gsea_c2$ensembl_gene, hs_gsea_c2$gs_name)
c3_pathway_to_genes <- split(hs_gsea_c3$ensembl_gene, hs_gsea_c3$gs_name)

# Taken from https://bioinformaticsbreakdown.com/how-to-gsea/
GSEA = function(gene_list, pathway_to_genes, pval) {
  set.seed(54321)
  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
  
  fgRes <- fgsea::fgsea(pathways = pathway_to_genes, 
                        stats = gene_list,
                        minSize=15,
                        maxSize=600) %>% 
    as.data.frame() %>% 
    dplyr::filter(padj < !!pval)
  
  ## Filter FGSEA by using gage results. Must be significant and in same direction to keep 
  gaRes = gage::gage(gene_list, gsets=pathway_to_genes, same.dir=TRUE, set.size =c(15,600))
  
  ups = as.data.frame(gaRes$greater) %>% 
    tibble::rownames_to_column("Pathway") %>% 
    dplyr::filter(!is.na(p.geomean) & q.val < pval ) %>%
    dplyr::select("Pathway")
  
  downs = as.data.frame(gaRes$less) %>% 
    tibble::rownames_to_column("Pathway") %>% 
    dplyr::filter(!is.na(p.geomean) & q.val < pval ) %>%
    dplyr::select("Pathway")
  
  keepups = data.table(fgRes[fgRes$NES > 0 & !is.na(match(fgRes$pathway, ups$Pathway)), ])
  keepdowns = data.table(fgRes[fgRes$NES < 0 & !is.na(match(fgRes$pathway, downs$Pathway)), ])
  
  fgRes <- fgRes %>%
    arrange(desc(NES))
  
  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  filtRes = rbind(head(fgRes, n = 15),
                  tail(fgRes, n = 15 ))
  g = ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_segment( aes(reorder(pathway, NES), xend=pathway, y=0, yend=NES)) +
    geom_point( size=5, aes( fill = Enrichment),
                shape=21, stroke=2) +
    scale_fill_manual(values = c("Down-regulated" = "dodgerblue",
                                 "Up-regulated" = "firebrick") ) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="GSEA - Biological Processes") + 
    theme_minimal()
  
  output = list("Results" = fgRes, "Plot" = g)
  return(output)
}

# divide by 3 for Bonferroni correction (3 pairwise comparisons)
res.long.v.short = GSEA(deframe(myTopHits.df.long.v.short[, 1:2]), c(c2_pathway_to_genes, c3_pathway_to_genes), pval = 0.05 / 3)
res.long.v.normal = GSEA(deframe(myTopHits.df.long.v.normal[, 1:2]), c(c2_pathway_to_genes, c3_pathway_to_genes), pval = 0.05 / 3)
res.short.v.normal = GSEA(deframe(myTopHits.df.short.v.normal[, 1:2]), c(c2_pathway_to_genes, c3_pathway_to_genes), pval = 0.05 / 3)


# Single experiment analysis ----
dir.create("results/single_experiments")
mart <- useEnsembl(biomart = "ensembl",
                   dataset = "celegans_gene_ensembl",
                   version = "100")

Tx.worm <- getBM(attributes=c('ensembl_transcript_id',
                              'ensembl_gene_id'),
                 mart = mart)

Tx.worm <- as_tibble(Tx.worm)
# we need to rename the two columns we just retrieved from biomart
Tx.worm <- dplyr::rename(Tx.worm, target_id = ensembl_transcript_id, 
                         gene_name = ensembl_gene_id)

for (dir in list.files("../query_download_process_data/ready_for_analysis/")) {
  if (!dir.exists(paste("../common_datastore/single_experiments/", dir, sep=''))) {
    path <- file.path("../query_download_process_data/ready_for_analysis", dir, list.files(paste("../query_download_process_data/ready_for_analysis/", dir, sep='')), 'abundance.tsv')
    
    Txi_gene <- tximport(path, 
                         type = "kallisto", 
                         tx2gene = Tx.worm, 
                         txOut = FALSE, # determines whether your data represented at transcript or gene level
                         ignoreTxVersion = FALSE)
    
    genes <- rownames(Txi_gene[["counts"]])
    samples = lapply(strsplit(path, "/"), '[', 5)
    colnames(Txi_gene[["counts"]]) <- samples
    dir.create(paste("../common_datastore/single_experiments/", dir, sep=''))
    file = paste(dir, '.csv', sep='')
    study_design_file <- labels_with_age[labels_with_age$Sample %in% samples, ]
    write.table(Txi_gene[["counts"]], file=paste('../common_datastore/single_experiments', 
                                                 dir, file, sep='/'), row.names=FALSE)
    write.table(study_design_file, file=paste('../common_datastore/single_experiments', dir, 'design_file.csv', sep='/')
                , row.names=FALSE)
  }
}

for (dir in list.files("../common_datastore/single_experiments/")) {
  gene_exp_data_file <- paste(dir, '.csv', sep='')
  gene_exp_data <- read.csv(paste('../common_datastore/single_experiments/', dir, '/', gene_exp_data_file ,sep=""), sep=' ')
  labels_with_age <- read.csv(paste('../common_datastore/single_experiments/', dir, '/', 'design_file.csv' ,sep=""), sep=' ')
  label_counts <- as.data.frame(table(labels_with_age$longevity))
  colnames(label_counts)[1] = 'Label'
  if (nrow(label_counts) == 1 | (1 %in% label_counts$Freq & nrow(label_counts) <= 2)) {
    next
  }
  rownames(gene_exp_data) <- genes
  
  myDGEList <- DGEList(gene_exp_data)
  cpm <- cpm(myDGEList)
  keepers <- rowSums(cpm>1)>=min(table(labels_with_age$longevity))
  myDGEList.filtered <- myDGEList[keepers,]
  
  myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
  
  group <- factor(labels_with_age$longevity)
  design <- model.matrix(~0 + group)
  
  if (all(c(0, 1, 2) %in% label_counts$Label)) {
    colnames(design) <- c('long', 'normal', 'short')
    variance_stabilized_adjusted <- voom(myDGEList.filtered.norm, design, plot = TRUE)
    fit <- lmFit(variance_stabilized_adjusted, design)
    contrast.matrix <- makeContrasts(long_vs_short = long - short,
                                     long_vs_normal = long - normal,
                                     short_vs_normal = short - normal,
                                     levels=design)
  }
  else if (all(c(0, 1) %in% label_counts$Label)) {
    colnames(design) <- c('long', 'normal')
    variance_stabilized_adjusted <- voom(myDGEList.filtered.norm, design, plot = TRUE)
    fit <- lmFit(variance_stabilized_adjusted, design)
    contrast.matrix <- makeContrasts(long_vs_normal = long - normal,
                                     levels=design)
  }
  else if (all(c(0, 2) %in% label_counts$Label)) {
    colnames(design) <- c('long', 'short')
    variance_stabilized_adjusted <- voom(myDGEList.filtered.norm, design, plot = TRUE)
    fit <- lmFit(variance_stabilized_adjusted, design)
    contrast.matrix <- makeContrasts(long_vs_short = long - short,
                                     levels=design)
  }
  else if (all(c(1, 2) %in% label_counts$Label)) {
    colnames(design) <- c('normal', 'short')
    variance_stabilized_adjusted <- voom(myDGEList.filtered.norm, design, plot = TRUE)
    fit <- lmFit(variance_stabilized_adjusted, design)
    contrast.matrix <- makeContrasts(short_vs_normal = short - normal,
                                     levels=design)
  }
  fits <- contrasts.fit(fit, contrast.matrix)
  ebFit <- eBayes(fits)
  results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.01, lfc=2)
  savefile <- 'test_results_matrix.csv'
  dir.create(paste('results/single_experiments/', dir, sep=''))
  write.csv(results, file=paste("results", "single_experiments", dir, savefile, sep='/'))
  
  if (all(c(0, 1, 2) %in% label_counts$Label)) {
    diffGenes_long_v_short <- data.frame(variance_stabilized_adjusted$E)[results[,1] != 0,]
    diffGenes_long_v_normal <- data.frame(variance_stabilized_adjusted$E)[results[,2] != 0,]
    diffGenes_short_v_normal <- data.frame(variance_stabilized_adjusted$E)[results[,3] != 0,]
    write.csv(diffGenes_long_v_short, file=paste("results", "single_experiments", dir,
                                                 'diffGenes_long_v_short.csv', sep='/'))
    write.csv(diffGenes_long_v_normal, file=paste("results", "single_experiments", dir,
                                                  'diffGenes_long_v_normal.csv', sep='/'))
    write.csv(diffGenes_short_v_normal, file=paste("results", "single_experiments", dir,
                                                   'diffGenes_short_v_normal.csv', sep='/'))
  }
  else if (all(c(0, 1) %in% label_counts$Label)) {
    diffGenes_long_v_normal <- data.frame(variance_stabilized_adjusted$E)[results[,1] != 0,]
    write.csv(diffGenes_long_v_normal, file=paste("results", "single_experiments", dir,
                                                  'diffGenes_long_v_normal.csv', sep='/'))
  }
  else if (all(c(0, 2) %in% label_counts$Label)) {
    diffGenes_long_v_short <- data.frame(variance_stabilized_adjusted$E)[results[,1] != 0,]
    write.csv(diffGenes_long_v_short, file=paste("results", "single_experiments", dir,
                                                 'diffGenes_long_v_short.csv', sep='/'))
  }
  else if (all(c(1, 2) %in% label_counts$Label)) {
    diffGenes_short_v_normal <- data.frame(variance_stabilized_adjusted$E)[results[,1] != 0,]
    write.csv(diffGenes_short_v_normal, file=paste("results", "single_experiments", dir,
                                                   'diffGenes_short_v_normal.csv', sep='/'))
  }
}

differentially_expressed_genes <- tibble()
colnames(differentially_expressed_genes) <- c('Bioproject', 'Gene', 'Comparison')

for (dir in list.files("../common_datastore/single_experiments/")) {
  for (file in list.files(paste("results", "single_experiments", dir, sep = '/'))) {
    if (file == 'diffGenes_long_v_normal.csv') {
      diff_genes <- read.csv(paste("results", "single_experiments", dir, "diffGenes_long_v_normal.csv",
                                   sep='/'))
      if(dim(diff_genes)[1] != 0) {
        colnames(diff_genes) <- c('Gene', colnames(diff_genes[2:length(colnames(diff_genes))]))
        diff_genes <- diff_genes %>%
          dplyr::select(Gene) %>%
          add_column(Bioproject=dir) %>%
          add_column(Comparison='long vs normal')
      }
    }
    else if (file == 'diffGenes_long_v_short.csv') {
      diff_genes <- read.csv(paste("results", "single_experiments", dir, "diffGenes_long_v_short.csv",
                                   sep='/'))
      if(dim(diff_genes)[1] != 0) {
        colnames(diff_genes) <- c('Gene', colnames(diff_genes[2:length(colnames(diff_genes))]))
        diff_genes <- diff_genes %>%
          dplyr::select(Gene) %>%
          add_column(Bioproject=dir) %>%
          add_column(Comparison='long vs short')
      }
    }
    else if (file == 'diffGenes_short_v_normal.csv') {
      diff_genes <- read.csv(paste("results", "single_experiments", dir, "diffGenes_short_v_normal.csv",
                                   sep='/'))
      if(dim(diff_genes)[1] != 0) {
        colnames(diff_genes) <- c('Gene', colnames(diff_genes[2:length(colnames(diff_genes))]))
        diff_genes <- diff_genes %>%
          dplyr::select(Gene) %>%
          add_column(Bioproject=dir) %>%
          add_column(Comparison='short vs normal')
      }
    }
    if (dim(diff_genes)[1] != 0) {
      differentially_expressed_genes <- bind_rows(differentially_expressed_genes, diff_genes)
    }
  }
}

differentially_expressed_genes <- distinct(differentially_expressed_genes)
comparison_counts <- differentially_expressed_genes %>% 
  group_by(Comparison) %>%
  distinct(Bioproject) %>%
  count()

colnames(comparison_counts) <- c('Comparison', 'Comparison_n')

differentially_expressed_genes <- differentially_expressed_genes %>%
  group_by(Comparison) %>%
  count(Gene)

differentially_expressed_genes_with_perc <- 
  inner_join(differentially_expressed_genes, comparison_counts) %>%
  mutate(percent = n / Comparison_n)

percent_25_genes <- differentially_expressed_genes_with_perc %>%
  filter(percent >= 0.25)

write.csv(percent_25_genes$Gene, 'results/single_experiments/percent_25_wb_ids.csv')

# Manually obtain names from https://wormbase.org/tools/mine/simplemine.cgi with Public Name and 
# "keep duplicate gene entries in results"
gene_names <- read.table('results/single_experiments/genes_25_percent.txt', sep='\t', header=TRUE)
percent_25_genes <- percent_25_genes %>%
  add_column(gene_names[2], .before=1)
colnames(percent_25_genes) <- c('Gene name', 'Comparison', 'WB_ID', 'n', 'Comparison_n', 'Percent')
percent_25_genes_long_vs_normal <- percent_25_genes %>%
  filter(Comparison == 'long vs normal')
write.table(percent_25_genes, row.names = FALSE, 'results/single_experiments/percent_25_gene_with_percentages.csv', )
write.table(percent_25_genes_long_vs_normal, row.names = FALSE, 'results/single_experiments/percent_25_gene_long_vs_normal.csv', )


percent_50_genes <- differentially_expressed_genes_with_perc %>%
  filter(percent >= 0.50)

write.csv(percent_50_genes$Gene, 'results/single_experiments/percent_50_wb_ids.csv')

# Manually obtain names
gene_names <- read.table('results/single_experiments/genes_50_percent.txt', sep='\t', header=TRUE)
percent_50_genes <- percent_50_genes %>%
  add_column(gene_names[2], .before=1)
colnames(percent_50_genes) <- c('Gene name', 'Comparison', 'WB_ID', 'n', 'Comparison_n', 'Percent')
write.table(percent_50_genes, row.names = FALSE, 'results/single_experiments/percent_50_gene_with_percentages.csv')
