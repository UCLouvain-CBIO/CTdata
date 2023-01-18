## Code to prepare `DAC_treated_cells_multimapping` dataset goes here

library(tidyverse)
library(DESeq2)
load("../../eh_data/GTEX_data.rda")
genes_in_gtex <- rownames(GTEX_data)

## RNAseq data from cells treated or not with 5-aza downloaded from SRA.
## Data was processed using a standard RNAseq pipeline including
## [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for the
## quality control of the raw data, and
## [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
## to remove low quality reads and trim the adapter from the sequences.
## [hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml) was used to align
## reads to grch38 genome.
## [featurecounts](https://rdrr.io/bioc/Rsubread/man/featureCounts.html) was
## used to assign reads to genes using Homo_sapiens.GRCh38.94.gtf.
## To allow for multi-mapped reads to be counted,
## hisat2 was run with -k 20 parameter (to report up to 20 alignments per read),
## and featurecounts was run with -M parameter (multi-mapping reads are counted).

load("../extdata/DAC_coldata.rda")
load("../extdata/DAC_raw_counts_multiM.rda")

dds <- DESeqDataSetFromMatrix(countData = raw_counts_with_MP,
                              colData = coldata,
                              design = ~ 1)
dds <- DESeq2::estimateSizeFactors(dds)
log1p_transformed <- log1p(counts(dds, normalize = TRUE))

## DESeq2 analysis
## Compare DAC treated cells to control cells separately each cell line
## to identify genes induced by DAC in at least one cell line.
cell_line <- unique(coldata$cell)[1]
coldata_by_cell_line <- coldata[coldata$cell == cell_line, ]
dds <- DESeqDataSetFromMatrix(
  countData = raw_counts_with_MP[, coldata_by_cell_line$sample],
  colData = coldata_by_cell_line,
  design = ~ treatment)
dds <- DESeq(dds)

res <- results(dds, name = "treatment_DAC_vs_CTL",
               independentFiltering = TRUE,
               altHypothesis = "greater")

res_all <- as_tibble(res, rownames = "ensembl_gene_id") %>%
  right_join(as_tibble(rowData(GTEX_data), rownames = "ensembl_gene_id") %>%
               dplyr::select(ensembl_gene_id, external_gene_name)) %>%
  dplyr::select(ensembl_gene_id, external_gene_name, log2FoldChange, padj) %>%
  mutate(log2FoldChange = round(log2FoldChange, 2)) %>%
  mutate(sign = case_when((!is.na(log2FoldChange) & log2FoldChange >= 2 &
                             !is.na(padj) & padj <= 0.05) ~ 1,
                          (is.na(log2FoldChange) | log2FoldChange < 2 |
                             is.na(padj) | padj > 0.05) ~ 0))

names(res_all) <- c("ensembl_gene_id", "external_gene_name",
                    paste0("logFC_", cell_line),
                    paste0("padj_", cell_line),
                    paste0("sign_", cell_line))

for(cell_line in unique(coldata$cell)[-1]) {

  coldata_by_cell_line <- coldata[coldata$cell == cell_line, ]
  dds <- DESeqDataSetFromMatrix(
    countData = raw_counts_with_MP[, coldata_by_cell_line$sample],
    colData = coldata_by_cell_line,
    design = ~ treatment)
  dds <- DESeq(dds)

  res <- results(dds, name = "treatment_DAC_vs_CTL",
                 independentFiltering = TRUE,
                 altHypothesis = "greater")
  res <- as_tibble(res, rownames = "ensembl_gene_id") %>%
    right_join(as_tibble(rowData(GTEX_data), rownames = "ensembl_gene_id") %>%
                 dplyr::select(ensembl_gene_id, external_gene_name)) %>%
    dplyr::select(ensembl_gene_id, external_gene_name, log2FoldChange, padj) %>%
    mutate(log2FoldChange = round(log2FoldChange, 2)) %>%
    mutate(sign = case_when((!is.na(log2FoldChange) & log2FoldChange >= 2 &
                               !is.na(padj) & padj <= 0.05) ~ 1,
                            (is.na(log2FoldChange) | log2FoldChange < 2 |
                               is.na(padj) | padj > 0.05) ~ 0))

  names(res) <- c("ensembl_gene_id", "external_gene_name",
                  paste0("logFC_", cell_line),
                  paste0("padj_", cell_line),
                  paste0("sign_", cell_line))

  res_all <- left_join(res_all, res)
}

## Tag genes significantly induced in at least one cell line
res_all$sign <- res_all %>%
  dplyr::select(starts_with("sign")) %>%
  rowSums()

res_all <- res_all %>%
  mutate(induced = ifelse(sign >= 1, TRUE, FALSE)) %>%
  dplyr::select(-starts_with("sign"))

res_all <- as.data.frame(res_all)
res_all <- column_to_rownames(res_all, "ensembl_gene_id")

## Create a summarizedExperiment
DAC_treated_cells_multimapping <- SummarizedExperiment(
  assays = list(log1p = log1p_transformed[genes_in_gtex, coldata$sample]),
  colData = coldata[, -1],
  rowData = res_all[genes_in_gtex, ])

save(DAC_treated_cells_multimapping, file = "../../eh_data/DAC_treated_cells_multimapping.rda",
     compress = "xz",
     compression_level = 9)
