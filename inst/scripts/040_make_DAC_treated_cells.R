## Code to prepare `DAC_treated_cells` dataset goes here

library(tidyverse)
library(DESeq2)
load("../../eh_data/GTEX_data.rda")
genes_in_gtex <- rownames(GTEX_data)

## RNAseq data from cells treated or not with 5-aza downloaded from SRA
## SRR3326020/SRR3326021	IMR5-75	CTL
## SRR3326022/SRR3326023	IMR5-75	DAC
## SRR9108737/SRR9108738	HCT116	CTL
## SRR9108739/SRR9108740	HCT116	DAC
## SRR1618781/SRR1618782	HEK293T	CTL
## SRR1618783/SRR1618784	HEK293T	DAC
## SRR3362409/SRR3362410	HMLER	DAC
## SRR3362411/SRR3362412	HMLER	CTL
## SRR12105788/SRR12105789	NCH612	CTL
## SRR12105790/SRR12105791	NCH612	DAC
## SRR12105792/SRR12105793	NCH1681	CTL
## SRR12105794/SRR12105795	NCH1681	DAC
## SRR12105780/SRR12105781	TS603	CTL
## SRR12105782/SRR12105783	TS603	DAC
## SRR5363797/SRR5363798	B2-1	CTL
## SRR5363799/SRR5363800	B2-1	DAC

## Data was processed using a standard RNAseq pipeline including
## [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for the
## quality control of the raw data, and
## [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
## to remove low quality reads and trim the adapter from the sequences.
## [hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml) was used to align
## reads to grch38 genome.
## [featurecounts](https://rdrr.io/bioc/Rsubread/man/featureCounts.html) was
## used to assign reads to genes using Homo_sapiens.GRCh38.105.gtf.

load("../../../CTdata_extdata/DAC_coldata.rda")
load("../../../CTdata_extdata/DAC_raw_counts.rda")

dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = coldata,
                              design = ~ 1)
dds <- DESeq2::estimateSizeFactors(dds)
normalised_counts <- counts(dds, normalize = TRUE)
log1p_transformed <- log1p(counts(dds, normalize = TRUE))

## DESeq2 analysis
## Compare DAC treated cells to control cells separately each cell line
## to identify genes induced by DAC in at least one cell line.

padj_thr <- 0.1
logFC_thr <- 2

## Define an expression threshold above which a gene will be considered as
## already activated in the cell line => Tag cells in which it makes no sense
## to test the activation upon DAC treatment
thr <- 15

cell_line <- unique(coldata$cell)[1]
coldata_by_cell_line <- coldata[coldata$cell == cell_line, ]
dds <- DESeqDataSetFromMatrix(
  countData = raw_counts[, coldata_by_cell_line$sample],
  colData = coldata_by_cell_line,
  design = ~ treatment)
dds <- DESeq(dds)

res <- results(dds, name = "treatment_DAC_vs_CTL",
               independentFiltering = FALSE,
               altHypothesis = "greater")

means_in_CTL <- enframe(rowMeans(counts(dds[, dds$treatment == "CTL"],
                                        normalize = TRUE)),
                        name = "ensembl_gene_id", value = "mean_in_CTL")

res_all <- as_tibble(res, rownames = "ensembl_gene_id") %>%
  right_join(as_tibble(rowData(GTEX_data), rownames = "ensembl_gene_id") %>%
               dplyr::select(ensembl_gene_id, external_gene_name)) %>%
  dplyr::select(ensembl_gene_id, external_gene_name, log2FoldChange,
                pvalue, padj) %>%
  mutate(log2FoldChange = round(log2FoldChange, 2)) %>%
  left_join(means_in_CTL) %>%
  mutate(expressed_in_CTL = ifelse(mean_in_CTL > thr, TRUE, FALSE)) %>%
  mutate(sign = case_when((!expressed_in_CTL & !is.na(log2FoldChange) &
                             log2FoldChange >= logFC_thr &
                             !is.na(padj) & padj <= padj_thr) ~ 1,
                          (is.na(log2FoldChange) | log2FoldChange < logFC_thr |
                             is.na(padj) | padj > padj_thr) ~ 0))

names(res_all) <- c("ensembl_gene_id", "external_gene_name",
                    paste0("logFC_", cell_line),
                    paste0("pval_", cell_line),
                    paste0("padj_", cell_line),
                    paste0("mean_CTL_", cell_line),
                    paste0("expressed_in_CTL_", cell_line),
                    paste0("sign_", cell_line))

for(cell_line in unique(coldata$cell)[-1]) {

  coldata_by_cell_line <- coldata[coldata$cell == cell_line, ]
  dds <- DESeqDataSetFromMatrix(
    countData = raw_counts[, coldata_by_cell_line$sample],
    colData = coldata_by_cell_line,
    design = ~ treatment)
  dds <- DESeq(dds)

  res <- results(dds, name = "treatment_DAC_vs_CTL",
                 independentFiltering = FALSE,
                 altHypothesis = "greater")

  means_in_CTL <- enframe(rowMeans(counts(dds[, dds$treatment == "CTL"],
                                          normalize = TRUE)),
                          name = "ensembl_gene_id", value = "mean_in_CTL")

  res <- as_tibble(res, rownames = "ensembl_gene_id") %>%
    right_join(as_tibble(rowData(GTEX_data), rownames = "ensembl_gene_id") %>%
                 dplyr::select(ensembl_gene_id, external_gene_name)) %>%
    dplyr::select(ensembl_gene_id, external_gene_name, log2FoldChange,
                  pvalue, padj) %>%
    mutate(log2FoldChange = round(log2FoldChange, 2)) %>%
    left_join(means_in_CTL) %>%
    mutate(expressed_in_CTL = ifelse(mean_in_CTL > thr, TRUE, FALSE)) %>%
    mutate(sign = case_when((!expressed_in_CTL & !is.na(log2FoldChange) &
                               log2FoldChange >= logFC_thr &
                               !is.na(padj) & padj <= padj_thr) ~ 1,
                            (is.na(log2FoldChange) |
                               log2FoldChange < logFC_thr |
                               is.na(padj) | padj > padj_thr) ~ 0))

  names(res) <- c("ensembl_gene_id", "external_gene_name",
                  paste0("logFC_", cell_line),
                  paste0("pval_", cell_line),
                  paste0("padj_", cell_line),
                  paste0("mean_CTL_", cell_line),
                  paste0("expressed_in_CTL_", cell_line),
                  paste0("sign_", cell_line))

  res_all <- left_join(res_all, res)
}



## Evaluate the number de cell lines already expressing the genes
res_all$expressed_in_n_CTLs <- res_all %>%
  dplyr::select(starts_with("expressed_in_CTL")) %>%
  rowSums(na.rm = TRUE)

## Tag genes significantly induced in at least one cell line
res_all$sign <- res_all %>%
  dplyr::select(starts_with("sign")) %>%
  rowSums(na.rm = TRUE)

## Induction cannot be evaluate if genes are already expressed in more
## than 7 out of 8 cell lines
res_all <- res_all %>%
  mutate(induced = case_when(expressed_in_n_CTLs == 8 ~ NA,
                             expressed_in_n_CTLs < 8 & sign >= 1 ~ TRUE,
                             expressed_in_n_CTLs < 8 & sign < 1 ~ FALSE)) %>%
  dplyr::select(-sign)

res_all <- as.data.frame(res_all)
res_all <- column_to_rownames(res_all, "ensembl_gene_id")

## Save as a SE
DAC_treated_cells <- SummarizedExperiment(
  assays = list(counts = normalised_counts[genes_in_gtex, coldata$sample],
                log1p = log1p_transformed[genes_in_gtex, coldata$sample]),
  colData = coldata[, -1],
  rowData = res_all[genes_in_gtex, ])

save(DAC_treated_cells, file = "../../eh_data/DAC_treated_cells.rda",
     compress = "xz",
     compression_level = 9)
