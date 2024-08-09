## Code to prepare `hESC_data` dataset goes here

library("tidyverse")
library("SummarizedExperiment")
load("../../eh_data/all_genes.rda")


## RNAseq data from embryonic stem cells were downloaded from Encode.
## H1 (total RNA)  (ENCFF892WVN.fastq.gz, ENCFF481BLH.fastq.gz)
## H7 (total RNA)  (ENCFF002DMP.fastq.gz, ENCFF002DMQ.fastq.gz)
## H9 (total RNA)  (ENCFF778LAE.fastq.gz, ENCFF650HFX.fastq.gz)
## HUES64 (polyA) (ENCFF098DFE.fastq.gz, ENCFF608FBD.fastq.gz)

## NB : No WGBS data for H7 but keeping it in expression anyway

## Data was processed using a standard RNAseq pipeline including
## [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for the
## quality control of the raw data, and
## [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
## to remove low quality reads and trim the adapter from the sequences.
## [hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml) was used to align
## reads to grch38 genome.
## [featurecounts](https://rdrr.io/bioc/Rsubread/man/featureCounts.html) was
## used to assign reads to genes using Homo_sapiens.GRCh38.105.gtf.

load("../../../CTdata_extdata/hESC_RNAseq_coldata.rda")
load("../../../CTdata_extdata/hESC_RNAseq_raw_counts.rda")

################################################################################
## Coldata creation
################################################################################


coldata <- as.data.frame(coldata) %>%
  mutate(cell_type = "hESC") %>%
  mutate(genotype = factor(coldata$genotype, levels = c("XX", "XY"))) %>%
  mutate(dataset = "Encode") %>%
  mutate(data = "Expression")

rownames(coldata) <- coldata$sample

################################################################################
## TPM normalisation
################################################################################

gene_lengths <-
  read_table("../../../CTdata_extdata/H1_featurecounts.tsv",
             skip = 1) %>%
  dplyr::select(Geneid, Length)
x1 <- counts / gene_lengths$Length * 1000
total <- colSums(x1)

TPM_values <- as_tibble(x1, rownames = "Geneid") %>%
  pivot_longer(names_to = "sample", values_to = "counts", -Geneid) %>%
  left_join(enframe(total) %>%
              dplyr::rename(sample = name, total = value)) %>%
  mutate(TPM = round(counts / total * 1000000, 2)) %>%
  dplyr::select(Geneid, sample, TPM) %>%
  pivot_wider(names_from = sample, values_from = TPM) %>%
  dplyr::rename(ensembl_gene_id = Geneid) %>%
  filter(ensembl_gene_id %in% all_genes$ensembl_gene_id)



################################################################################
## Create hESC_data SE
################################################################################

rowdata <- TPM_values %>%
  dplyr::select(ensembl_gene_id) %>%
  left_join(dplyr::select(all_genes, ensembl_gene_id, external_gene_name)) %>%
  column_to_rownames("ensembl_gene_id")


hESC_data <- SummarizedExperiment(assays = column_to_rownames(TPM_values,
                                                              "ensembl_gene_id"),
                                  colData = coldata[, -1],
                                  rowData = rowdata)

save(hESC_data,
     file = "../../eh_data/hESC_data.rda",
     compress = "xz",
     compression_level = 9)

