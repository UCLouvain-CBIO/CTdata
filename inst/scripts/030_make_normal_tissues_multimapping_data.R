## Code to prepare `normal_tissues_multimapping_data` dataset goes here

library("tidyverse")
library("SummarizedExperiment")
load("../../eh_data/GTEX_data.rda")

## RNAseq data from a set of normal tissues were downloaded from Encode.
## testis (ENCFF140UYT.fastq.gz, ENCFF794EAB.fastq.gz)
## thyroid_gland (ENCFF151GUG.fastq.gz, ENCFF628TMU.fastq.gz)
## gastrocnemius_medialis (ENCFF004CNM.fastq.gz, ENCFF086LCO.fastq.gz)
## adrenal_gland (ENCFF911BTP.fastq.gz, ENCFF904VHM.fastq.gz)
## thoracic_aorta (ENCFF661PVB.fastq.gz, ENCFF152GUK.fastq.gz)
## subcutaneous_adipose_tissue (ENCFF667CWY.fastq.gz, ENCFF360KZB.fastq.gz)
## stomach (ENCFF741NGG.fastq.gz, ENCFF582ILA.fastq.gz)
## upper_lobe_of_left_lung (ENCFF719YBM.fastq.gz, ENCFF801ZKX.fastq.gz)
## suprapubic_skin (ENCFF398KGB.fastq.gz, ENCFF058JYK.fastq.gz)
## breast_epithelium (ENCFF767QVV.fastq.gz, ENCFF050GYP.fastq.gz)
## lower_leg_skin (ENCFF431RAQ.fastq.gz, ENCFF008OVI.fastq.gz)
## sigmoid_colon (ENCFF153ULW.fastq.gz, ENCFF182OWD.fastq.gz)
## transverse_colon (ENCFF411UIT.fastq.gz, ENCFF992NAN.fastq.gz)
## spleen (ENCFF567ORO.fastq.gz, ENCFF743QDT.fastq.gz)
## gastroesophageal_sphincter (ENCFF250BPC.fastq.gz, ENCFF842DCO.fastq.gz)
## tibial_nerve (ENCFF534AYT.fastq.gz, ENCFF935LBC.fastq.gz)
## esophagus_muscularis_mucosa (ENCFF091UZU.fastq.gz, ENCFF163DLM.fastq.gz)
## omental_fat_pad (ENCFF745PTG.fastq.gz, ENCFF824ZLA.fastq.gz)
## esophagus_squamous_epithelium (ENCFF585TOW.fastq.gz, ENCFF230GLF.fastq.gz)

## Data was processed using a standard RNAseq pipeline including
## [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for the
## quality control of the raw data, and
## [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
## to remove low quality reads and trim the adapter from the sequences.
## [hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml) was used to align
## reads to grch38 genome.
## [featurecounts](https://rdrr.io/bioc/Rsubread/man/featureCounts.html) was
## used to assign reads to genes using Homo_sapiens.GRCh38.94.gtf.

## Two different pipelines were run in order to remove or not multi-mapping reads.
## When multi-mapping was allowed, hisat2 was run with -k 20 parameter (reports
## up to 20 alignments per read), and featurecounts was run with -M parameter
## (multi-mapping reads are counted).

################################################################################
## Data generated without allowing multi-mapping
################################################################################

load("../extdata/normal_tissues_RNAseq_coldata.rda")
load("../extdata/normal_tissues_RNAseq_raw_counts.rda")

## TPM normalisation
## Keep only genes present in GTEx database
gene_lengths <-
  read_table("../extdata/adrenal_gland_featurecounts.tsv",
             skip = 1) %>%
  dplyr::select(Geneid, Length)
x1 <- raw_counts / gene_lengths$Length * 1000
total <- colSums(x1)

TPM_matrix_no_multimapping <- as_tibble(x1, rownames = "Geneid") %>%
  gather(sample, counts, -Geneid) %>%
  left_join(enframe(total) %>%
              dplyr::rename(sample = name, total = value)) %>%
  mutate(TPM = round(counts / total * 1000000, 2)) %>%
  dplyr::select(Geneid, sample, TPM) %>%
  spread(sample, TPM) %>%
  dplyr::rename(ensembl_gene_id = Geneid) %>%
  right_join(as_tibble(rowData(GTEX_data), rownames = "ensembl_gene_id") %>%
               dplyr::select(ensembl_gene_id, external_gene_name)) %>%
  dplyr::select(ensembl_gene_id, external_gene_name, everything())

mat_no_multimapping <- as.matrix(TPM_matrix_no_multimapping[, -c(1:2)])
rownames(mat_no_multimapping) <- TPM_matrix_no_multimapping$ensembl_gene_id

################################################################################
## Data generated when allowing multi-mapping
################################################################################

load(file = "../extdata/normal_tissues_RNAseq_raw_counts_multiM.rda")

## Normalise in TPM
## Keep only genes present in GTEx database
x1 <- raw_counts_with_MP / gene_lengths$Length * 1000
total <- colSums(x1)

TPM_matrix_with_multimapping <- as_tibble(x1, rownames = "Geneid") %>%
  gather(sample, counts, -Geneid) %>%
  left_join(enframe(total) %>%
              dplyr::rename(sample = name, total = value)) %>%
  mutate(TPM = round(counts / total * 1000000, 2)) %>%
  dplyr::select(Geneid, sample, TPM) %>%
  spread(sample, TPM) %>%
  dplyr::rename(ensembl_gene_id = Geneid) %>%
  right_join(as_tibble(rowData(GTEX_data), rownames = "ensembl_gene_id") %>%
               dplyr::select(ensembl_gene_id, external_gene_name)) %>%
  dplyr::select(ensembl_gene_id, external_gene_name, everything())

mat_with_multimapping <- as.matrix(TPM_matrix_with_multimapping[, -c(1:2)])
rownames(mat_with_multimapping) <- TPM_matrix_with_multimapping$ensembl_gene_id

################################################################################
## Assess testis-specificity of genes flagged as "lowly_expressed" in
## GTEX_category from GTEX_data (when their TPM < 1 in testis)
################################################################################
## Many Cancer-Testis genes belong to gene families from which members
## have identical or nearly identical sequences. This is likely the reason
## why these genes are not detected in GTEx database, as GTEX processing
## pipeline specifies that overlapping intervals between genes are excluded
## from all genes for counting.
## Testis specificity of genes flagged as "lowly_expressed" in GTEX_category
## could however be assessed by comparing expression values obtained by
## counting or not multi-mapped reads in a set of normal tissues.
## Some of these genes are indeed only detected when multimapping reads are
## not discarded.
## Genes flagged as "testis_specific" in `multimapping_analysis` column must be
## detectable in testis (TPM >= 1) when multimapping is allowed, their TPM
## value must have increased when multimapping is allowed (ratio > 5), and
## these genes must also have a TPM value (obtained by allowing multimapping)
## at least 10 times higher in testis than in any other somatic tissue (where
## the maximum expression always has to be below 1 TPM).

ratio_multi_not_multi <- as_tibble(mat_no_multimapping,
                                   rownames = "ensembl_gene_id") %>%
  dplyr::select(ensembl_gene_id, testis) %>%
  dplyr::rename(TPM_testis_no_multi = testis) %>%
  left_join(as_tibble(mat_with_multimapping, rownames = "ensembl_gene_id") %>%
              dplyr::select(ensembl_gene_id, testis) %>%
              dplyr::rename(TPM_testis_when_multi = testis)) %>%
  mutate(ratio = TPM_testis_when_multi / TPM_testis_no_multi)

max_in_somatic <- rowMax(
  mat_with_multimapping[, -which(colnames(mat_with_multimapping) == "testis")])
max_not_testis <- tibble(ensembl_gene_id = rownames(mat_with_multimapping),
                         max_in_somatic = max_in_somatic)

genes_testis_specific_in_multimapping <- ratio_multi_not_multi %>%
  left_join(max_not_testis) %>%
  mutate(ratio_testis_other = TPM_testis_when_multi / max_in_somatic) %>%
  left_join(as_tibble(rowData(GTEX_data), rownames = "ensembl_gene_id")) %>%
  mutate(multimapping_analysis = case_when(
    GTEX_category != "lowly_expressed" ~ "not_analysed",
    GTEX_category == "lowly_expressed" & TPM_testis_when_multi >= 1 &
      ratio >= 5 & ratio_testis_other >= 10 & max_in_somatic <= 1 ~ "testis_specific",
    GTEX_category == "lowly_expressed" &
      (TPM_testis_when_multi < 1 | ratio < 5 | ratio_testis_other < 10) ~
      "not_testis_specific"))
rowdata <-
  tibble(ensembl_gene_id = TPM_matrix_with_multimapping$ensembl_gene_id,
         external_gene_name =
           TPM_matrix_with_multimapping$external_gene_name) %>%
  left_join(genes_testis_specific_in_multimapping %>%
              dplyr::select(ensembl_gene_id, GTEX_category,
                            multimapping_analysis)) %>%
  mutate(lowly_expressed_in_GTEX = case_when(
    GTEX_category == "lowly_expressed" ~ TRUE,
    GTEX_category != "lowly_expressed" ~ FALSE)) %>%
  dplyr::select(ensembl_gene_id, external_gene_name,
                lowly_expressed_in_GTEX, multimapping_analysis)

rowdata$multimapping_analysis[rowdata$multimapping_analysis == "not_analysed"] <- NA

rowdata <- column_to_rownames(rowdata, "ensembl_gene_id")

################################################################################
## Create normal_tissues_multimapping SE
################################################################################

normal_tissues_multimapping_data <- SummarizedExperiment(
  assays = list(TPM_no_multimapping = mat_no_multimapping,
                TPM_with_multimapping = mat_with_multimapping),
  rowData = rowdata)

save(normal_tissues_multimapping_data,
     file = "../../eh_data/normal_tissues_multimapping_data.rda",
     compress = "xz",
     compression_level = 9)

