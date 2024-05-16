## code to prepare `GTEX_data` dataset goes here

library(tidyverse)
library(SummarizedExperiment)
library(BiocFileCache)
library(biomaRt)

bfc <- BiocFileCache(cache = "../BiocFileCache",
                     ask = FALSE)

if (length(bfcquery(bfc, "GTEX")$rid) == 0) {
  url <- "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
  bfcadd(bfc, "GTEX", fpath = url)
}

GTEX_data <- as_tibble(read.table(bfc[[bfcquery(bfc, "GTEX")$rid]],
                                  skip = 2, header = TRUE, sep = "\t"))

##########################################################################
## Clean and simplify GTEX data
##########################################################################

## Add external_gene_name to GTEX database based on ensembl_gene_id
## (Gene named given in `Description` column is not always external_gene_name)
## ! external_gene_names are sometimes associated to several ensembl_gene_id !
## To keep the more relevant external_gene_names, select the ones associated
## to canonical transcripts
ensembl <- biomaRt::useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
attributes_vector <- c("ensembl_gene_id",
                       "external_gene_name",
                       "ensembl_transcript_id",
                       "external_transcript_name",
                       "chromosome_name",
                       "transcript_biotype",
                       "transcript_is_canonical")
transcripts_infos <- as_tibble(biomaRt::getBM(attributes = attributes_vector,
                                              mart = ensembl))
canonical_transcripts <- transcripts_infos %>%
  dplyr::filter(transcript_is_canonical == 1) %>%
  dplyr::filter(external_transcript_name != "") %>%
  dplyr::filter(chromosome_name %in% c(1:22, "X", "Y", "MT")) %>%
  # no transcripts on assembly patches
  dplyr::filter(transcript_biotype == "protein_coding" |
                  transcript_biotype == "lncRNA")

GTEX_data <- GTEX_data %>%
  dplyr::select(-Description) %>%
  dplyr::rename(ensembl_gene_id = Name) %>%
  dplyr::mutate(ensembl_gene_id = sub(pattern = ".\\d+$",
                                      x = ensembl_gene_id,
                                      replacement = '')) %>%
  left_join(canonical_transcripts %>%
              dplyr::select(ensembl_gene_id, external_gene_name)) %>%
  dplyr::filter(!is.na(external_gene_name))

## Clean tissue names
## Pool same categories of tissues and calculate mean of TPM in pooled tissues
GTEX_data <- GTEX_data %>%
  dplyr::rename("Small Intestine...Terminal Ileum" =
                  "Small.Intestine...Terminal.Ileum") %>%
  dplyr::rename("Salivary Gland" = "Minor.Salivary.Gland") %>%
  dplyr::rename("Adrenal Gland" = "Adrenal.Gland") %>%
  dplyr::rename("Blood" = "Whole.Blood") %>%
  dplyr::rename("Fallopian Tube" = "Fallopian.Tube") %>%
  dplyr::rename("Fibroblasts" = "Cells...Cultured.fibroblasts") %>%
  dplyr::rename("Lymphocytes_EBV" = "Cells...EBV.transformed.lymphocytes") %>%
  pivot_longer(names_to = "tissue", values_to = "tpm",
               -c(ensembl_gene_id, external_gene_name)) %>%
  dplyr::mutate(tissue = sub(pattern = "\\.\\.\\..*$",
                             x = tissue, replacement = '')) %>%
  dplyr::group_by(ensembl_gene_id, tissue) %>%
  dplyr::mutate(TPM = mean(tpm)) %>%
  dplyr::select(-tpm) %>%
  unique() %>%
  pivot_wider(names_from = tissue, values_from = TPM)

ordered_tissues <- sort(colnames(GTEX_data)[-c(1:2)])
GTEX_data <- GTEX_data %>%
  dplyr::select("ensembl_gene_id", "external_gene_name",
                all_of(ordered_tissues))

##########################################################################
## Evaluate maximum TPM in somatic tissues
## Evaluate quantile 75% of TPM in somatic tissues
## Calculate ratio_testis_somatic
## (defined as expression in testis / max_TPM_somatic)
##########################################################################

GTEX_data <- GTEX_data %>%
  dplyr::select("ensembl_gene_id", "external_gene_name", "Testis",
                "Ovary", everything()) %>%
  rowwise() %>%
  dplyr::mutate(max_TPM_somatic = max(c_across(c(Adipose:Vagina)))) %>%
  dplyr::mutate(q75_TPM_somatic = quantile(c_across(c(Adipose:Vagina)),
                                           0.75)) %>%
  dplyr::mutate(ratio_testis_somatic = Testis / max_TPM_somatic) %>%
  ungroup()

## Tag genes that are not detected in GTEx database
## (or genes very lowly expressed in all tissues)
lowly_expressed <- GTEX_data %>%
  dplyr::filter(Testis < 1 & max_TPM_somatic < 1 & Ovary < 1) %>%
  pull(external_gene_name)

## Tag testis_specific genes
testis_specific <- GTEX_data %>%
  dplyr::filter(!external_gene_name %in% lowly_expressed) %>%
  dplyr::filter(Testis >= 1) %>%
  dplyr::filter(max_TPM_somatic < 0.5) %>%
  dplyr::filter(ratio_testis_somatic >= 10) %>%
  pull(external_gene_name)

## Tag testis_preferential genes
## Compare to testis-specific genes, testis-preferential genes can have
## an expression level in a somatic tissue > 0.5 but this should
## only occur in a minority of somatic tissues (filter en q75_TPM_somatic)
testis_preferential <- GTEX_data %>%
  dplyr::filter(!external_gene_name %in% lowly_expressed) %>%
  dplyr::filter(!external_gene_name %in% testis_specific)  %>%
  dplyr::filter(Testis >= 1) %>%
  dplyr::filter(ratio_testis_somatic >= 10) %>%
  dplyr::filter(q75_TPM_somatic < 0.5) %>%
  pull(external_gene_name)

## Add GTEX_category column summarizing the category
## assigned to each gene using GTEx database.
GTEX_data <- GTEX_data %>%
  mutate(GTEX_category = case_when(
    external_gene_name %in% lowly_expressed ~ "lowly_expressed",
    external_gene_name %in% testis_specific ~ "testis_specific",
    external_gene_name %in% testis_preferential ~ "testis_preferential"))
GTEX_data[is.na(GTEX_data$GTEX_category), "GTEX_category"] <- "other"

Gtex_mat <- as.matrix(GTEX_data %>%
                        dplyr::select(-c(ensembl_gene_id, external_gene_name,
                                         max_TPM_somatic, ratio_testis_somatic,
                                         GTEX_category)))
rownames(Gtex_mat) <- GTEX_data$ensembl_gene_id
rowdata <- as.data.frame(GTEX_data %>%
                           dplyr::select(c(ensembl_gene_id, external_gene_name,
                                           GTEX_category,
                                           max_TPM_somatic )))
rowdata <- column_to_rownames(rowdata, "ensembl_gene_id")
GTEX_data <- SummarizedExperiment(assays = list(TPM = Gtex_mat),
                                  rowData = rowdata)

save(GTEX_data, file = "../../eh_data/GTEX_data.rda",
     compress = "xz",
     compression_level = 9)
