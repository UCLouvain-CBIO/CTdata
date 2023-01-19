## Code to prepare `CCLE_data` dataset goes here

library(tidyverse)
library(SummarizedExperiment)
library(BiocFileCache)
library(biomaRt)

bfc <- BiocFileCache(cache = "../BiocFileCache",
                     ask = FALSE)

if (length(bfcquery(bfc, "CCLE_data")$rid) == 0) {
  url <- "https://ndownloader.figshare.com/files/34989922"
  bfcadd(bfc, "CCLE_data", fpath = url)
}

if (length(bfcquery(bfc, "CCLE_cell_metadata")$rid) == 0) {
  url <- "https://ndownloader.figshare.com/files/35020903"
  bfcadd(bfc, "CCLE_cell_metadata", fpath = url)
}

CCLE_full_data <- read_csv(bfc[[bfcquery(bfc, "CCLE_data")$rid]])
names(CCLE_full_data)[1] <- "DepMap_ID"
CCLE_metadata <- read_csv(bfc[[bfcquery(bfc, "CCLE_cell_metadata")$rid]])

## Select tumor cell lines and clean "primary_disease" category
coldata <- CCLE_metadata %>%
  filter(DepMap_ID %in% CCLE_full_data$DepMap_ID) %>%
  filter(primary_disease != "Unknown") %>%
  filter(primary_disease != "Fibroblast") %>%
  filter(primary_disease != "Engineered") %>%
  filter(primary_disease != "Non-Cancerous") %>%
  mutate(type = sub(pattern = " Cancer", x = primary_disease, replacement = '')) %>%
  mutate(type = sub(pattern = "Colon/", x = type, replacement = '')) %>%
  mutate(type = sub(pattern = "Endometrial/", x = type, replacement = '')) %>%
  mutate(type = gsub(pattern = " ", x = type, replacement = '_')) %>%
  filter(!is.na(stripped_cell_line_name)) %>%
  as.data.frame()

## Keep cancer types with more than 30 cell lines
tumor_types <- coldata %>%
  dplyr::group_by(type) %>%
  dplyr::summarize(n = n()) %>%
  filter(n >= 30) %>%
  pull(type)

coldata <- coldata[coldata$type %in% tumor_types, ]
rownames(coldata) <- coldata$stripped_cell_line_name

## Reformat CCLE data keeping only genes present in GTEx data
TPM <- CCLE_full_data %>%
  pivot_longer(names_to = "gene", values_to = "expression", -DepMap_ID) %>%
  mutate(ensembl_gene_id = str_extract(gene, "ENSG\\d{11}")) %>%
  dplyr::select(-gene) %>%
  filter(DepMap_ID %in% coldata$DepMap_ID) %>%
  filter(ensembl_gene_id %in% rownames(GTEX_data)) %>%
  left_join(as_tibble(rowData(GTEX_data), rownames = "ensembl_gene_id") %>%
              dplyr::select(ensembl_gene_id, external_gene_name)) %>%
  left_join(coldata %>%
              dplyr::select(DepMap_ID, stripped_cell_line_name)) %>%
  dplyr::select(-DepMap_ID) %>%
  pivot_wider(names_from = stripped_cell_line_name, values_from = expression)

TPM_mat <- TPM %>%
  dplyr::select(-ensembl_gene_id, -external_gene_name) %>%
  as.matrix()
rownames(TPM_mat) <- TPM$ensembl_gene_id

## Values are in log2(TPM + 1), convert to TPM
TPM_mat <- 2^(TPM_mat) - 1

# Estimate frequencies of activation of each CT gene
# in all selected cell lines.
# Gene are considered as "activated" if TPM >= TPM_thr
TPM_thr <- 10
activation_frequencies <- tibble(ensembl_gene_id = rownames(TPM_mat)) %>%
  left_join(as_tibble(rowData(GTEX_data), rownames = "ensembl_gene_id") %>%
              dplyr::select(ensembl_gene_id, external_gene_name))
binary <- ifelse(TPM_mat >= TPM_thr, 1, 0)
tmp <- rowSums(binary) / ncol(binary) * 100
tmp <- enframe(tmp, name = "ensembl_gene_id",
               value = "percent_of_positive_CCLE_cell_lines")
activation_frequencies <- left_join(activation_frequencies, tmp)

# Estimate the percent of cell line in which the genes
# are not expressed (TPM <= TPM_low_thr)
TPM_low_thr <- 0.1
binary <- ifelse(TPM_mat <= TPM_low_thr, 1, 0)
tmp <- rowSums(binary) / ncol(binary) * 100
repression_frequencies <- enframe(tmp, name = "ensembl_gene_id",
                                  value = "percent_of_negative_CCLE_cell_lines")

# Max expression (TPM) in a cell line
max_TPM <- tibble(ensembl_gene_id = rownames(TPM_mat),
                  max_TPM_in_CCLE = rowMax(TPM_mat))

rowdata <- left_join(activation_frequencies, repression_frequencies) %>%
  left_join(max_TPM) %>%
  mutate(percent_of_positive_CCLE_cell_lines =
           round(percent_of_positive_CCLE_cell_lines, 1)) %>%
  mutate(percent_of_negative_CCLE_cell_lines =
           round(percent_of_negative_CCLE_cell_lines, 1)) %>%
  mutate(CCLE_category =
           case_when(is.na(percent_of_positive_CCLE_cell_lines) ~ "not_in_CCLE",
                     percent_of_negative_CCLE_cell_lines < 20 ~ "leaky",
                     percent_of_negative_CCLE_cell_lines >= 20 &
                       percent_of_positive_CCLE_cell_lines > 0 ~ "activated",
                     percent_of_negative_CCLE_cell_lines >= 20 &
                       percent_of_positive_CCLE_cell_lines == 0 ~ "not_activated")) %>%
  column_to_rownames("ensembl_gene_id")

CCLE_data <- SummarizedExperiment(assays = list(TPM = TPM_mat),
                                  rowData = rowdata,
                                  colData = coldata[colnames(TPM_mat), ])

save(CCLE_data, file = "../../eh_data/CCLE_data.rda",
     compress = "xz",
     compression_level = 9)


