## Code to prepare `TCGA_TPM` dataset goes here

library("TCGAbiolinks")
library("tidyverse")
library("SummarizedExperiment")
library("BiocFileCache")
library("org.Hs.eg.db")

load(file = "../../eh_data/GTEX_data.rda")

bfc <- BiocFileCache(cache = "../BiocFileCache",
                     ask = FALSE)

# Load and save TCGA expression data
for(tumor_code in c("SKCM", "LUAD", "LUSC", "COAD", "ESCA",
                    "BRCA", "HNSC")) {

  rname <- paste0("TCGA_", tumor_code, "_STAR")
  if(length(bfcquery(bfc, rname)$rid) == 0) {

    savepath <- bfcnew(bfc, rname, ext=".RData")

    # Query data
    query_gene_counts_harmonized <-
      GDCquery(project = paste0("TCGA-", tumor_code),
               data.category = "Transcriptome Profiling",
               data.type = "Gene Expression Quantification",
               workflow.type = "STAR - Counts",
               legacy = FALSE)

    # Download files
    GDCdownload(query_gene_counts_harmonized, method = "api",
                directory = paste0("./tmp/TCGA_", tumor_code, "_STAR"),
                files.per.chunk = 20)

    # Create a Summarized-Experiment object
    GDCprepare(query_gene_counts_harmonized,
               directory = paste0("./tmp/TCGA_", tumor_code, "_STAR"),
               save = TRUE,
               save.filename = savepath,
               remove.files.prepared = FALSE)
  }
}

prepare_data <- function(tum) {
  rname <- paste0("TCGA_", tum, "_STAR")
  load(bfc[[bfcquery(bfc, rname)$rid]])

  assays(data) <- assays(data)["tpm_unstrand"]
  names(assays(data)) <- "TPM"

  data <- data[, colData(data)$shortLetterCode == "NT" |
                 colData(data)$shortLetterCode == "TP" |
                 colData(data)$shortLetterCode == "TM"]

  ## Remove duplicated ensembl ref associated to PAR_Y, and duplicated genes
  data <- data[grep("_PAR_", x = rownames(data), invert = TRUE, value = TRUE), ]
  data <- data[rowData(data)$gene_id[!duplicated(rowData(data)$gene_name)], ]
  ## Keep only genes present in Gtex data
  rownames(data) <- substr(rownames(data), 1, 15)
  data <- data[rownames(data) %in% rownames(GTEX_data), ]
  rowData(data) <- rowData(GTEX_data)[rownames(data), "external_gene_name"]
  names(rowData(data)) <- "external_gene_name"

  # Add frequencies of activation (TPM >= TPM_thr) of each gene in all tumor
  # samples.
  tumors_only <- data[ , colData(data)$shortLetterCode != 'NT']
  TPM_thr <- 1
  rowdata <- as_tibble(rowData(tumors_only), rownames = "ensembl_gene_id")
  binary <- ifelse(assay(tumors_only) >= TPM_thr, 1, 0)
  tmp <- rowSums(binary) / ncol(binary) * 100
  tmp <- enframe(tmp, name = "ensembl_gene_id",
                 value = paste0("percent_pos_", tum))
  rowdata <- left_join(rowdata, tmp)

  # Estimate the percent of tumors in which genes are repressed
  # (TPM < TPM_low_thr)
  TPM_low_thr <- 0.5
  binary <- ifelse(assay(tumors_only) <= TPM_low_thr, 1, 0)
  tmp <- rowSums(binary) / ncol(binary) * 100
  tmp <- enframe(tmp, name = "ensembl_gene_id",
                 value = paste0("percent_neg_", tum))
  rowdata <- rowdata %>%
    filter(ensembl_gene_id %in% rownames(GTEX_data)) %>%
    left_join(tmp)

  # Estimate the q75 expression in normal peritumoral tissues
  # (only when at least 20 NT samples =>
  # exclude SKCM and ESCA where not enough NT samples)
  # Will be used to filter out genes that are expressed in a high proportion of
  # peritumoral tissues => not strictly testis-specific
  NT_only <- data[ , colData(data)$shortLetterCode == 'NT']
  if (ncol(NT_only) >= 20) {
    q75_in_NT <-  tibble(ensembl_gene_id = rownames(NT_only),
                         q75 = rowQuantiles(assay(NT_only), probs = 0.75))
    names(q75_in_NT)[-1] <- paste0('TPM_q75_in_NT_', tum)
    rowdata <- rowdata %>%
      left_join(q75_in_NT)
  }

  # Max expression (TPM) in a tumor
  max_TPM <- tibble(ensembl_gene_id = rownames(tumors_only),
                    max_TPM = rowMax(assay(tumors_only)))
  names(max_TPM) <- c("ensembl_gene_id", paste0("max_TPM_in_", tum))

  rowData(data) <- rowdata %>%
    left_join(max_TPM)
  return(assign(x = paste0(tum, "_TPM"), value = data))
}

SKCM <- prepare_data(tum = "SKCM")
LUAD <- prepare_data(tum = "LUAD")
LUSC <- prepare_data(tum = "LUSC")
COAD <- prepare_data(tum = "COAD")
ESCA <- prepare_data(tum = "ESCA")
BRCA <- prepare_data(tum = "BRCA")
HNSC <- prepare_data(tum = "HNSC")

# # Sanity check
# all(rownames(SKCM) == rownames(LUAD),
#     rownames(SKCM) == rownames(LUSC),
#     rownames(SKCM) == rownames(COAD),
#     rownames(SKCM) == rownames(BRCA),
#     rownames(SKCM) == rownames(ESCA),
#     rownames(SKCM) == rownames(HNSC))

coldata_common_variables <-
  colnames(colData(SKCM))[
    colnames(colData(SKCM)) %in% colnames(colData(LUAD)) &
      colnames(colData(SKCM)) %in% colnames(colData(LUSC)) &
      colnames(colData(SKCM)) %in% colnames(colData(COAD)) &
      colnames(colData(SKCM)) %in% colnames(colData(ESCA)) &
      colnames(colData(SKCM)) %in% colnames(colData(BRCA)) &
      colnames(colData(SKCM)) %in% colnames(colData(HNSC))]

TPM <- cbind(assay(SKCM), assay(LUAD), assay(LUSC), assay(COAD),
             assay(ESCA), assay(BRCA), assay(HNSC))

rowdata <- as_tibble(rowData(SKCM)) %>%
  left_join(as_tibble(rowData(LUAD))) %>%
  left_join(as_tibble(rowData(LUSC))) %>%
  left_join(as_tibble(rowData(COAD))) %>%
  left_join(as_tibble(rowData(ESCA))) %>%
  left_join(as_tibble(rowData(BRCA))) %>%
  left_join(as_tibble(rowData(HNSC)))

coldata <- rbind(colData(SKCM)[, coldata_common_variables],
                 colData(LUAD)[, coldata_common_variables],
                 colData(LUSC)[, coldata_common_variables],
                 colData(COAD)[, coldata_common_variables],
                 colData(ESCA)[, coldata_common_variables],
                 colData(BRCA)[, coldata_common_variables],
                 colData(HNSC)[, coldata_common_variables])

## Add to colData `global hypomethylation levels` from paper:
## DNA methylation loss promotes immune evasion of tumours with high
## mutation and copy number load. Jang et al., Nature Commun 2019
## Keep also `CD8 T cells` and `Proliferation score` columns
global_hypo <- readxl::read_xlsx(
  "../extdata/41467_2019_12159_MOESM4_ESM.xlsx", skip = 3)
names(global_hypo) <- c("project_id", "Sample", "global_methylation",
                        "CD8_T_cells", "proliferation_score")
global_hypo$project_id <- paste0("TCGA-", global_hypo$project_id)
coldata$Sample <- substr(coldata$sample, 1, 15)
coldata <- as_tibble(coldata) %>%
  left_join(global_hypo) %>%
  as.data.frame()
coldata <- column_to_rownames(coldata, "barcode")
TCGA_TPM <- SummarizedExperiment(assays = list(TPM = TPM),
                                 colData = coldata,
                                 rowData = rowdata)

# Add frequencies of activation (TPM >= TPM_thr) of each gene in all types of
# tumor samples.
tumors_only <- TCGA_TPM[, colData(TCGA_TPM)$shortLetterCode != 'NT']
TPM_thr <- 1
binary <- ifelse(assay(tumors_only) >= TPM_thr, 1, 0)
tmp <- rowSums(binary) / ncol(binary) * 100
tmp <- enframe(tmp, name = "ensembl_gene_id", value = "percent_pos_tum")
rowdata <- rowdata %>%
  left_join(tmp)

# Estimate the percent of tumors in which genes are repressed
# (TPM < TPM_low_thr)
TPM_low_thr <- 0.5
binary <- ifelse(assay(tumors_only) <= TPM_low_thr, 1, 0)
tmp <- rowSums(binary) / ncol(binary) * 100
tmp <- enframe(tmp, name = "ensembl_gene_id", value = "percent_neg_tum")
rowdata <- left_join(rowdata, tmp)

rowdata$max_TPM_in_TCGA <- rowMax(assay(tumors_only))
rowdata$max_q75_in_NT <- rowMax(as.matrix(
  rowdata %>% dplyr::select(starts_with("TPM_q75"))))

rowdata <- rowdata %>%
  mutate(TCGA_category = case_when(
    percent_neg_tum < 20 ~ "leaky",
    percent_neg_tum >= 20 & percent_pos_tum >= 1  &
      max_TPM_in_TCGA >= 5 ~ "activated",
    percent_neg_tum >= 20 & percent_pos_tum < 1 ~ "not_activated",
    percent_neg_tum >= 20 & percent_pos_tum >= 1 &
      max_TPM_in_TCGA < 5 ~ "lowly_activated"))


rowdata <- as.data.frame(rowdata)
rowdata <- column_to_rownames(rowdata, "ensembl_gene_id")

rowData(TCGA_TPM) <- rowdata

save(TCGA_TPM, file = "../../eh_data/TCGA_TPM.rda",
     compress = "xz",
     compression_level = 9)
