## Code to prepare `HPA_cell_type_specificity` dataset goes here

library("tidyverse")
library("SummarizedExperiment")
load("../../eh_data/GTEX_data.rda")

# Data from the Human Protein Atlas (version 23.0) in tab-separated format

proteinatlas <- read_tsv("../extdata/proteinatlas.tsv")


# Use the column `RNA single cell type specific nTPM` to flag gene
# as "not_testis_specific" genes when they are specific of a somatic cell type
HPA_single_cell_specificities <-
  proteinatlas$`RNA single cell type specific nTPM`[
    !is.na(proteinatlas$`RNA single cell type specific nTPM`)]

HPA_single_cell_specificities_list <- strsplit(
  HPA_single_cell_specificities, split = ";")
names(HPA_single_cell_specificities_list) <- proteinatlas[
  !is.na(proteinatlas$`RNA single cell type specific nTPM`),]$Gene

detected_specificities <- unique(gsub(
  pattern = ": .*", x = unlist(strsplit(
    HPA_single_cell_specificities, split = ";")), ''))

germ_cell_types <- c("Spermatocytes", "Spermatogonia", "Early spermatids",
                     "Late spermatids", "Oocytes")

accepted_cell_types <- c("Spermatocytes", "Spermatogonia", "Early spermatids",
                         "Late spermatids", "Oocytes", "Syncytiotrophoblasts",
                         "Cytotrophoblasts", "Extravillous trophoblasts")

somatic_cell_types <- detected_specificities[!detected_specificities %in%
                                               accepted_cell_types]


x <- unlist(strsplit(HPA_single_cell_specificities_list[[1]], split = ";"))
cells <- gsub(pattern = ":.*$", x = x, replacement = '')
values <- gsub(pattern = "^.*: ", x = x, replacement = '')
max_som <- tibble(cell_type = cells,
                  nTPM = as.numeric(values)) %>%
  arrange(desc(nTPM)) %>%
  filter(cell_type %in% somatic_cell_types) %>%
  head(1) %>%
  pull(nTPM)
max_germ <- tibble(cell_type = cells,
                   nTPM = as.numeric(values)) %>%
  arrange(desc(nTPM)) %>%
  filter(cell_type %in% germ_cell_types) %>%
  head(1) %>%
  pull(nTPM)
HPA_nTPM <- tibble(external_gene_name = names(HPA_single_cell_specificities_list[1]),
       max_HPA_somatic = ifelse(length(max_som) > 0, max_som, 0),
       max_HPA_germcell = ifelse(length(max_germ) > 0, max_germ, 0))

for (i in seq_along(HPA_single_cell_specificities_list)[-1]) {
  x <- unlist(strsplit(HPA_single_cell_specificities_list[[i]], split = ";"))
  cells <- gsub(pattern = ":.*$", x = x, replacement = '')
  values <- gsub(pattern = "^.*: ", x = x, replacement = '')
  max_som <- tibble(cell_type = cells,
                    nTPM = as.numeric(values)) %>%
    arrange(desc(nTPM)) %>%
    filter(cell_type %in% somatic_cell_types) %>%
    head(1) %>%
    pull(nTPM)
  max_germ <- tibble(cell_type = cells,
                     nTPM = as.numeric(values)) %>%
    arrange(desc(nTPM)) %>%
    filter(cell_type %in% germ_cell_types) %>%
    head(1) %>%
    pull(nTPM)
  tmp <- tibble(external_gene_name = names(HPA_single_cell_specificities_list[i]),
                max_HPA_somatic = ifelse(length(max_som) > 0, max_som, 0),
                max_HPA_germcell = ifelse(length(max_germ) > 0, max_germ, 0))
  HPA_nTPM <- rbind(HPA_nTPM, tmp)
}

HPA_nTPM <- HPA_nTPM[!duplicated(HPA_nTPM$external_gene_name),]

HPA_cell_type_specificities <- tibble(
  ensembl_gene_id = rownames(GTEX_data),
  external_gene_name = rowData(GTEX_data)$external_gene_name) %>%
  left_join(
    proteinatlas %>%
      dplyr::rename(ensembl_gene_id = Ensembl) %>%
      dplyr::rename(HPA_scRNAseq_celltype_specific_nTPM = `RNA single cell type specific nTPM`) %>%
      dplyr::select(ensembl_gene_id, HPA_scRNAseq_celltype_specific_nTPM)) %>%
  left_join(HPA_nTPM) %>%
  mutate(not_detected_in_somatic_HPA = case_when(
    max_HPA_somatic == 0 ~ TRUE,
    max_HPA_somatic > 0 ~ FALSE)) %>%
  mutate(HPA_ratio_germ_som = max_HPA_germcell / max_HPA_somatic)

save(HPA_cell_type_specificities,
     file = "../../eh_data/HPA_cell_type_specificities.rda",
     compress = "xz",
     compression_level = 9)
