## Code to prepare `HPA_cell_type_specificity` dataset goes here

library("tidyverse")
library("SummarizedExperiment")
load("../../eh_data/GTEX_data.rda")

# Data from the Human Protein Atlas (version 23.0) in tab-separated format

proteinatlas <- read_tsv("../extdata/proteinatlas.tsv")


# Use the column `RNA single cell type specific nTPM` to flag gene
# as "not_testis_specific" genes if specific of a cell type other than a germ
# cell or a placental cell and as "testis_specific" genes if specific of a germ
# cell type or a placental cell

cell_type_specificities_list <- strsplit(
  proteinatlas$`RNA single cell type specific nTPM`, split = ";")
names(cell_type_specificities_list) <- proteinatlas$Gene
detected_specificities <- unique(gsub(
  pattern = ": .*", x = unlist(strsplit(
    proteinatlas$`RNA single cell type specific nTPM`, split = ";")), ''))

accepted_cell_types <- c("Spermatocytes", "Spermatogonia", "Early spermatids",
                         "Late spermatids", "Oocytes", "Syncytiotrophoblasts",
                         "Cytotrophoblasts", "Extravillous trophoblasts", NA)

somatic_cell_types <- detected_specificities[!detected_specificities %in%
                                               accepted_cell_types]

suspect <- names(grep(somatic_cell_types[1], cell_type_specificities_list,
                      value = TRUE))
for (tissue in somatic_cell_types[-1]){
  tmp <- names(grep(tissue, cell_type_specificities_list, value = TRUE))
  suspect <- c(suspect, tmp)
}
suspect_HPA <- unique(suspect)

HPA_cell_type_specificities <- tibble(
  ensembl_gene_id = rownames(GTEX_data),
  external_gene_name = rowData(GTEX_data)$external_gene_name) %>%
  mutate(HPA_category = case_when(
    !external_gene_name %in% proteinatlas$Gene ~ NA,
    external_gene_name %in% suspect_HPA ~ "not_testis_specific",
    external_gene_name %in% proteinatlas$Gene &
      !external_gene_name %in% suspect_HPA ~ "testis_specific")) %>%
  left_join(
    proteinatlas %>%
      dplyr::rename(ensembl_gene_id = Ensembl) %>%
      dplyr::select(ensembl_gene_id, `RNA single cell type specific nTPM`)) %>%
  dplyr::rename(HPA_RNA_single_cell_type_specific_nTPM =
                  `RNA single cell type specific nTPM`)

save(HPA_cell_type_specificities,
     file = "../../eh_data/HPA_cell_type_specificities.rda",
     compress = "xz",
     compression_level = 9)
