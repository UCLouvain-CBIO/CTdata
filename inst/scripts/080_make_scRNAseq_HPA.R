## Code to prepare `scRNAseq_HPA` dataset goes here

library("tidyverse")
library("SummarizedExperiment")
library("SingleCellExperiment")

## scRNAseq from normal tissues downloaded from the Human Protein Atlas:
## - rna_single_cell_type.tsv downloaded directly from
## https://www.proteinatlas.org/about/download
## - The cell types group associated to each cell type were copied manually from
## https://www.proteinatlas.org/humanproteome/single+cell+type

data <- read_tsv("../extdata/rna_single_cell_type.tsv")
data <- data %>%
  dplyr::rename(ensembl_gene_id = Gene, external_gene_name = "Gene name") %>%
  pivot_wider(names_from = "Cell type", values_from = "nTPM")
mat <- as.matrix(data[, -(1:2)])
rownames(mat) <- data$ensembl_gene_id
colnames(mat) <- str_to_sentence(colnames(mat))

rowdata <- data.frame("external_gene_name" = data$external_gene_name)
rownames(rowdata) <- data$ensembl_gene_id

cell_types <- c("Adipocyte_cells", "Blood_immune_cells", "Endocrine_cells",
                "Endothelial_cells", "Germ_cells", "Glandular_cells",
                "Glial_cells", "Mesenchymal_cells", "Pigment_cells",
                "Specialized_epithelial_cells", "Trophoblast_cells")

cell_types_by_group <- read_tsv(paste0("../extdata/", cell_types[1]), col_select = 1,
                                col_names = c("Cell_type")) %>%
  mutate(group = cell_types[1])

for (cell in cell_types[-1]){
  tmp <- read_tsv(paste0("../extdata/", cell), col_select = 1,
                  col_names = c("Cell_type")) %>%
    mutate(group = cell)
  cell_types_by_group <- rbind(cell_types_by_group, tmp)
}
cell_types_by_group$Cell_type <- str_to_sentence(cell_types_by_group$Cell_type)

# Rm "Hofbauer cells" (placental cells => confusing)
cell_types_by_group <- cell_types_by_group %>%
  filter(Cell_type != "Hofbauer cells")

scRNAseq_HPA <- SingleCellExperiment(assays = list(TPM = mat[, cell_types_by_group$Cell_type]),
                                     colData = cell_types_by_group,
                                     rowData = rowdata)

not_somatic_group <- c("Germ_cells", "Trophoblast_cells")
somatic_groups <- unique(cell_types_by_group$group[
  !cell_types_by_group$group %in% not_somatic_group])

##########################################################################
## Add 3 columns to the rowData:
## - `values_in_somatics` giving the maximum expression value found in a
## somatic cell type
## - `max_in_germcells_group` giving the maximum expression value found in
## a germ cell type
## - `Higher_in_somatic_cell_type` column, specifying if a somatic cell type
## was found to express the gene at a higher level than any germ cell type
##########################################################################

values_in_somatics <- as_tibble(assay(scRNAseq_HPA),
                                rownames = "ensembl_gene_id") %>%
  pivot_longer(names_to = "Cell_type", values_to = "counts",
               -ensembl_gene_id) %>%
  left_join(as_tibble(colData(scRNAseq_HPA))) %>%
  filter(group %in% somatic_groups) %>%
  group_by(ensembl_gene_id, group) %>%
  summarise(max_TPM_in_a_somatic_cell_type = max(counts)) %>%
  filter(max_TPM_in_a_somatic_cell_type == max(max_TPM_in_a_somatic_cell_type)) %>%
  dplyr::select(-group) %>%
  unique()

values_in_germCells <- as_tibble(assay(scRNAseq_HPA),
                                 rownames = "ensembl_gene_id") %>%
  pivot_longer(names_to = "Cell_type", values_to = "counts",
               -ensembl_gene_id) %>%
  left_join(as_tibble(colData(scRNAseq_HPA))) %>%
  filter(group == "Germ_cells") %>%
  group_by(ensembl_gene_id) %>%
  summarise(max_in_germcells_group = max(counts)) %>%
  unique()

rowdata <- as_tibble(rowData(scRNAseq_HPA), rownames = "ensembl_gene_id") %>%
  left_join(values_in_somatics) %>%
  left_join(values_in_germCells) %>%
  mutate(Higher_in_somatic_cell_type =
           case_when((max_TPM_in_a_somatic_cell_type > 0 |
                        max_in_germcells_group > 0) &
                       max_TPM_in_a_somatic_cell_type > max_in_germcells_group
                     ~ TRUE,
                     (max_TPM_in_a_somatic_cell_type > 0 |
                        max_in_germcells_group > 0) &
                       max_TPM_in_a_somatic_cell_type <= max_in_germcells_group
                     ~ FALSE,
                     (max_TPM_in_a_somatic_cell_type == 0 &
                        max_in_germcells_group == 0) ~ NA)) %>%
  as.data.frame()

rownames(rowdata) <- rowdata$ensembl_gene_id
rowData(scRNAseq_HPA) <- rowdata[, -1]

save(scRNAseq_HPA, file = "../../eh_data/scRNAseq_HPA.rda",
     compress = "xz",
     compression_level = 9)
