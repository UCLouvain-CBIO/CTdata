## Code to prepare `CT_genes` dataset goes here

library("tidyverse")
library("SingleCellExperiment")

load("../../eh_data/CT_genes.rda")
load("../../eh_data/scRNAseq_HPA.rda")

##########################################################################
## Add the `Higher_in_somatic_cell_type` column (based on scRNAseq data
## of normal tissues from the Human Protei Atlas),
## specifying if some somatic cell types express the genes at a higher
## level than the level found in germ cell types, and filtered the CT_genes
## table accordingly.
##########################################################################
CT_genes <- CT_genes %>%
  left_join(as_tibble(rowData(scRNAseq_HPA)) %>%
              dplyr::select(external_gene_name, Higher_in_somatic_cell_type)) %>%
  filter(is.na(Higher_in_somatic_cell_type) | Higher_in_somatic_cell_type == FALSE)

save(CT_genes, file = "../../eh_data/CT_genes.rda",
     compress = "xz",
     compression_level = 9)



