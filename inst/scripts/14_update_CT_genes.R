## Code to prepare `CT_genes` dataset goes here

library("tidyverse")
library("SingleCellExperiment")

load("../../eh_data/CT_genes.rda")
load("../../eh_data/testis_sce.rda")

################################################################################
## Add `testis_cell_type` column, specifying the testis cell-type showing the
## highest median expression of each gene.
################################################################################
CT_genes <- CT_genes %>%
  left_join(as_tibble(rowData(testis_sce), rownames = "external_gene_name"))
table(CT_genes$testis_cell_type, useNA = 'ifany')

save(CT_genes, file = "../../eh_data/CT_genes.rda",
     compress = "xz",
     compression_level = 9)



