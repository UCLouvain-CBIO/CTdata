## Code to prepare `CT_genes` dataset goes here

library("tidyverse")
library("SingleCellExperiment")

load("../../eh_data/CT_genes.rda")
load("../../eh_data/testis_sce.rda")

##########################################################################
## Add 3 columns based on testis scRNAseq data:
## - `percent_pos_testis_germcells` giving the percent of testis germ cells
## in which the genes are detected (count > 0)
## - `percent_pos_testis_somatic` giving the percent of testis somatic cells
## in which the genes are detected (count > 0)
## - `testis_cell_type` column, specifying the testis cell-type showing
## the highest mean expression of each gene.
##########################################################################

CT_genes <- CT_genes %>%
  left_join(as_tibble(rowData(testis_sce)))

save(CT_genes, file = "../../eh_data/CT_genes.rda",
     compress = "xz",
     compression_level = 9)



