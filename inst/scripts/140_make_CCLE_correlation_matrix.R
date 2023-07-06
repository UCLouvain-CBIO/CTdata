## Code to prepare `CCLE_correlation_matrix` dataset goes here

library("tidyverse")
load("../../eh_data/CCLE_data.rda")
load("../../eh_data/CT_genes.rda")

x <- as_tibble(assay(CCLE_data), rownames = "ensembl_gene_id") %>%
  pivot_longer(names_to = "cell_line", values_to = "TPM", -ensembl_gene_id) %>%
  mutate(TPM = log1p(TPM)) %>%
  pivot_wider(names_from = ensembl_gene_id, values_from = TPM) %>%
  dplyr::select(-cell_line) %>%
  as.matrix()

CCLE_correlation_matrix <- cor(x, method = "pearson")

CCLE_correlation_matrix <- CCLE_correlation_matrix[CT_genes$ensembl_gene_id, ]

save(CCLE_correlation_matrix, file = "../../eh_data/CCLE_correlation_matrix.rda",
     compress = "xz",
     compression_level = 9)
