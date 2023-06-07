## Code to prepare `testis_sce` dataset goes here


library(tidyverse)
library(SingleCellExperiment)
library(scater)

## Data from The adult human testis transcriptional cell atlas (Guo et al. 2018)
## `GSE112013_Combined_UMI_table.txt.gz` was downloaded from GEO (accession:
## GSE11201). `TableS1` (cell metadata) comes from the paper's supplemental data.

## ! `TableS1` could only be downoaded as a pdf file!
## => Copied the pdf content and pasted it in data/Guo_et_al_TableS1 text file.
command <- "cat ../extdata/Guo_et_al_TableS1 | tr ' ' '\n' > ../extdata/Guo_et_al_TableS1_tmp"
system(command)

tableS1 <- read_table(file = "../extdata/Guo_et_al_TableS1_tmp",
                      col_names = FALSE, skip = 5)

CellID <- vector()
nGene <- vector()
nUMI <- vector()
Final_clusters <- vector()

for (i in seq(1, nrow(tableS1), by = 4)) {
  CellID <- c(CellID, tableS1$X1[i])
  nGene <- c(nGene, tableS1$X1[i+1])
  nUMI <- c(nUMI, tableS1$X1[i+2])
  Final_clusters <- c(Final_clusters, tableS1$X1[i+3])
}

metadata <- tibble(CellID = CellID,
                   nGene = nGene,
                   nUMI = nUMI,
                   clusters = Final_clusters) %>%
  mutate(type = case_when(clusters == 1 ~ "SSC",
                          clusters == 2 ~ "Spermatogonia",
                          clusters == 3 ~ "Early_spermatocyte",
                          clusters == 4 ~ "Late_spermatocyte",
                          clusters == 5 ~ "Round_spermatid",
                          clusters == 6 ~ "Elongated_spermatid",
                          clusters == 7 ~ "Sperm1",
                          clusters == 8 ~ "Sperm2",
                          clusters == 9 ~ "Macrophage",
                          clusters == 10 ~ "Endothelial",
                          clusters == 11 ~ "Myoid",
                          clusters == 12 ~ "Sertoli",
                          clusters == 13 ~ "Leydig")) %>%
  mutate(Donor = sub(pattern = "-\\w*-\\d$", x = CellID, replacement = ''))

metadata$type <- factor(metadata$type,
                        levels = c("SSC", "Spermatogonia", "Early_spermatocyte",
                                   "Late_spermatocyte", "Round_spermatid",
                                   "Elongated_spermatid", "Sperm1", "Sperm2",
                                   "Macrophage", "Endothelial", "Myoid",
                                   "Sertoli", "Leydig"))

metadata$clusters <- as.factor(metadata$clusters)

counts <- read_tsv(file = "../extdata/GSE112013_Combined_UMI_table.txt")

mat <- as.matrix(counts[, -1])
rownames(mat) <- counts$Gene

coldata <- data.frame(metadata[, -1], row.names = metadata$CellID)

testis_sce <- SingleCellExperiment(assays = list(counts = mat[, rownames(coldata)]),
                                   colData = coldata)

# keep <- rowSums(assay(testis_sce)) > 0
# testis_sce <- testis_sce[keep,]

testis_sce <- logNormCounts(testis_sce)

##########################################################################
## Add a `testis_cell_type` column in the rowData, specifying the testis
## cell-type showing the highest mean expression of each gene.
##########################################################################

testis_cell_type <- as_tibble(logcounts(testis_sce),
                              rownames = "external_gene_name") %>%
  pivot_longer(names_to = "CellID", values_to = "logCounts", -external_gene_name) %>%
  left_join(as_tibble(colData(testis_sce), rownames = "CellID") %>%
              select(CellID, type)) %>%
  group_by(external_gene_name, type) %>%
  summarise(mean_exp = mean(logCounts)) %>%
  filter(mean_exp > 0) %>%
  filter(mean_exp == max(mean_exp)) %>%
  dplyr::rename(testis_cell_type = type) %>%
  select(external_gene_name, testis_cell_type) %>%
  unique()

rowData(testis_sce) <- tibble(external_gene_name = rownames(testis_sce)) %>%
  left_join(testis_cell_type)

save(testis_sce, file = "../../eh_data/testis_sce.rda",
     compress = "xz",
     compression_level = 9)
