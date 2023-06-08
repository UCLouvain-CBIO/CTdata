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

testis_sce <- logNormCounts(testis_sce)

##########################################################################
## Add 3 columns to the rowData:
## - `percent_pos_testis_germcells` giving the percent of testis germ cells
## in which the genes are detected (count > 0)
## - `percent_pos_testis_somatic` giving the percent of testis somatic cells
## in which the genes are detected (count > 0)
## - `testis_cell_type` column, specifying the testis cell-type showing
## the highest mean expression of each gene.
## Genes must be detected (count > 0) in at least 1% of the cells within
## a specific testis cell-type to be assigned to it.
##########################################################################

germ_cells <- c("SSC", "Spermatogonia", "Early_spermatocyte",
                "Late_spermatocyte","Round_spermatid", "Elongated_spermatid",
                "Sperm1", "Sperm2")
somatic_cells <- c("Macrophage", "Endothelial", "Myoid", "Sertoli", "Leydig")
n_germ_cells <- dim(testis_sce[, testis_sce$type %in% germ_cells])[2]
n_somatic_cells <- dim(testis_sce[, testis_sce$type %in% somatic_cells])[2]

percent_pos_germcells <- as_tibble(logcounts(testis_sce),
                                   rownames = "external_gene_name") %>%
  pivot_longer(names_to = "CellID", values_to = "logCounts", -external_gene_name) %>%
  left_join(as_tibble(colData(testis_sce), rownames = "CellID") %>%
              select(CellID, type)) %>%
  filter(type %in% germ_cells) %>%
  group_by(external_gene_name) %>%
  summarise(n_pos = count(logCounts > 0)) %>%
  mutate(n_germ_cells = n_germ_cells, percent_pos_testis_germcells = n_pos / n_germ_cells * 100)

percent_pos_somatic <- as_tibble(logcounts(testis_sce),
                                 rownames = "external_gene_name") %>%
  pivot_longer(names_to = "CellID", values_to = "logCounts", -external_gene_name) %>%
  left_join(as_tibble(colData(testis_sce), rownames = "CellID") %>%
              select(CellID, type)) %>%
  filter(type %in% somatic_cells) %>%
  group_by(external_gene_name) %>%
  summarise(n_pos = count(logCounts > 0)) %>%
  mutate(n_germ_cells = n_somatic_cells, percent_pos_testis_somatic = n_pos / n_germ_cells * 100)

testis_cell_type <- as_tibble(logcounts(testis_sce),
                              rownames = "external_gene_name") %>%
  pivot_longer(names_to = "CellID", values_to = "logCounts", -external_gene_name) %>%
  left_join(as_tibble(colData(testis_sce), rownames = "CellID") %>%
              select(CellID, type)) %>%
  group_by(external_gene_name, type) %>%
  summarise(mean_exp = mean(logCounts), n_pos = count(logCounts > 0)) %>%
  left_join(enframe(table(testis_sce$type), name = "type", value = "n_cells")) %>%
  mutate(n_cells = as.vector(n_cells), percent_pos = n_pos / n_cells * 100) %>%
  filter(percent_pos > 1) %>%
  filter(mean_exp == max(mean_exp)) %>%
  dplyr::rename(testis_cell_type = type) %>%
  select(external_gene_name, testis_cell_type) %>%
  unique()

rowData(testis_sce) <- tibble(external_gene_name = rownames(testis_sce)) %>%
  left_join(percent_pos_germcells %>% select(external_gene_name, percent_pos_testis_germcells)) %>%
  left_join(percent_pos_somatic %>% select(external_gene_name, percent_pos_testis_somatic)) %>%
  left_join(testis_cell_type)

save(testis_sce, file = "../../eh_data/testis_sce.rda",
     compress = "xz",
     compression_level = 9)
