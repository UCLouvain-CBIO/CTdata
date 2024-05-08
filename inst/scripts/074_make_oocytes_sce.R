## Code to prepare `oocytes_sce` dataset goes here

library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(biomaRt)

## Data from paper: Decoding dynamic epigenetic landscapes in human oocytes
## using single-cell multi-omics sequencing (Yan et al. Cell Stem Cell 2021)

## `GSE154762_hO_scChaRM_count_matix.txt.gz` was downloaded from GEO (accession:
## GSE154762).

table <- read_table(file = "../extdata/GSE154762_hO_scChaRM_count_matix.txt.gz",
                    skip = 1, col_names = F)
header <- read_table(file = "../extdata/GSE154762_hO_scChaRM_count_matix.txt.gz",
                     col_names = F,  n_max = 1)
colnames(table) <- c("gene", header)

metadata <- tibble(cell = colnames(table)[-1],
                   type = gsub("RNA_", "", cell)) %>%
  mutate(type = gsub("_sc.*$", "", type)) %>%
  mutate(type = case_when(type %in% c("GO1", "GO2") ~ "Growing oocytes",
                          type == "FGO" ~ "Fully grown oocytes",
                          type == "MI" ~ "Metaphase I",
                          type == "MII" ~ "Metaphase II",
                          !type %in% c("GO1", "GO2") ~ type)) %>%
  mutate(type = factor(type, levels =
                         c("Growing oocytes", "Fully grown oocytes",
                           "Metaphase I", "Metaphase II", "Granulosa",
                           "Immune", "StromaC1", "StromaC2"))) %>%
  mutate(sex = "F", stage = "meiotic") %>%
  mutate(germcell = ifelse(
    type %in% c("Growing oocytes", "Fully grown oocytes",
                "Metaphase I", "Metaphase II"), TRUE, FALSE))

# Some gene names are not the official ones!
ensembl <- biomaRt::useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
attributes_vector <- c("ensembl_gene_id",
                       "external_gene_name",
                       "external_synonym",
                       "ensembl_transcript_id",
                       "external_transcript_name",
                       "chromosome_name",
                       "transcript_biotype",
                       "transcript_is_canonical")
transcripts_infos <- as_tibble(biomaRt::getBM(attributes = attributes_vector,
                                              mart = ensembl))
canonical_transcripts <- transcripts_infos %>%
  dplyr::filter(transcript_is_canonical == 1) %>%
  dplyr::filter(external_transcript_name != "") %>%
  dplyr::filter(chromosome_name %in% c(1:22, "X", "Y", "MT")) %>%
  dplyr::filter(transcript_biotype == "protein_coding" |
                  transcript_biotype == "lncRNA")

counts_correct_gene_names <- table %>%
  filter(gene %in% canonical_transcripts$external_gene_name)

counts_incorrect_gene_names <- table %>%
  filter(!gene %in% canonical_transcripts$external_gene_name)

# Change incorrect genes names using official ones, when possible
counts_incorrect_gene_names <- counts_incorrect_gene_names %>%
  left_join(canonical_transcripts %>%
              dplyr::select(ensembl_gene_id, external_gene_name,
                            external_synonym),
            by = c("gene" = "external_synonym")) %>%
  dplyr::select(gene, ensembl_gene_id, external_gene_name, everything()) %>%
  filter(!is.na(external_gene_name)) %>%
  filter(!external_gene_name %in% counts_correct_gene_names$gene) %>%
  filter(!duplicated(external_gene_name)) %>%
  dplyr::select(-gene, -ensembl_gene_id) %>%
  dplyr::rename(gene = external_gene_name)

counts <- rbind(counts_correct_gene_names, counts_incorrect_gene_names)

mat <- as.matrix(counts[, -1])
rownames(mat) <- counts$gene
coldata <- data.frame(metadata[,-1], row.names = metadata$cell)

oocytes_sce <- SingleCellExperiment(assays =
                                    list(counts = mat[, rownames(coldata)]),
                                   colData = coldata)
oocytes_sce <- logNormCounts(oocytes_sce)

save(oocytes_sce, file = "../../eh_data/oocytes_sce.rda",
     compress = "xz",
     compression_level = 9)
