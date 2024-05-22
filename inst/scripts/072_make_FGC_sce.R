## Code to prepare `FGC_sce` dataset goes here

library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(biomaRt)

## Data from paper: Single-cell roadmap of human gonadal development
## (Garcia-Alonso, Nature 2022)
## (https://www.nature.com/articles/s41586-022-04918-4)

## ee58527e-e1e4-465d-8dc8-800ee40f14f2.rds data dowloaded from
## https://cellxgene.cziscience.com/collections/661a402a-2a5a-4c71-9b05-b346c57bc451Data

seur <- readRDS("../extdata/ee58527e-e1e4-465d-8dc8-800ee40f14f2.rds")

FGC_sce <- SingleCellExperiment(
  assays = list(counts = as.matrix(seur@assays$RNA$counts)),
  colData = seur@meta.data)

## Normalization.
FGC_sce <- logNormCounts(FGC_sce)
levels(FGC_sce$sex) <- c("F", "M")

FGC_sce$type <- paste0(FGC_sce$sex, "_", FGC_sce$celltype)

# Keep coherently annotated cells
FGC_sce <- FGC_sce[, FGC_sce$type %in%
             c("F_PGC", "F_GC", "F_GC_mitotic",  "F_oogonia_STRA8",
               "F_oogonia_meiotic", "F_pre_oocyte", "F_oocyte",
               "M_PGC", "M_GC", "M_GC_mitotic", "M_pre_spermatogonia")]

## Merge cell types into 5 major germ cells states.
## See Supplementary Note 3 "Characterisation of germ cells" from paper.
FGC_sce$type[FGC_sce$type %in% c("F_GC", "F_GC_mitotic")] <- "F_GC"
FGC_sce$type[FGC_sce$type %in% c("M_GC", "M_GC_mitotic")] <- "M_GC"

FGC_sce$type[FGC_sce$type %in% c("F_oogonia_STRA8", "F_oogonia_meiotic")] <-
  "F_oogonia"
FGC_sce$type[FGC_sce$type %in% c("F_pre_oocyte", "F_oocyte")] <- "F_oocyte"

FGC_sce$type <- factor(FGC_sce$type, levels =
                     c("F_PGC", "F_GC","F_oogonia", "F_oocyte",
                       "M_PGC", "M_GC", "M_pre_spermatogonia"))

FGC_sce$stage <- "pre-meiotic"
FGC_sce[, FGC_sce$type %in%
      c("F_oocyte")]$stage <- "meiotic"
FGC_sce$germcell <- TRUE

## Change FGC_sce rownames (from ensembl_id to gene names,
## keeping only official gene names)
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

FGC_sce <- FGC_sce[rownames(FGC_sce) %in%
                     canonical_transcripts$ensembl_gene_id, ]
ensembl_gene <- as.data.frame(canonical_transcripts[,1:2]) %>% unique()
rownames(ensembl_gene) <- ensembl_gene$ensembl_gene_id
rownames(FGC_sce) <- ensembl_gene[rownames(FGC_sce),]$external_gene_name

save(FGC_sce, file = "../../eh_data/FGC_sce.rda",
     compress = "xz",
     compression_level = 9)
