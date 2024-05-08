## Code to prepare `TCGA_CT_methylation` dataset goes here

library("TCGAbiolinks")
library("tidyverse")
library("SummarizedExperiment")
library("GenomicRanges")
library("BiocFileCache")
library("org.Hs.eg.db")

load("../extdata/all_genes.rda")

bfc <- BiocFileCache(cache = "../BiocFileCache",
                     ask = FALSE)

for(tumor_code in c("SKCM", "LUAD", "LUSC", "COAD", "ESCA", "BRCA", "HNSC")) {
  rname <- paste0("TCGA_", tumor_code, "_methylation")

  if(length(bfcquery(bfc, rname)$rid) == 0) {

    savepath <- bfcnew(bfc, rname, ext=".RData")

    # Query data
    query_gene_counts_harmonized <- GDCquery(project = paste0("TCGA-", tumor_code),
                                             data.category = "DNA Methylation",
                                             data.type = "Methylation Beta Value",
                                             platform = "Illumina Human Methylation 450",
                                             legacy = FALSE)

    # Download files
    GDCdownload(query_gene_counts_harmonized, method = "api",
                directory = paste0("./tmp/TCGA_", tumor_code, "_methylation"),
                files.per.chunk = 10)

    # Create a Summarized-Experiment object
    GDCprepare(query_gene_counts_harmonized,
               directory = paste0("./tmp/TCGA_", tumor_code, "_methylation"),
               save = TRUE,
               save.filename = savepath,
               remove.files.prepared = FALSE)
  }
}

## Promoter region is defined as `nt_up` nucleotides upstream TSS
## and `nt_down` nucleotides downstream TSS
nt_up <- 1000
nt_down <- 200

promoter_gr <- makeGRangesFromDataFrame(
  all_genes %>%
    dplyr::select(ensembl_gene_id, external_gene_name, external_transcript_name,
                  chromosome_name, strand, transcription_start_site) %>%
    mutate(chromosome_name = paste0("chr", chromosome_name)) %>%
    mutate(strand = if_else(strand == 1, '+', '-')) %>%
    mutate(start = case_when(strand == '+' ~ transcription_start_site - nt_up,
                             strand == '-' ~ transcription_start_site - nt_down)) %>%
    mutate(stop = case_when(strand == '+' ~ transcription_start_site + nt_down,
                            strand == '-' ~ transcription_start_site + nt_up)),
  keep.extra.columns = TRUE,
  seqnames.field = "chromosome_name",
  start.field = "start",
  end.field = "stop")

## SKCM_methylation
load(bfc[[bfcquery(bfc, "TCGA_SKCM_methylation")$rid]])
SKCM_methylation <- subsetByOverlaps(data, promoter_gr)

## LUAD_methylation
load(bfc[[bfcquery(bfc, "TCGA_LUAD_methylation")$rid]])
LUAD_methylation <- subsetByOverlaps(data, promoter_gr)

## LUSC_methylation
load(bfc[[bfcquery(bfc, "TCGA_LUSC_methylation")$rid]])
LUSC_methylation <- subsetByOverlaps(data, promoter_gr)

## COAD_methylation
load(bfc[[bfcquery(bfc, "TCGA_COAD_methylation")$rid]])
COAD_methylation <- subsetByOverlaps(data, promoter_gr)

## ESCA_methylation
load(bfc[[bfcquery(bfc, "TCGA_ESCA_methylation")$rid]])
ESCA_methylation <- subsetByOverlaps(data, promoter_gr)

## BRCA_methylation
load(bfc[[bfcquery(bfc, "TCGA_BRCA_methylation")$rid]])
BRCA_methylation <- subsetByOverlaps(data, promoter_gr)

## HNSC_methylation
load(bfc[[bfcquery(bfc, "TCGA_HNSC_methylation")$rid]])
HNSC_methylation <- subsetByOverlaps(data, promoter_gr)

met <- cbind(assay(SKCM_methylation), assay(LUAD_methylation),
             assay(LUSC_methylation), assay(COAD_methylation),
             assay(ESCA_methylation), assay(BRCA_methylation),
             assay(HNSC_methylation))

colData(SKCM_methylation)$sample <- substr(colData(SKCM_methylation)$samples, 1, 16)
colData(SKCM_methylation)$project_id <- "TCGA-SKCM"
colData(LUAD_methylation)$sample <- substr(colData(LUAD_methylation)$samples, 1, 16)
colData(LUAD_methylation)$project_id <- "TCGA-LUAD"
colData(LUSC_methylation)$sample <- substr(colData(LUSC_methylation)$samples, 1, 16)
colData(LUSC_methylation)$project_id <- "TCGA-LUSC"
colData(COAD_methylation)$sample <- substr(colData(COAD_methylation)$samples, 1, 16)
colData(COAD_methylation)$project_id <- "TCGA-COAD"
colData(ESCA_methylation)$sample <- substr(colData(ESCA_methylation)$samples, 1, 16)
colData(ESCA_methylation)$project_id <- "TCGA-ESCA"
colData(BRCA_methylation)$sample <- substr(colData(BRCA_methylation)$samples, 1, 16)
colData(BRCA_methylation)$project_id <- "TCGA-BRCA"
colData(HNSC_methylation)$sample <- substr(colData(HNSC_methylation)$samples, 1, 16)
colData(HNSC_methylation)$project_id <- "TCGA-HNSC"

coldata <- rbind(colData(SKCM_methylation),
                 colData(LUAD_methylation),
                 colData(LUSC_methylation),
                 colData(COAD_methylation),
                 colData(ESCA_methylation),
                 colData(BRCA_methylation),
                 colData(HNSC_methylation))

TCGA_methylation <- SummarizedExperiment(assays = list(methylation = met),
                                            colData = coldata,
                                            rowRanges = rowRanges(SKCM_methylation))

save(TCGA_methylation, file = "../../eh_data/TCGA_methylation.rda",
     compress = "xz",
     compression_level = 9)
