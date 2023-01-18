## Code to prepare `TCGA_CT_methylation` dataset goes here

library("TCGAbiolinks")
library("tidyverse")
library("SummarizedExperiment")
library("GenomicRanges")
library("BiocFileCache")
library("org.Hs.eg.db")

load("../extdata/CT_list.rda")

bfc <- BiocFileCache(cache = "/home/users/aloriot/.cache/BiocFileCache",
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

CT_promoter_gr <- makeGRangesFromDataFrame(
  CT_list %>%
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

## SKCM_CT_methylation
load(bfc[[bfcquery(bfc, "TCGA_SKCM_methylation")$rid]])
SKCM_CT_methylation <- subsetByOverlaps(data, CT_promoter_gr)

## LUAD_CT_methylation
load(bfc[[bfcquery(bfc, "TCGA_LUAD_methylation")$rid]])
LUAD_CT_methylation <- subsetByOverlaps(data, CT_promoter_gr)

## LUSC_CT_methylation
load(bfc[[bfcquery(bfc, "TCGA_LUSC_methylation")$rid]])
LUSC_CT_methylation <- subsetByOverlaps(data, CT_promoter_gr)

## COAD_CT_methylation
load(bfc[[bfcquery(bfc, "TCGA_COAD_methylation")$rid]])
COAD_CT_methylation <- subsetByOverlaps(data, CT_promoter_gr)

## ESCA_CT_methylation
load(bfc[[bfcquery(bfc, "TCGA_ESCA_methylation")$rid]])
ESCA_CT_methylation <- subsetByOverlaps(data, CT_promoter_gr)

## BRCA_CT_methylation
load(bfc[[bfcquery(bfc, "TCGA_BRCA_methylation")$rid]])
BRCA_CT_methylation <- subsetByOverlaps(data, CT_promoter_gr)

## HNSC_CT_methylation
load(bfc[[bfcquery(bfc, "TCGA_HNSC_methylation")$rid]])
HNSC_CT_methylation <- subsetByOverlaps(data, CT_promoter_gr)

met <- cbind(assay(SKCM_CT_methylation), assay(LUAD_CT_methylation),
             assay(LUSC_CT_methylation), assay(COAD_CT_methylation),
             assay(ESCA_CT_methylation), assay(BRCA_CT_methylation),
             assay(HNSC_CT_methylation))

colData(SKCM_CT_methylation)$sample <- substr(colData(SKCM_CT_methylation)$samples, 1, 16)
colData(SKCM_CT_methylation)$project_id <- "TCGA-SKCM"
colData(LUAD_CT_methylation)$sample <- substr(colData(LUAD_CT_methylation)$samples, 1, 16)
colData(LUAD_CT_methylation)$project_id <- "TCGA-LUAD"
colData(LUSC_CT_methylation)$sample <- substr(colData(LUSC_CT_methylation)$samples, 1, 16)
colData(LUSC_CT_methylation)$project_id <- "TCGA-LUSC"
colData(COAD_CT_methylation)$sample <- substr(colData(COAD_CT_methylation)$samples, 1, 16)
colData(COAD_CT_methylation)$project_id <- "TCGA-COAD"
colData(ESCA_CT_methylation)$sample <- substr(colData(ESCA_CT_methylation)$samples, 1, 16)
colData(ESCA_CT_methylation)$project_id <- "TCGA-ESCA"
colData(BRCA_CT_methylation)$sample <- substr(colData(BRCA_CT_methylation)$samples, 1, 16)
colData(BRCA_CT_methylation)$project_id <- "TCGA-BRCA"
colData(HNSC_CT_methylation)$sample <- substr(colData(HNSC_CT_methylation)$samples, 1, 16)
colData(HNSC_CT_methylation)$project_id <- "TCGA-HNSC"

coldata <- rbind(colData(SKCM_CT_methylation),
                 colData(LUAD_CT_methylation),
                 colData(LUSC_CT_methylation),
                 colData(COAD_CT_methylation),
                 colData(ESCA_CT_methylation),
                 colData(BRCA_CT_methylation),
                 colData(HNSC_CT_methylation))

TCGA_CT_methylation <- SummarizedExperiment(assays = list(methylation = met),
                                            colData = coldata,
                                            rowRanges = rowRanges(SKCM_CT_methylation))

save(TCGA_CT_methylation, file = "../../eh_data/TCGA_CT_methylation.rda",
     compress = "xz",
     compression_level = 9)
