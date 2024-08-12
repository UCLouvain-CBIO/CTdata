## Code to prepare `mean_methylation_in_hESC` dataset goes here

library(GenomicRanges)
library(tidyverse)
library(SummarizedExperiment)

load("../../eh_data/all_genes.rda")
load("../../eh_data/methylation_in_hESC.rda")
load("../../../CTdata_extdata/hESC_RNAseq_coldata.rda")

## Promoter region is defined as `nt_up` nucleotides upstream TSS
## and `nt_down` nucleotides downstream TSS
nt_up <- 1000
nt_down <- 200


## Calculate mean methylation of each promoter in cell types
## and store CpG number by promoter
prom_mean_met_in_hESC <- tibble(cell = c(colnames(methylation_in_hESC)))

for (gene in all_genes$external_gene_name) {
  print(which(all_genes$external_gene_name == gene))
  TSS <- all_genes %>%
    filter(external_gene_name == gene) %>%
    pull(transcription_start_site)

  chr <- all_genes %>%
    filter(external_gene_name == gene) %>%
    pull(chr)

  strand <- all_genes %>%
    filter(external_gene_name == gene) %>%
    pull(strand)

  if (strand == 1) { # Analyse region at +/- nt_up and nt_down around TSS
    promoter_gr <- GRanges(seqnames = paste0("chr", chr),
                           strand = '+',
                           ranges = IRanges(start = TSS - nt_up,
                                            end = TSS + nt_down))
    promoter_gr$TSS <- TSS
  }

  if (strand == -1) { # Analyse region at +/- nt_up and nt_down around TSS
    promoter_gr <- GRanges(seqnames = paste0("chr", chr),
                           strand = '-',
                           ranges = IRanges(start = TSS - nt_down,
                                            end = TSS + nt_up))
    promoter_gr$TSS <- TSS
  }

  promoter_methylation <- subsetByOverlaps(methylation_in_hESC,
                                           promoter_gr)
  tmp <- enframe(colMeans(assay(promoter_methylation), na.rm = TRUE),
                 name = "cell", value = gene)

  prom_mean_met_in_hESC <- left_join(prom_mean_met_in_hESC, tmp)
}


## Store mean methylation level by cell type by promoter
prom_mean_met_in_hESC <- prom_mean_met_in_hESC %>%
  pivot_longer(names_to = "external_gene_name", values_to = "mean_methylation",
               -cell) %>%
  mutate(mean_methylation = as.numeric(mean_methylation)) %>%
  pivot_wider(names_from = cell, values_from = mean_methylation)

################################################################################
## Add coldata and save as a SE
################################################################################

mat <- prom_mean_met_in_hESC %>%
  dplyr::select(-external_gene_name) %>%
  as.matrix()
rownames(mat) <- prom_mean_met_in_hESC$external_gene_name

coldata <- as.data.frame(coldata) %>%
  mutate(cell_type = "hESC") %>%
  mutate(genotype = factor(genotype, levels = c("XX", "XY"))) %>%
  mutate(dataset = "Encode") %>%
  mutate(data = "Methylation") %>%
  dplyr::select(-files) %>%
  column_to_rownames("sample")

rowdata <- all_genes %>%
  dplyr::select(external_gene_name, ensembl_gene_id) %>%
  column_to_rownames("external_gene_name")

mean_methylation_in_hESC <- SummarizedExperiment(assays = mat,
                                                 colData =
                                                   coldata[colnames(mat), ],
                                                 rowData = rowdata)

save(mean_methylation_in_hESC,
     file = "../../eh_data/mean_methylation_in_hESC.rda",
     compress = "xz",
     compression_level = 9)

