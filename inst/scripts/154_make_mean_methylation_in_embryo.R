## Code to prepare `mean_methylation_in_embryo` dataset goes here

library("tidyverse")
library("GenomicRanges")
library("SummarizedExperiment")

load("../extdata/all_promoter_regions_hg19.rda")
load("../../eh_data/methylation_in_embryo.rda")
load("../../eh_data/all_genes.rda")



################################################################################
# Define promoter regions
################################################################################

## Promoter region is defined as `nt_up` nucleotides upstream TSS
## and `nt_down` nucleotides downstream TSS
nt_up <- 1000
nt_down <- 500

all_promoter_regions <- as_tibble(all_promoter_regions_hg19) %>%
  dplyr::select(external_gene_name, chr, TSS_liftover, strand) %>%
  mutate(start = case_when(strand == "+" ~ TSS_liftover - nt_up,
                           strand == "-" ~ TSS_liftover - nt_down)) %>%
  mutate(stop = case_when(strand == "+" ~ TSS_liftover + nt_down,
                          strand == "-" ~ TSS_liftover + nt_up)) %>%
  mutate(chromosome_name = paste0("chr", chr))

all_prom_GR <- makeGRangesFromDataFrame(all_promoter_regions,
                                        keep.extra.columns = TRUE,
                                        seqnames.field = "chromosome_name",
                                        start.field = "start",
                                        end.field = "stop")

gene_ids <- all_genes %>%
  dplyr::select("ensembl_gene_id", "external_gene_name") %>%
  filter(external_gene_name %in% all_prom_GR$external_gene_name)

all_prom_GR$ensembl_gene_id <- gene_ids$ensembl_gene_id


################################################################################
# Compute mean meth for all promoters
################################################################################

mean_meth <- tibble(cell = colnames(methylation_in_embryo))

for (i in 1:length(all_prom_GR)) {
  print(i)
  tested_gene_GR <- subsetByOverlaps(methylation_in_embryo, all_prom_GR[i])
  mat <- assay(tested_gene_GR)
  tmp <- enframe(colMeans(mat, na.rm = TRUE), name = "cell",
                 value = all_prom_GR[i]$external_gene_name)
  mean_meth <- suppressMessages(left_join(mean_meth, tmp))
}


mat <- t(as.matrix(mean_meth[,-1]))
colnames(mat) <- mean_meth$cell

################################################################################
# Save it as a RangedSE
################################################################################

mean_methylation_in_embryo <-
  SummarizedExperiment(assays = mat ,
                       colData = colData(methylation_in_embryo),
                       rowRanges = all_prom_GR)

save(mean_methylation_in_embryo,
     file = "../../eh_data/mean_methylation_in_embryo.rda",
     compress = "xz",
     compression_level = 9)
