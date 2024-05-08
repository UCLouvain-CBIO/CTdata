## Code to prepare `mean_methylation_in_tissues` dataset goes here

library(GenomicRanges)
library(tidyverse)
library(SummarizedExperiment)

load("../extdata/all_genes.rda")
load("../../eh_data/methylation_in_tissues.rda")

## Promoter region is defined as `nt_up` nucleotides upstream TSS
## and `nt_down` nucleotides downstream TSS
nt_up <- 1000
nt_down <- 200

## Calculate mean methylation of each promoter in tissues
## and store CpG number by promoter
i <- 0
prom_mean_met_in_tissues <- tibble(tissue =
                                     c(colnames(methylation_in_tissues),
                                       "CpG_number"))

for (gene in all_genes$external_gene_name) {
  i <- i + 1
  print(i)
  TSS <- all_genes %>%
    filter(external_gene_name == gene) %>%
    pull(transcription_start_site)

  chr <- all_genes %>%
    filter(external_gene_name == gene) %>%
    pull(chromosome_name)

  strand <- all_genes %>%
    filter(external_gene_name == gene) %>%
    pull(strand)

  TSS <- all_genes %>%
    filter(external_gene_name == gene) %>%
    pull(transcription_start_site)

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

  promoter_methylation <- subsetByOverlaps(methylation_in_tissues,
                                           promoter_gr)
  tmp <- enframe(colMeans(assay(promoter_methylation), na.rm = TRUE),
                 name = "tissue", value = gene)
  # Store the number of CpG in the promoter region
  tmp <- rbind(tmp, c("CpG_number", dim(promoter_methylation)[1]))
  prom_mean_met_in_tissues <- left_join(prom_mean_met_in_tissues, tmp)
}

## Store CpG number by promoter
CT_CpG_number <- prom_mean_met_in_tissues %>%
  filter(tissue == "CpG_number") %>%
  pivot_longer(names_to = "external_gene_name", values_to = "CpG_number",
               -tissue) %>%
  dplyr::select(-tissue) %>%
  mutate(CpG_number = as.integer(CpG_number))

## Store mean methylation level by tissue by promoter
prom_mean_met_in_tissues <- prom_mean_met_in_tissues %>%
  filter(tissue != "CpG_number") %>%
  pivot_longer(names_to = "external_gene_name", values_to = "mean_methylation",
               -tissue) %>%
  mutate(mean_methylation = as.numeric(mean_methylation)) %>%
  pivot_wider(names_from = tissue, values_from = mean_methylation)
mat <- prom_mean_met_in_tissues %>%
  dplyr::select(-external_gene_name) %>%
  as.matrix()
rownames(mat) <- prom_mean_met_in_tissues$external_gene_name

## Calculate CpG densities and ratios of methylation in somatic tissues vs sperm
methylation_analysis <- tibble(
  external_gene_name = prom_mean_met_in_tissues$external_gene_name,
  somatic_met_level = prom_mean_met_in_tissues %>%
    dplyr::select(-c(external_gene_name, placenta, testis, sperm)) %>%
    dplyr::rowwise() %>%
    rowMeans(na.rm = TRUE),
  sperm_met_level = prom_mean_met_in_tissues %>%
    dplyr::select(sperm) %>%
    pull(sperm))

## CT Genes controlled by methylation should have somatic_methylation == TRUE
## and sperm_methylation == FALSE
methylation_analysis <- methylation_analysis %>%
  left_join(all_genes %>%
              dplyr::select(external_gene_name, ensembl_gene_id)) %>%
  mutate(ratio_somatic_sperm = somatic_met_level / sperm_met_level) %>%
  left_join(CT_CpG_number) %>%
  mutate(CpG_density = CpG_number / (nt_up + nt_down) * 100) %>%
  mutate(CpG_promoter = case_when(CpG_density < 2 ~ "low",
                                  CpG_density >= 2 &
                                    CpG_density < 4 ~ "intermediate",
                                  CpG_density >= 4 ~ "high")) %>%
  mutate(somatic_methylation =
           case_when(somatic_met_level < 50 ~ FALSE,
                     somatic_met_level >= 50 ~ TRUE)) %>%
  mutate(germline_methylation =
           case_when(ratio_somatic_sperm > 2 ~ FALSE,
                     ratio_somatic_sperm <= 2 ~ TRUE)) %>%
  column_to_rownames("external_gene_name") %>%
  dplyr::select(ensembl_gene_id, CpG_density,
                CpG_promoter, somatic_met_level, sperm_met_level,
                somatic_methylation, germline_methylation)

mean_methylation_in_tissues <-
  SummarizedExperiment(assays = mat,
                       rowData = methylation_analysis)

save(mean_methylation_in_tissues, file = "../../eh_data/mean_methylation_in_tissues.rda",
     compress = "xz",
     compression_level = 9)
load("../../eh_data/mean_methylation_in_tissues.rda")
