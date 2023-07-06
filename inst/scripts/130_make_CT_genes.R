## Code to prepare `CT_genes` dataset goes here

library("tidyverse")
library("SummarizedExperiment")

load("../extdata/CT_list.rda")
load(file = "../../eh_data/TCGA_CT_methylation.rda")
load(file = "../../eh_data/CT_mean_methylation_in_tissues.rda")
load(file = "../../eh_data/CT_methylation_in_tissues.rda")

################################################################################
## Add CpG densities and promoter methylation analysis in normal tissues
################################################################################
CT_genes <- CT_list %>%
  left_join(as_tibble(rowData(CT_mean_methylation_in_tissues)))

################################################################################
## Based on DAC induction and on methylation levels in normal tissues
## (if available).For some genes, methylation analysis was not possible due
## to multimapping issues. In this case, genes are still considered as regulated
## by methylation if they show a strong activation in cells treated with
## 5-Aza-2′-Deoxycytidine.
################################################################################

CT_genes <- CT_genes %>%
  mutate(regulated_by_methylation =
           case_when(DAC_induced == FALSE ~ FALSE,
                     DAC_induced == TRUE & is.na(germline_methylation) ~ TRUE,
                     DAC_induced == TRUE &
                       somatic_methylation == TRUE &
                       germline_methylation == FALSE ~ TRUE,
                     DAC_induced == TRUE &
                       (somatic_methylation == FALSE |
                          germline_methylation == TRUE) ~ FALSE))

################################################################################
## Add X_linked information
################################################################################

CT_genes <- CT_genes %>%
  mutate(X_linked = ifelse(chromosome_name == "X", TRUE, FALSE))

################################################################################
## Flag genes as "oncogenic" or "tumor suppressor" using
## [Cancermine](http://bionlp.bcgsc.ca/cancermine/), a literature-mined database
## of drivers, oncogenes and tumor suppressors in cancer.
################################################################################
if (!file.exists("../extdata/cancermine_collated.tsv")) {
  download.file(url = "https://zenodo.org/record/7537824/files/cancermine_collated.tsv?download=1",
                destfile = "../extdata/cancermine_collated.tsv")
}
cancermine <- read_tsv("../extdata/cancermine_collated.tsv")
oncogenes <- cancermine %>%
  filter(role == "Oncogene") %>%
  pull(gene_normalized)
TS <- cancermine %>%
  filter(role == "Tumor_Suppressor") %>%
  pull(gene_normalized)
CT_genes <- CT_genes %>%
  mutate(oncogene = case_when(external_gene_name %in% oncogenes ~ "oncogene")) %>%
  mutate(tumor_suppressor = case_when(external_gene_name %in% TS ~ "tumor_suppressor"))

################################################################################
# Reorder the final table
################################################################################
CT_genes <- CT_genes %>%
  dplyr::rename(chr = chromosome_name) %>%
  dplyr::select("ensembl_gene_id", "external_gene_name", "family",
                "chr", "strand", "transcription_start_site", "X_linked",
                "TPM_testis", "max_TPM_somatic",
                "GTEX_category", "lowly_expressed_in_GTEX",
                "multimapping_analysis", "testis_specificity",
                "testis_cell_type", "Higher_in_somatic_cell_type",
                "percent_of_positive_CCLE_cell_lines",
                "percent_of_negative_CCLE_cell_lines", "max_TPM_in_CCLE",
                "CCLE_category", "percent_pos_tum", "percent_neg_tum",
                "max_TPM_in_TCGA", "TCGA_category", "DAC_induced",
                "somatic_met_level", "sperm_met_level", "somatic_methylation",
                "germline_methylation", "regulated_by_methylation", "CpG_density",
                "CpG_promoter", "external_transcript_name",
                "ensembl_transcript_id", "transcript_biotype", "oncogene",
                "tumor_suppressor")

save(CT_genes, file = "../../eh_data/CT_genes.rda",
     compress = "xz",
     compression_level = 9)



