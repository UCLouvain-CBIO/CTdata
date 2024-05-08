## Code to prepare `all_genes` and `CT_genes` dataset goes here

library("tidyverse")
library("SummarizedExperiment")

load("../extdata/all_genes_prelim.rda")
load(file = "../../eh_data/TCGA_methylation.rda")
load(file = "../../eh_data/mean_methylation_in_tissues.rda")
load(file = "../../eh_data/methylation_in_tissues.rda")
load(file = "../../eh_data/DAC_treated_cells_multimapping.rda")

################################################################################
## Add DAC_induced column specifying if genes are induced by
## 5-Aza-2′-Deoxycytidine
################################################################################
induced <- as_tibble(rowData(DAC_treated_cells_multimapping),
                     rownames = "ensembl_gene_id") %>%
  filter(induced) %>%
  pull(external_gene_name)

not_induced <- as_tibble(rowData(DAC_treated_cells_multimapping),
                         rownames = "ensembl_gene_id") %>%
  filter(!induced) %>%
  pull(external_gene_name)

unclear <- as_tibble(rowData(DAC_treated_cells_multimapping),
                     rownames = "ensembl_gene_id") %>%
  filter(is.na(induced)) %>%
  pull(external_gene_name)

all_genes <- all_genes_prelim %>%
  mutate(DAC_induced = case_when(external_gene_name %in% induced ~ TRUE,
                                 external_gene_name %in% not_induced ~ FALSE,
                                 external_gene_name %in% unclear ~ NA))


################################################################################
## Add CpG densities and promoter methylation analysis in normal tissues
################################################################################
all_genes <- all_genes %>%
  left_join(as_tibble(rowData(mean_methylation_in_tissues)))

################################################################################
## Based on DAC induction and on methylation levels in normal tissues
## (if available).For some genes, methylation analysis was not possible due
## to multimapping issues. In this case, genes are still considered as regulated
## by methylation if they show a strong activation in cells treated with
## 5-Aza-2′-Deoxycytidine.
################################################################################

all_genes <- all_genes %>%
  mutate(regulated_by_methylation =
           case_when((somatic_methylation | is.na(somatic_methylation)) &
                       DAC_induced ~ TRUE,
                     !somatic_methylation ~ FALSE,
                     !DAC_induced ~ FALSE,
                     is.na(DAC_induced) ~ NA))

################################################################################
## Add X_linked information
################################################################################

all_genes <- all_genes %>%
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
all_genes <- all_genes %>%
  mutate(oncogene = case_when(external_gene_name %in% oncogenes ~
                                "oncogene")) %>%
  mutate(tumor_suppressor = case_when(external_gene_name %in% TS ~
                                        "tumor_suppressor"))

################################################################################
# Reorder the final table
################################################################################
all_genes <- all_genes %>%
  dplyr::rename(chr = chromosome_name) %>%
  dplyr::select("ensembl_gene_id", "external_gene_name", "CT_gene_type",
                "testis_specificity", "regulated_by_methylation", "X_linked",
                "chr", "strand", "transcript_start", "transcript_end",
                "transcription_start_site", everything())

save(all_genes, file = "../../eh_data/all_genes.rda",
     compress = "xz",
     compression_level = 9)

CT_genes <- all_genes %>%
  filter(CT_gene_type != "other") %>%
  arrange(desc(CT_gene_type))

save(CT_genes, file = "../../eh_data/CT_genes.rda",
     compress = "xz",
     compression_level = 9)



