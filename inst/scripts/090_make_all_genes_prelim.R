## Code to prepare `all_genes_prelim` dataset goes here

library("tidyverse")
library("SummarizedExperiment")
library("biomaRt")

load(file = "../../eh_data/GTEX_data.rda")
load(file = "../../eh_data/CCLE_data.rda")
load(file = "../../eh_data/normal_tissues_multimapping_data.rda")
load(file = "../../eh_data/HPA_cell_type_specificities.rda")
load(file = "../../eh_data/DAC_treated_cells_multimapping.rda")
load(file = "../../eh_data/TCGA_TPM.rda")

################################################################################
## Start from rowData(GTEX), summarizing the tissue specificity category
## assigned to each gene using GTEx database
################################################################################
all_genes_prelim <- as_tibble(rowData(GTEX_data), rownames = "ensembl_gene_id")
all_genes_prelim$TPM_testis <- assay(GTEX_data)[, "Testis"]


################################################################################
## Add multimapping_analysis column from normal_tissues_multimapping_data,
## assessing testis-specificity of genes flagged as "lowly_expressed" in
## GTEX_category from GTEX_data (when their TPM < 1 in testis)
################################################################################
all_genes_prelim <- all_genes_prelim %>%
  left_join(as_tibble(rowData(normal_tissues_multimapping_data),
                      rownames = "ensembl_gene_id"))


################################################################################
## Add HPA_cell_type_specificity (based on data from the Human Protein Atlas) to
## add cell type specificities and flag genes expressed in a somatic cell type
################################################################################
all_genes_prelim <- all_genes_prelim %>%
  left_join(HPA_cell_type_specificities)


################################################################################
## Add info from rowData(CCLE_data), summarising the analysis of CCLE
## database. In CCLE_category, genes are tagged as "activated" when
## at least 20 % of cell lines are negative (TPM <= 0.5)
## and at least 1% of cell line is highly positive (TPM >= 1)
## and and expressed at TPM > 5 in at least one cell line
## (to avoid leaky genes)
################################################################################
all_genes_prelim <- all_genes_prelim %>%
  left_join(as_tibble(rowData(CCLE_data), rownames = "ensembl_gene_id"))
all_genes_prelim[is.na(all_genes_prelim$CCLE_category),
                 "CCLE_category"] <- "not_in_CCLE"


################################################################################
## Add info from rowData(TCGA_TPM), summarizing the analysis of TCGA
## tumor samples. In TCGA_category, genes are tagged as "activated" when
## at least 20 % of tumors are negative (TPM <= 0.5), at least 1% of tumors
## are positive (TPM >= 1) and expressed at TPM > 5 in at least one tumor
## (to avoid leaky genes). Genes that were flagged as testis-specific in
## `multimapping analysis` are flagged as "mulimapping_issue" in TCGA_category,
## as most of them  are not detected in TCGA database
## (we cannot filter on their expression in TCGA dataset).
################################################################################
all_genes_prelim <- all_genes_prelim %>%
  left_join(as_tibble(rowData(TCGA_TPM), rownames = "ensembl_gene_id") %>%
              dplyr::select(ensembl_gene_id, external_gene_name,
                            percent_pos_tum, percent_neg_tum, max_TPM_in_TCGA,
                            max_q75_in_NT, TCGA_category))

all_genes_prelim[all_genes_prelim$lowly_expressed_in_GTEX == TRUE ,
                 "TCGA_category"] <- "multimapping_issue"

################################################################################
## Add testis_specificity summarizing testis-specificity analysis from
## GTEX, multimapping, HPA_cell_type_specificities, CCLE and TCGA categories
## Testis-specific genes must be:
## - classified as "testis-specific" in GTEX or multimappig analysis
## - must not be detected in any somatic cell type (in HPA scRNAseq analysis)
## - can not be classified as "leaky" in CLLE and TCGA categories
## - can not be detected (TPM > 0.5) in more than 75% of normal peritumoral
## tissues from a TCGA tumor type
## Testis_preferential genes must be:
## - classified as "testis_preferential" in GTEX or multimappig analysis
## - or initially classified as "testis-specific" in GTEX or multimappig
## analysis but downgraded to "testis_preferential because they were detected
## in a somatic cell type (in HPA scRNAseq analysis)
## - can be detected in some somatic cell type but the level has to be at
## least 10 times lower than the level detected in a germ cell type
## (in HPA scRNAseq analysis)
## - must not be classified as "leaky" in CLLE and TCGA categories
################################################################################

all_genes_prelim <- all_genes_prelim %>%
  mutate(testis_specificity = case_when(
    (GTEX_category == "testis_specific" |
       multimapping_analysis == "testis_specific") &
      (is.na(not_detected_in_somatic_HPA) | not_detected_in_somatic_HPA) &
      CCLE_category != "leaky" &
      TCGA_category != "leaky" &
      (TCGA_category == "multimapping_issue" |
         max_q75_in_NT < 0.5) ~ 'testis_specific',
    # Downgrade TS to TP based on HPA (if detected in somatic cell type but only
    # if it is detected at a level lower than 10 times than in a germ cell type)
    (GTEX_category == "testis_specific" |
       multimapping_analysis == "testis_specific") &
      !not_detected_in_somatic_HPA &
      HPA_ratio_germ_som >= 10 &
      CCLE_category != "leaky" &
      TCGA_category != "leaky" ~ 'testis_preferential',
    (GTEX_category == "testis_preferential" |
       multimapping_analysis == "testis_preferential") &
      (is.na(HPA_ratio_germ_som) | HPA_ratio_germ_som >= 10) &
      CCLE_category != "leaky" &
      TCGA_category != "leaky"  ~ "testis_preferential"))

all_genes_prelim$testis_specificity[
  is.na(all_genes_prelim$testis_specificity)] <- "not_testis_specific"

################################################################################
## create CT_gene_type column, specifiying if a gene is
## - a "Cancer testis gene" (CT_gene): testis-specific genes activated in CCLE
## cell lines and TCGA tumor samples.
## - or a "Cancer testis-preferential gene" (CTP_gene) (testis-preferential
## genes activated in CCLE cell lines and TCGA tumor samples).
################################################################################
all_genes_prelim <- all_genes_prelim %>%
  mutate(CT_gene_type = case_when(
    testis_specificity == "testis_specific" & CCLE_category == "activated" &
      (TCGA_category == "multimapping_issue" | TCGA_category == "activated")
    ~ "CT_gene",
    testis_specificity == "testis_preferential" & CCLE_category == "activated" &
      (TCGA_category == "multimapping_issue" | TCGA_category == "activated")
    ~ "CTP_gene"))

all_genes_prelim$CT_gene_type[is.na(all_genes_prelim$CT_gene_type)] <- "other"


################################################################################
## Associate each gene to its most biologically relevant transcript,
## selecting the canonical transcript from ensembl database.
## Add the TSS coordinates
################################################################################
ensembl <- biomaRt::useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
attributes_vector <- c("ensembl_gene_id",
                       "external_gene_name",
                       "ensembl_transcript_id",
                       "external_transcript_name",
                       "chromosome_name",
                       "strand",
                       "transcript_start",
                       "transcript_end",
                       "transcription_start_site",
                       "transcript_length",
                       "transcript_biotype",
                       "transcript_is_canonical")
transcripts_infos <- as_tibble(biomaRt::getBM(attributes = attributes_vector,
                                              mart = ensembl))
canonical_transcripts <- transcripts_infos %>%
  filter(transcript_is_canonical == 1) %>%
  filter(chromosome_name %in% c(1:22, "X", "Y", "MT")) %>%
  filter(transcript_biotype == "protein_coding" |
           transcript_biotype == "lncRNA")

all_genes_prelim <- all_genes_prelim %>%
  left_join(canonical_transcripts %>%
              dplyr::select(ensembl_gene_id,
                            external_transcript_name, ensembl_transcript_id,
                            chromosome_name, strand, transcript_start,
                            transcript_end, transcription_start_site,
                            transcript_length, transcript_biotype))

## Rm the few duplicated external genes (with more than one canonical
## transcript)
all_genes_prelim <- all_genes_prelim %>%
  filter(!duplicated(external_gene_name))


################################################################################
## Check the backbone of all selected CT genes by visualising RNA-seq reads
## from a testis sample using the Genome Integrative Viewer (IGV).
## The aim was initially to identify precisely the transcription start site of
## each gene, but unexpectedly we observed that for some genes, the reads were
## not properly aligned on exons, but were instead spread across a wide genomic
## region spanning the genes. These genes, flagged as "unclear" in IGV_backbone
## column, were removed from the CT_gene category
## as their expression values in GTEX, TCGA and CCLE might reflect a poorly
## defined transcription in these regions and are hence likely unreliable.
################################################################################

all_genes_prelim$IGV_backbone <- NA
all_genes_prelim[all_genes_prelim$external_gene_name == "MBD3L2",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "NKX2-4",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "FOXR2",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "PSG3",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "DMP1",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "ZNF679",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "OR8G5",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "VRTN",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "MBD3L5",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "SPANXN1",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "ROR1-AS1",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC01344",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "SMYD3-AS1",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC01913",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "KIAA2012-AS1",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC01980",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC01811",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "IGF2BP2-AS1",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC02020",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "SNHG27",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC02492",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC02434",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC02113",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC02522",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "SIM1-AS1",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC00326",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "ST7-OT4",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC01517",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC00867",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC02742",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "DLG2-AS2",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC02378",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "MRPS35-DT",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC02156",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC00392",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC02293",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "G2E3-AS1",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC00221",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC02152",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC02184",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC01925",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC01895",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC01029",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "ZBTB46-AS1",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC01203",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "FAM197Y7",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "FAM197Y6",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC00112",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "ARHGAP29-AS1",
                 "IGV_backbone"] <- "unclear"


all_genes_prelim[all_genes_prelim$external_gene_name == "LINC01804",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "DLX2-DT",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC02050",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC01208",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC01019",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC02211",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC02109",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC02531",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "FIBCD1-AS1",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC02646",
                 "IGV_backbone"] <- "unclear"

all_genes_prelim[all_genes_prelim$external_gene_name == "LINC02700",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC02327",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC01583",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "TP53TG3D",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "SMIM47",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "KRTAP20-4",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC01665",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "XKR3",
                 "IGV_backbone"] <- "unclear"





################################################################################
## idem but for CTP_genes
################################################################################
all_genes_prelim[all_genes_prelim$external_gene_name == "IL5",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "MYCNOS",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LEF1-AS1",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC02367",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC02076",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC00664",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "UMODL1",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "P2RX3",
                 "IGV_backbone"] <- "unclear"

all_genes_prelim[all_genes_prelim$external_gene_name == "GS1-24F4.2",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "USP2-AS1",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LINC00943",
                 "IGV_backbone"] <- "unclear"
all_genes_prelim[all_genes_prelim$external_gene_name == "LKAAEAR1",
                 "IGV_backbone"] <- "unclear"

all_genes_prelim[!is.na(all_genes_prelim$IGV_backbone) &
                   all_genes_prelim$IGV_backbone == "unclear",
                 "CT_gene_type"] <- "other"

################################################################################
## Add gene family column
################################################################################

all_genes_prelim$family <- NA
all_genes_prelim[grep(pattern = "MAGE", x = all_genes_prelim$external_gene_name,
                      value = FALSE), "family"] <- "MAGE"
all_genes_prelim[grep(pattern = "CT45A",
                      x = all_genes_prelim$external_gene_name,
                      value = FALSE), "family"] <- "CT45A"
all_genes_prelim[grep(pattern = "CTAG", x = all_genes_prelim$external_gene_name,
                      value = FALSE), "family"] <- "CTAG"
all_genes_prelim[grep(pattern = "GAGE", x = all_genes_prelim$external_gene_name,
                      value = FALSE), "family"] <- "GAGE"
all_genes_prelim[grep(pattern = "XAGE", x = all_genes_prelim$external_gene_name,
                      value = FALSE), "family"] <- "GAGE"
all_genes_prelim[grep(pattern = "PAGE", x = all_genes_prelim$external_gene_name,
                      value = FALSE), "family"] <- "PAGE"
all_genes_prelim[grep(pattern = "SPANX",
                      x = all_genes_prelim$external_gene_name,
                      value = FALSE), "family"] <- "SPANX"
all_genes_prelim[grep(pattern = "SSX", x = all_genes_prelim$external_gene_name,
                      value = FALSE), "family"] <- "SSX"
all_genes_prelim[grep(pattern = "VCX", x = all_genes_prelim$external_gene_name,
                      value = FALSE), "family"] <- "VCX"
all_genes_prelim[grep(pattern = "TSPY", x = all_genes_prelim$external_gene_name,
                      value = FALSE), "family"] <- "TSPY"
all_genes_prelim[grep(pattern = "RBMY", x = all_genes_prelim$external_gene_name,
                      value = FALSE), "family"] <- "RBMY"

save(all_genes_prelim, file = "../extdata/all_genes_prelim.rda",
     compress = "xz",
     compression_level = 9)
