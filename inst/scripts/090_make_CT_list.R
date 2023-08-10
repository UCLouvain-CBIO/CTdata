## Code to prepare `CT_list` dataset goes here

library("tidyverse")
library("SummarizedExperiment")
library("biomaRt")

load(file = "../../eh_data/GTEX_data.rda")
load(file = "../../eh_data/CCLE_data.rda")
load(file = "../../eh_data/normal_tissues_multimapping_data.rda")
load(file = "../../eh_data/testis_sce.rda")
load(file = "../../eh_data/scRNAseq_HPA.rda")
load(file = "../../eh_data/DAC_treated_cells_multimapping.rda")
load(file = "../../eh_data/TCGA_TPM.rda")

################################################################################
## Start from rowData(GTEX), summarizing the tissue specificity category
## assigned to each gene using GTEx database
################################################################################
all_genes <- as_tibble(rowData(GTEX_data), rownames = "ensembl_gene_id")
all_genes$TPM_testis <- assay(GTEX_data)[, "Testis"]

################################################################################
## Add multimapping_analysis column from normal_tissues_multimapping_data,
## assessing testis-specificity of genes flagged as "lowly_expressed" in
## GTEX_category from GTEX_data (when their TPM < 1 in testis)
################################################################################
all_genes <- all_genes %>%
  left_join(as_tibble(rowData(normal_tissues_multimapping_data),
                      rownames = "ensembl_gene_id"))

################################################################################
## Add testis_specificity summarizing testis-specificity analysis from
## GTEX and multimapping
################################################################################
all_genes <- all_genes %>%
  mutate(testis_specificity = case_when(
    GTEX_category == "testis_specific" |
      multimapping_analysis == "testis_specific" ~ "testis_specific",
    GTEX_category == "testis_preferential" ~ "testis_preferential"))

##########################################################################
## Add the `testis_cell_type` column (based on testis scRNAseq data),
## specifying the testis cell-type showing the highest mean expression
## of each gene.
## Remove genes for which testis_cell_type corresponds to a testis somatic
## cell type.
##########################################################################
all_genes <- all_genes %>%
  left_join(as_tibble(rowData(testis_sce)) %>%
              dplyr::select(external_gene_name, testis_cell_type))

##########################################################################
## Add the `Higher_in_somatic_cell_type` column (based on scRNAseq data
## of normal tissues from the Human Protein Atlas),
## specifying if some somatic cell types express the genes at a higher
## level than the level found in germ cell types, and filtered the CT_genes
## table accordingly.
##########################################################################
all_genes <- all_genes %>%
  left_join(as_tibble(rowData(scRNAseq_HPA), rownames = "ensembl_gene_id") %>%
              dplyr::select(ensembl_gene_id, external_gene_name,
                            Higher_in_somatic_cell_type))

################################################################################
## Add info from rowData(CCLE_data), summarising the analysis of CCLE
## database. In CCLE_category, genes are tagged as "activated" when
## at least 20 % of cell lines are negative (TPM <= 0.1)
## and at least one cell line is highly positive (TPM >= 10)
################################################################################
all_genes <- all_genes %>%
  left_join(as_tibble(rowData(CCLE_data), rownames = "ensembl_gene_id"))
all_genes[is.na(all_genes$CCLE_category), "CCLE_category"] <- "not_in_CCLE"

################################################################################
## Add info from rowData(TCGA_TPM), summarising the analysis of TCGA
## tumor samples. In TCGA_category, genes are tagged as "activated" when
## at least 20 % of tumors are negative (TPM <= 0.1)
## and at least one tumor is highly positive (TPM >= 10). Genes that were
## flagged as testis-specific in `multimapping analysis` are flagged as
## "mulimapping_issue" in TCGA_category, as most of them  are not detected
## in TCGA database.
################################################################################

all_genes <- all_genes %>%
  left_join(as_tibble(rowData(TCGA_TPM), rownames = "ensembl_gene_id") %>%
              dplyr::select(ensembl_gene_id, external_gene_name, percent_pos_tum,
                            percent_neg_tum, max_TPM_in_TCGA, TCGA_category))

all_genes[all_genes$lowly_expressed_in_GTEX == TRUE &
            all_genes$multimapping_analysis == "testis_specific",
          "TCGA_category"] <- "multimapping_issue"

################################################################################
## Add DAC column specifying if genes are induced by 5-Aza-2′-Deoxycytidine
################################################################################
induced <- as_tibble(rowData(DAC_treated_cells_multimapping),
                     rownames = "ensembl_gene_id") %>%
  filter(induced == TRUE) %>%
  pull(external_gene_name)

all_genes <- all_genes %>%
  mutate(DAC_induced = case_when(external_gene_name %in% induced ~ TRUE,
                                 !external_gene_name %in% induced ~ FALSE))

################################################################################
## Associate each gene to its likely most biologically relevant transcript,
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
  filter(transcript_biotype == "protein_coding" | transcript_biotype == "lncRNA")
all_genes <- all_genes %>%
  left_join(canonical_transcripts %>%
              dplyr::select(ensembl_gene_id,
                            external_transcript_name, ensembl_transcript_id,
                            chromosome_name, strand, transcript_start,
                            transcript_end, transcription_start_site,
                            transcript_length, transcript_biotype))

################################################################################
## create CT_list keeping testis-specific genes activated in CCLE cell lines
## and TCGA tumors.
################################################################################

CT_list <- all_genes %>%
  filter(testis_specificity == "testis_specific" |
           testis_specificity == "testis_preferential") %>%
  filter(!testis_cell_type %in% c( "Macrophage", "Endothelial", "Myoid",
                                   "Sertoli", "Leydig")) %>%
  filter(is.na(Higher_in_somatic_cell_type) | Higher_in_somatic_cell_type == FALSE) %>%
  filter(CCLE_category == "activated") %>%
  filter(TCGA_category == "activated" | TCGA_category == "multimapping_issue")

################################################################################
## For each CT gene, the most relevant transcript was validated manually, by
## visualizing on IGV RNA-seq alignments from a testis sample.
## When reads didn't fit to a referenced transcript in the testis sample,
## external_transcript_name was set to NA.
## When the testis-transcript is not the one activated in tumor cells lines,
## external_transcript_name was set to NA.
################################################################################

CT_list[CT_list$external_gene_name == "FAM81B", "external_transcript_name"] <- "FAM81B-205"
CT_list[CT_list$external_gene_name == "LIN28B", "external_transcript_name"] <- "LIN28B-202"
CT_list[CT_list$external_gene_name == "CHRNB4", "external_transcript_name"] <- "CHRNB4-204"
CT_list[CT_list$external_gene_name == "OR3A2", "external_transcript_name"] <- "OR3A2-205"
CT_list[CT_list$external_gene_name == "TAF7L", "external_transcript_name"] <- "TAF7L-203"
CT_list[CT_list$external_gene_name == "MBD3L2", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "NKX2-4", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "FOXR2", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "MMP20", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "CCDC196", "external_transcript_name"] <- "CCDC196-203"
CT_list[CT_list$external_gene_name == "CIB3", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "PSG3", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "UBE2L5", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "DMP1", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "ADAD1", "external_transcript_name"] <- "ADAD1-204"
CT_list[CT_list$external_gene_name == "IL5", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "ZNF679", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "RNF148", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "OR8G5", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "UBE2L5", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "VRTN", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "MBD3L5", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "SPANXN1", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "ROR1-AS1", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC01344", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "SMYD3-AS1", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "MYCNOS", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC01913", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC01087", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "KIAA2012-AS1", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC01986", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC01980", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC01811", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "IGF2BP2-AS1", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC02020", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LEF1-AS1", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "SNHG27", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC02492", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC02434", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC02241", "external_transcript_name"] <- "LINC02241-226"
CT_list[CT_list$external_gene_name == "LINC02113", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC02522", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC01016", "external_transcript_name"] <- "LINC01016-206"
CT_list[CT_list$external_gene_name == "TDRG1", "external_transcript_name"] <- "TDRG1-201"
CT_list[CT_list$external_gene_name == "SIM1-AS1", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC00326", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC01010", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "ST7-OT4", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC02669", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC01517", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC00867", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC02742", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "DLG2-AS2", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC02367", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC02378", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "MRPS35-DT", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC02156", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC01198", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC00392", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC02293", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "G2E3-AS1", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC00221", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC01193", "external_transcript_name"] <- "LINC01193-201"
CT_list[CT_list$external_gene_name == "EWSAT1", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "INSYN1-AS1", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC02152", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC02184", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC02076", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC02074", "external_transcript_name"] <- "LINC02074-201"
CT_list[CT_list$external_gene_name == "LINC01925", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC01895", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC01029", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC00664", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "MIR663AHG", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "ZBTB46-AS1", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC01203", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "MAGEA4-AS1", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "FAM197Y7", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "FAM197Y6", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "UMODL1", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "LINC00112", "external_transcript_name"] <- NA

# Genes that are not in the list anymore but were verified at one point
# CT_list[CT_list$external_gene_name == "LKAAEAR1", "external_transcript_name"] <- NA
# CT_list[CT_list$external_gene_name == "RTL9", "external_transcript_name"] <- "RTL9-202"
# CT_list[CT_list$external_gene_name == "LALBA", "external_transcript_name"] <- NA
# CT_list[CT_list$external_gene_name == "NEUROG2-AS1", "external_transcript_name"] <- NA
# CT_list[CT_list$external_gene_name == "MIR3681HG", "external_transcript_name"] <- NA
# CT_list[CT_list$external_gene_name == "LINC02912", "external_transcript_name"] <- NA
# CT_list[CT_list$external_gene_name == "LINC00452", "external_transcript_name"] <- NA
# CT_list[CT_list$external_gene_name == "LINC01531", "external_transcript_name"] <- "LINC01531-202"

# Testis transcript is not the one found in tumors
CT_list[CT_list$external_gene_name == "SUN3", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "SLC7A11-AS1", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "HORMAD1", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "P2RX3", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "ARHGAP29-AS1", "external_transcript_name"] <- NA
CT_list[CT_list$external_gene_name == "TUBB8B", "external_transcript_name"] <- NA


## Change transcripts infos according to selected transcript
CT_list <- CT_list %>%
  dplyr::select(-c(ensembl_transcript_id, chromosome_name, strand, transcript_start,
                   transcription_start_site, transcript_end, transcript_length,
                   transcript_biotype)) %>%
  left_join(transcripts_infos %>%
              dplyr::select(ensembl_gene_id,
                            external_transcript_name, ensembl_transcript_id,
                            chromosome_name, strand, transcription_start_site,
                            transcript_length, transcript_biotype))
## Remove genes with TSS set to NA
CT_list <- CT_list %>%
  filter(!is.na(transcription_start_site))

################################################################################
## Add gene family column
################################################################################

CT_list$family <- NA
CT_list[grep(pattern = "MAGE", x = CT_list$external_gene_name,
             value = FALSE), "family"] <- "MAGE"
CT_list[grep(pattern = "CT45A", x = CT_list$external_gene_name,
             value = FALSE), "family"] <- "CT45A"
CT_list[grep(pattern = "CTAG", x = CT_list$external_gene_name,
             value = FALSE), "family"] <- "CTAG"
CT_list[grep(pattern = "GAGE", x = CT_list$external_gene_name,
             value = FALSE), "family"] <- "GAGE"
CT_list[grep(pattern = "XAGE", x = CT_list$external_gene_name,
             value = FALSE), "family"] <- "GAGE"
CT_list[grep(pattern = "PAGE", x = CT_list$external_gene_name,
             value = FALSE), "family"] <- "PAGE"
CT_list[grep(pattern = "SPANX", x = CT_list$external_gene_name,
             value = FALSE), "family"] <- "SPANX"
CT_list[grep(pattern = "SSX", x = CT_list$external_gene_name,
             value = FALSE), "family"] <- "SSX"
CT_list[grep(pattern = "VCX", x = CT_list$external_gene_name,
             value = FALSE), "family"] <- "VCX"
CT_list[grep(pattern = "TSPY", x = CT_list$external_gene_name,
             value = FALSE), "family"] <- "TSPY"
CT_list[grep(pattern = "RBMY", x = CT_list$external_gene_name,
             value = FALSE), "family"] <- "RBMY"

save(CT_list, file = "../extdata/CT_list.rda",
     compress = "xz",
     compression_level = 9)


