#' Genes expression in GTEX
#'
#' @description Gene expression data in normal tissues from GTEX database.
#'
#' @format A `SummarizedExperiment` object with 24359 rows and 30 columns
#' - Rows correspond to genes (ensembl_gene_id)
#' - Columns correspond to tissues
#' - Expression data from the assay are TPM values
#'
#' @details The RowData contains
#' - A column named `GTEX_category`, specifying the tissue specificity
#' category ("testis_specific", "testis-preferential", "lowly_expressed" or
#' "other") assigned to each gene using expression values in testis and in
#' somatic tissues, has been added to the rowData.
#' - A column named `TPM_testis` giving the expression level in testis.
#' - A column named `max_TPM_somatic` giving the maximum expression level
#' found in a somatic tissue.
#' - A column named `q75_TPM_somatic` giving the quantile 75% expression level
#' found in a somatic tissue.
#'
#' @source Downloaded from
#' https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz.
#' Some categories of tissues were pooled (mean expression values are
#' given in pooled tissues). (see inst/scripts/make_GTEX_data.R for details)
"GTEX_data"

#' Genes expression data in CCLE
#'
#' @description Gene expression data in cancer cell lines from CCLE
#'
#' @format A SummarizedExperiment object with 24327 rows and 1229 columns
#' - Rows correspond to genes (ensembl_gene_id)
#' - Columns correspond to CCLE cell lines
#' - Expression data from the assay are TPM values
#' - Cell lines metadata are stored in colData
#'
#' @details The RowData contains
#' - a column `percent_of_positive_CCLE_cell_lines` that gives the
#' percent of CCLE cell lines (all cell lines combined) expressing
#' the gene at a highly level (TPM >= 10).
#' - a column `percent_of_negative_CCLE_cell_lines` that gives the percent
#' of CCLE cell lines (all cell lines combined) in which genes are
#' completely repressed (TPM < 0.1)
#' - a column `CCLE_category` gives the category ("activated",
#' "not_activated", "leaky") assigned to each gene.
#' "activated" category corresponds to genes expressed in at least
#' one cell line and repressed in at least 20% of cell lines.
#'
#' @source TPM values downloaded using depmap bioconductor package.
#' (see inst/scripts/make_CCLE_data.R for details)
"CCLE_data"

#' Gene expression values in normal tissues with or without multimapping
#'
#' @description Gene expression values (TPM) in a set of normal tissues
#' obtained by counting or not multi-mapped reads. Many CT genes belong to gene
#' families from which members have identical or nearly identical sequences.
#' Some CT can only be detected in RNAseq data in which multimapping reads are
#' not discared.
#'
#' @format A SummarizedExperiment object with 24359 rows and 18 columns
#' - Rows correspond to genes (ensembl_gene_id)
#' - Columns correspond to normal tissues
#' - `TPM_no_multimapping` assay gives TPM expression values obtained
#' when discarding multimapped reads
#' - `TPM_with_multimapping` assay gives TPM expression values obtained
#' by counting multimapped reads
#'
#' @details A column named `multimapping_analysis` has been added to the
#' rowData. It summarizes the testis specificity analysis of genes flagged
#' as "lowly_expressed" in GTEX_data.
#'
#' @source RNAseq fastq files were downloaded from Encode database.
#' (see inst/scripts/make_normal_tissues_multimapping.R for details)
"normal_tissues_multimapping_data"

#' RNseq expression values and differential expression analysis of
#' genes in cell lines treated or not with a demethylating agent.
#'
#' @description Gene expression values in a set of cell lines treated
#' or not with 5-Aza-2′-Deoxycytidine (DAC), a demethylating agent.
#'
#' @format A SummarizedExperiment object with 24359 rows and 32 columns
#' - Rows correspond to genes (ensembl reference)
#' - Columns correspond to samples
#' - Expression data are normalised counts log-transformed (log1p)
#' - The colData contains the SRA references from which fastq files
#' were downloaded
#' - Results of a differential expression evaluating the DAC treatment
#' effect are stored in the rowData.
#'
#' @details Differential expression analysis was done using DESeq2_1.36.0,
#' using as design = ~ cell + treatment
#' (see inst/scripts/make_DAC_treated_cells.R for details).
#'
#' @source RNAseq fastq files were downloaded from Encode database.
#' SRA reference of samples are stored in the colData.
"DAC_treated_cells"

#' RNseq expression values (multi-mapping was allowed) and differential
#' expression analysis of genes in cell lines treated or not with a
#' demethylating agent.
#'
#' @description Gene expression values in a set of cell lines treated
#' or not with 5-Aza-2′-Deoxycytidine (DAC), a demethylating agent.
#' Many CT genes belong to gene families from which members have
#' identical or nearly identical sequences. Some CT can only be detected
#' in RNAseq data in which multimapping reads are not discared.
#'
#' @format A SummarizedExperiment object with 24359 rows and 32 columns
#' - Rows correspond to genes (ensembl_gene_id)
#' - Columns correspond to samples
#' - Expression data are normalised counts log-transformed (log1p)
#' - The colData contains the SRA references from which fastq files
#' were downloaded
#' - Results of a differential expression evaluating the DAC treatment
#' effect are stored in the rowData.
#'
#' @details Differential expression analysis was done using DESeq2_1.36.0,
#' using as design = ~ cell + treatment
#' (see inst/scripts/make_DAC_treated_cells_multimapping.R for details).
#'
#' @source RNAseq fastq files were downloaded from Encode database.
#' SRA reference of samples are stored in the colData.
"DAC_treated_cells_multimapping"

#' Promoter methylation of Cancer-Testis Genes in normal tissues
#'
#' @description Methylation values of CpG located within Cancer-Testis
#' (CT) promoters in a set of normal tissues
#'
#' @format A RangedSummarizedExperiment object with 53578 rows and 14 columns
#' - Rows correspond to CpG positions (located within CT genes promoters)
#' - Columns correspond to normal tissues
#' - Methylation values from WGBS data
#'
#' @source WGBS methylation data was downloaded from Encode and from GEO databases.
#' (see inst/scripts/make_CT_methylation_in_tissues.R for details)
"CT_methylation_in_tissues"

#' Promoter mean methylation of Cancer Testis Genes in normal tissues
#'
#' @description Mean methylation values of CpG located within Cancer-Testis
#' (CT) promoters in a set of normal tissues
#'
#' @format A SummarizedExperiment object with 306 rows and 14 columns
#' - Rows correspond to CT genes
#' - Mean methylation levels in normal tissues are stored in columns
#' - CpG densities and results of methylation analysis are stored
#' in rowData
#'
#' @details The RowData contains:
#' - A column named `CpG_density`, gives the density of CpG within each promoter
#' (number of CpG / promoter length * 100).
#' - A column `CpG_promoter` that classifies the promoters according to their
#' CpG densities: "CpG_low" (CpG_density < 2), "CpG_intermediate"
#' (CpG_density >= 2 & CpG_density < 4), and "CpG_high (CpG_density >= 4).
#' - A column `methylation_in_tissues` summarising the methylation analysis
#' in normal tissues.
#'
#' @source WGBS methylation data was downloaded from Encode and from GEO
#' databases. Mean methylation levels are evaluated using methylation values
#' of CpG located in promoter region (defined as 1000 nt upstream TSS and
#' 200 nt downstream TSS) (see inst/scripts/make_CT_mean_methylation_in_tissues.R
#' for details)
"CT_mean_methylation_in_tissues"

#' Gene expression in TCGA samples
#'
#' @description Gene expression data in TCGA samples (tumor and peritumoral samples)
#'
#' @format A SummarizedExperiment object with 24350 rows and 4087 columns
#' - Rows correspond to CT genes (ensembl reference)
#' - Columns correspond to samples
#' - Expression data from the assay are TPM values
#' - Clinical information are stored in colData
#' - Genes information are stored in RowData
#'
#' @source SKCM, LUAD, LUSC, COAD, ESCA, BRCA and HNSC expression data were downloaded
#' with TCGAbiolinks. Global hypomethylation levels from paper `DNA methylation
#' loss promotes immune evasion of tumours with high mutation and copy number load`
#' from Jang et al., Nature Commun 2019 were added to colData.
#' (see inst/scripts/make_TCGA_TPM.R for details)
"TCGA_TPM"

#' Methylation of CT promoters in TCGA samples
#'
#' @description Methylation values of probes located within Cancer-Testis
#' (CT) promoters in samples from TCGA (tumor and peritumoral samples)
#'
#' @format A RangedSummarizedExperiment object with 689 rows and 3423 columns
#' - Rows correspond to Infinium 450k probes
#' - Columns correspond to samples
#' - Methylation data from the assay are Beta values
#' - Clinical information are stored in colData
#' - probe information (hg38 coordinates) are stored in RowData
#'
#' @source SKCM, LUAD, LUSC, COAD, ESCA, BRCA and HNSC methylation data were
#' downloaded with TCGAbiolinks and subsetted to select probes located in CT genes
#' promoter regions (see inst/scripts/make_TCGA_CT_methylation.R for details)
"TCGA_CT_methylation"

#' CT genes description table
#'
#' @description Cancer-Testis (CT) genes description
#'
#' @format A `tibble` object with 306 rows and 33 columns.
#' - Rows correspond to CT genes
#' - Columns give CT genes characteristics
#'
#' @details CT_genes characteristics column:
#' - Column `family` gives the gene family name.
#' - Column `X_linked` classifies genes in "X-linked" or "not_X_linked"
#' - Column `GTEX_category` gives the category ("testis_specific",
#' "testis_preferential" or "lowly_expressed") assigned to each gene
#' using GTEx database
#' - Column `TPM_testis` gives the gene expression level in testis
#' (using GTEX database)
#' - Column `max_TPM_somatic` gives the maximum expression level
#' found in a somatic tissue (using GTEX database)
#' - Column `q75_TPM_somatic` gives the quantile 75% expression level
#' found in a somatic tissue (using GTEX database)
#' - Column `multimapping_analysis` informs if the gene was found to
#' be testis-specific when multi-mapped reads were counted for gene
#' expression in normal tissues
#' - Column `testis_specificity` gives the testis-specificity of genes
#' assigned to each gene using `GTEX_category` and `multimapping_analysis`.
#' - Column `percent_of_positive_CCLE_cell_lines` gives the percent of
#' CCLE cancer cell lines in which genes are expressed (genes were
#' considered as expressed if TPM >= 10)
#' - Column `percent_of_negative_CCLE_cell_lines` gives the percent of
#' CCLE cancer cell lines in which genes are repressed (TPM <= 0.1)
#' - Column `max_TPM_in_CCLE` gives the highest expression level of genes
#' in CCLE cell lines.
#' - Column `CCLE_category` gives the category assigned to each gene
#' using CCLE data. "Activated" category corresponds to genes expressed
#' in at least one cell line (TPM >= 10) and repressed in at least 20% of cell lines
#' - Column `percent_pos_tum` gives the percent of
#' TCGA cancer samples in which genes are expressed (genes were
#' considered as expressed if TPM >= 10)
#' - Column `percent_neg_tum` gives the percent of
#' TCGA cancer samples in which genes are repressed (TPM <= 0.1)
#' - Column `max_TPM_in_TCGA` gives the highest expression level of genes
#' in TCGA cancer sample
#' - Column `TCGA_category` gives the category assigned to each gene
#' using TCGA data. "Activated" category corresponds to genes expressed
#' in at least one tumor (TPM >= 10) and repressed in at least 20% of samples
#' - Column `DAC` summarises the results ("induced" or "not_induced") of a
#' differential expression evaluating gene induction upon DAC treatment
#' in a series of cell lines
#' - Column `methylation_in_tissues` summarises the analysis of gene
#' promoter methylation in somatic and germline normal tissues
#' ("methylated_in_somatic_unmethylated_in_germline", "methylated_in_somatic_and_germline"
#' or "unmethylated_in_somatic")
#' - Column `regulation` summarises the regulation category ("methylation"
#' or "not_methylation") that was assigned to genes based on DAC induction
#' and on promoter methylation levels in normal tissues (when available)
#' - Column `met_exp_corr_TCGA` gives the coefficient correlation between promoter
#' methylation and expression in TCGA samples
#' - Column named `CpG_density`, gives the density of CpG within each promoter
#' (number of CpG / promoter length * 100)
#' - Column `CpG_promoter` classifies the promoters according to their
#' CpG densities: "CpG_low" (CpG_density < 2), "CpG_intermediate"
#' (CpG_density >= 2 & CpG_density < 4), and "CpG_high (CpG_density >= 4)
#' - Columns `external_transcript_name`, `ensembl_transcript_id`, `chromosome_name`,
#' `strand`, `transcription_start_site`, `transcript_length` and `transcript_biotype`
#' give the references and informations about the most biologically relevant transcript
#' associated to each gene.
#' - Columns `oncogene` and `tumor_suppressor` informs if oncogenic and tumor-suppressor
#' functions have been associated to genes (source: [Cancermine](http://bionlp.bcgsc.ca/cancermine/))
#'
#' @source (see inst/scripts/make_CT_genes.R for details)
"CT_genes"

#' CCLE correlation_matrix
#'
#' @description CCLE correlation matrix
#'
#' @format A `matrix` object with 306 rows and 24327 columns.
#' - Rows correspond to CT genes
#' - Rows correspond to all genes from CCLE database
#'
#' @details Correlation coefficients between CT genes and all other genes
#' are given in the matrix. These correlation coefficients were calculated
#' from using log transformed expression values from `CCLE_data`.
#'
#' @source (see inst/scripts/make_CCLE_correlation_matrix.R for details)
"CCLE_correlation_matrix"


