#' Genes expression in GTEX
#'
#' @description
#'
#' Gene expression data in normal tissues from GTEx database.
#'
#' @format
#'
#' A `SummarizedExperiment` object with 24504 rows and 32 columns
#'
#' - Rows correspond to genes (ensembl_gene_id as rownames)
#' - Columns correspond to tissues
#' - Expression data from the assay are TPM values
#'
#' @details
#'
#' The rowData contains
#'
#' - A column named `GTEX_category`, specifying the tissue specificity
#'   category ("testis_specific", "testis-preferential",
#'   "lowly_expressed" or "other") assigned to each gene using
#'   expression values in testis and in somatic tissues, has been
#'   added to the rowData. "testis_specific" genes are expressed
#'   exclusively in testis (expression in testis >= 1 TPM, highest
#'   expression in somatic tissues < 0.5 TPM, and expressed at least
#'   10x more in testis than in any somatic
#'   tissue). "testis-preferential" genes are genes expressed in
#'   testis but also in a few somatic tissues (expression in testis >=
#'   1 TPM, and allowed in a minority of somatic tissues
#'   (q75_TPM_somatic < 0.5) and expressed at least 10x more in testis than in
#'   any somatic tissue). "lowly_expressed" genes are genes undetectable in GTEX
#'   database probably due to multi-mapping issues (expression in all
#'   GTEX tissues < 1 TPM).
#'
#' - A column named `q75_TPM_somatic` giving the quantile 75% of TPM
#'   in a somatic tissue.
#' - A column named `max_TPM_somatic` giving the maximum expression
#'   level found in a somatic tissue.
#' - A column named `ratio_testis_somatic` giving the ratio between the TPM
#'   in testis and the max TPM in a somatic tissue
#'
#' @source
#'
#' Downloaded from
#' https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz.
#' Some categories of tissues were pooled (mean expression values are
#' given in pooled tissues) (see `scripts/make_GTEX_data.R` for
#' details).
#'
#' @name GTEX_data
#'
#' @docType data
NULL

#' Testis scRNAseq data
#'
#' @description
#'
#' Testis single cell RNAseq data from
#' `The adult human testis transcriptional cell atlas` (Guo et al. 2018)
#'
#' @format
#'
#' A `SingleCellExperiment` object with 20891 rows and 6490 columns
#'
#' - Rows correspond to genes (gene names as rownames)
#' - Columns correspond to testis cells
#'
#' @details
#'
#' Description of the colData:
#'
#' - Column `nGene` gives the number of distinct genes detected per cell.
#'
#' - Column `nUMI` gives the total UMI number per cell.
#'
#' - Column `clusters` gives cluster number defined in the Guo's paper.
#'
#' - Column `type` gives the testis cell type associated to the cluster number.
#'
#' - Column `Donor` gives the Donor origin of the cells.
#'
#' Description of the rowData:
#'
#' - Column `percent_pos_testis_germcells` gives the percent of testis germ cells
#'   in wich the genes are detected (count > 0) (based on testis scRNAseq data).
#'
#' - Column `percent_pos_testis_somatic` gives the percent of testis somatic cells
#'   in wich the genes are detected (count > 0) (based on testis scRNAseq data).
#'
#' - Column`testis_cell_type` specifies the testis cell-type showing the highest
#'  mean expression of each gene (based on testis scRNAseq data).
#'
#' The rowData contains the `testis_cell_type` column, specifying the testis
#' cell-type showing the highest mean expression of each gene.
#'
#' @source
#'
#' The count matrix `GSE112013_Combined_UMI_table.txt.gz` was downloaded from
#' GEO (accession: GSE112013). Metadata correspond to `TableS1` from the paper's
#' supplemental data. The data were converted in a SingleCellExperiment
#' (see `scripts/13_make_testis_sce.R` for details).
#'
#' @name testis_sce
#'
#' @docType data
NULL

#' Gene expression in human cell types
#'
#' @description
#'
#' Gene expression profiles in different human cell types based on scRNAseq data
#' obtained from the Human Protein Atlas (https://www.proteinatlas.org)
#'
#' @format
#'
#' A `SingleCellExperiment` object with 20082 rows and 66 columns
#'
#' - Rows correspond to genes (ensembl gene id as rownames)
#' - Columns correspond to cell types
#' - Expression values correspond to transcripts per million protein coding genes
#' (pTPM)
#'
#' @details
#'
#' Description of the colData:
#'
#' - Column `Cell_type` gives cell type.
#'
#' - Column `group` gives the cell type group (defined in the Human Protein Atlas).
#'
#' Description of the rowData:
#'
#' - Column `max_TPM_in_a_somatic_cell_type` gives the maximum expression value
#' found in a somatic cell type
#'
#' - Column `max_in_germcells_group` gives the maximum expression value found
#' in a germ cell type
#'
#' - Column `Higher_in_somatic_cell_type` specifies if a somatic cell type
## was found to express the gene at a higher level than any germ cell type
#'
#' @source
#'
#' Gene expression values in cell types, based on multiple scRNAseq datasets
#' obtained from the Human Protein Atlas
#' (https://www.proteinatlas.org/about/download)
#' The data were converted in a SummarizedExperiment
#' (see `scripts/14_make_scRNAseq_HPA.R` for details).
#'
#' @name scRNAseq_HPA
#'
#' @docType data
NULL

#' Genes expression data in CCLE
#'
#' @description
#'
#' Gene expression data in cancer cell lines from CCLE
#'
#' @format
#'
#' A `SummarizedExperiment` object with 24473 rows and 1229 columns
#'
#' - Rows correspond to genes (ensembl_gene_id)
#' - Columns correspond to CCLE cell lines
#' - Expression data from the assay are TPM values
#' - Cell lines metadata are stored in colData
#'
#' @details
#'
#' The rowData contains
#'
#' - A column `percent_of_positive_CCLE_cell_lines` that gives the
#'   percentage of CCLE cell lines (all cell lines combined)
#'   expressing the gene (TPM >= 1).
#'
#' - A column `percent_of_negative_CCLE_cell_lines` that gives the
#'   percent of CCLE cell lines (all cell lines combined) in which
#'   genes are repressed (TPM < 0.5)
#'
#' - A column `max_TPM_in_CCLE` that gives the maximal expression (in
#'   TPM) found in all cell lines.
#'
#' - A column `CCLE_category` gives the category ("activated",
#'   "not_activated", "leaky") assigned to each gene.  "activated"
#'   category corresponds to genes expressed (TPM >= 1) in at
#'   least 1% of cell lines, repressed (TPM <= 0.5) in at least 20% of
#'   cell lines with a maximal expression higher than 5 TPM.
#'   "not_activated" category corresponds to genes
#'   repressed (TPM <= 0.5) in at least 20% of cell lines but
#'   expressed (TPM >= 1) less than 1%.  "leaky" category
#'   corresponds to genes repressed (TPM <= 0.5) in less than 20% of
#'   cell lines. "lowly_expressed" corresponds to genes repressed (TPM <= 0.5)
#'   in at least 20%, expressed (TPM >= 1) in more than 1 % of cell lines, with
#'   a maximum expression lower than 5 TPM.
#'
#'
#' @source
#'
#' TPM values downloaded using depmap bioconductor package (see
#' `scripts/make_CCLE_data.R` for details).
#'
#' @name CCLE_data
#'
#' @docType data
NULL

#' Gene expression values in normal tissues
#'
#' @description
#'
#' Gene expression values (TPM) in a set of normal tissues obtained by
#' counting or not multi-mapped reads. Many CT genes belong to gene
#' families from which members have identical or nearly identical
#' sequences.  Some CT can only be detected in RNAseq data in which
#' multimapping reads are not discarded.
#'
#' @format
#'
#' A `SummarizedExperiment` object with 24504 rows and 18 columns
#'
#' - Rows correspond to genes (ensembl_gene_id)
#'
#' - Columns correspond to normal tissues.
#'
#' - First assay, `TPM_no_multimapping`, gives TPM expression values
#'   obtained when discarding multimapped reads.
#'
#' - Second assay, `TPM_with_multimapping`, gives TPM expression
#'   values obtained by counting multimapped reads.
#'
#' @details
#'
#' A column named `multimapping_analysis` has been added to the
#' rowData. It summarizes the testis specificity analysis of genes
#' flagged as "lowly_expressed" in GTEX_data. Genes are considered
#' "testis_specific" when, with multimapping allowed, they are
#' detectable in testis (TPM >= 1), their TPM value has increased
#' compared to without multimapping (ratio > 5), and their TPM value
#' is at least 10 times higher in testis than in any other somatic
#' tissue (where the maximum expression always has to be below 1 TPM).
#' Genes are considered "testis_preferential" when, with multimapping allowed,
#' they are detectable in testis (TPM >= 1), their TPM value has increased
#' compared to without multimapping (ratio > 5), and their TPM value
#' is at least 10 times higher in testis than in any other somatic
#' tissue (where the maximum expression is above 1 TPM).
#'
#' @source
#'
#' RNAseq fastq files were downloaded from Encode database (see
#' `scripts/make_normal_tissues_multimapping.R` for details).
#'
#' @name normal_tissues_multimapping_data
#'
#' @docType data
NULL


#' DE genes with/without demethylating agent
#'
#' @description
#'
#' Gene expression values in a set of cell lines treated or not with
#' 5-Aza-2'-Deoxycytidine (DAC), a demethylating agent.
#'
#' @format
#'
#' A `SummarizedExperiment` object with 24516 rows and 32 columns
#'
#' - Rows correspond to genes (ensembl_gene_id).
#'
#' - Columns correspond to samples.
#'
#' - Expression data correspond to counts that have been normalised
#'   (by DESeq2 method) and log-transformed (log1p).
#'
#' - The colData contains the SRA references of the fastq files that
#'   were downloaded, and informations about the cell lines and the
#'   DAC treatment.
#'
#' - The rowData contains the results of a differential expression
#'   evaluating the DAC treatment effect. For each each cell line, the
#'   log2FC between treated and control cells is given, as well as the
#'   p-adjusted value. The column `induced` flags genes significantly
#'   induced (log2FoldChange >= 2 and padj <= 0.1) in at least one
#'   cell line. The threshold is not too stringent as DAC is expected to induce
#'   low expression levels (demethylation doesn't necessarily occurs in
#'   all treated cells...).
#'   When all cells lines already express the gene before DAC treatment, no
#'   assessment of induction was done.
#'
#' @details
#'
#' Differential expression analysis was done using DESeq2_1.36.0,
#' using as `design = ~ treatment` (see
#' `scripts/make_DAC_treated_cells.R` for details).
#'
#' @source RNAseq
#'
#' fastq files were downloaded from Encode database.  SRA reference of
#' samples are stored in the colData.
#'
#' @name DAC_treated_cells
#'
#' @docType data
NULL

#' DE genes treated or not with a demethylating agent
#'
#' @description
#'
#' Gene expression values in a set of cell lines treated or not with
#' 5-Aza-2'-Deoxycytidine (DAC), a demethylating agent.  Many CT genes
#' belong to gene families from which members have identical or nearly
#' identical sequences. Some CT can only be detected in RNAseq data in
#' which multimapping reads are not discarded.
#'
#' @format
#'
#' A `SummarizedExperiment` object with 24516 rows and 32 columns
#'
#' - Rows correspond to genes (ensembl_gene_id).
#'
#' - Columns correspond to samples.
#'
#' - Expression data correspond to counts that have been normalised
#'   (by DESeq2 method) and log-transformed (log1p).
#'
#' - The colData contains the SRA references of the fastq files that
#'   were downloaded, and informations about the cell lines and the
#'   DAC treatment.
#'
#' - The rowData contains the results of a differential expression
#'   evaluating the DAC treatment effect. For each each cell line, the
#'   log2FC between treated and control cells is given, as well as the
#'   p-adjusted value. The column `induced` flags genes significantly
#'   induced (log2FoldChange >= 2 and padj <= 0.1) in at least one
#'   cell line. The threshold is not too stringent as DAC is expected to induce
#'   low expression levels (demethylation doesn't necessarily occurs in
#'   all treated cells...).
#'   When all cells lines already express the gene before DAC treatment, no
#'   assessment of induction was done.
#'
#' @details
#'
#' Differential expression analysis was done using DESeq2_1.36.0,
#' using as `design = ~ treatment` (see
#' `scripts/make_DAC_treated_cells_multimapping.R` for details).
#'
#' @source
#'
#' RNAseq fastq files were downloaded from Encode database.  SRA
#' reference of samples are stored in the colData.
#'
#' @name DAC_treated_cells_multimapping
#'
#' @docType data
NULL

#' Methylation of CpGs within all genes promoters
#'
#' @description
#'
#' Methylation values of CpGs located within all genes
#' promoters in a set of normal tissues.
#'
#' @format
#'
#' A `RangedSummarizedExperiment` object with 4280327 rows and 14 columns
#'
#' - Rows correspond to CpGs (located within all genes promoters)
#'
#' - Columns correspond to normal tissues
#'
#' - Methylation values from WGBS data
#'
#' - rowRanges correspond to CpG positions
#'
#' @source
#'
#' WGBS methylation data was downloaded from Encode and from GEO
#' databases (see `scripts/make_methylation_in_tissues.R` for
#' details).
#'
#' @name methylation_in_tissues
#'
#' @docType data
NULL


#' All genes' promoters mean methylation
#'
#' @description
#'
#' Mean methylation values of all CpGs located within all genes
#' promoters in a set of normal tissues
#'
#' @format
#'
#' A `SummarizedExperiment` object with 24502 rows and 14 columns
#'
#' - Rows correspond to all genes (gene names as rownames)
#'
#' - Mean methylation levels in normal tissues are stored in columns
#'
#' - CpG densities and results of methylation analysis are stored in
#'   rowData
#'
#' @details
#'
#' The rowData contains:
#'
#' - A column named `CpG_density`, gives the density of CpG within
#'   each promoter (number of CpG / promoter length * 100).
#'
#' - A column `CpG_promoter` that classifies the promoters according
#'   to their CpG densities: "low" (CpG_density < 2), "intermediate"
#'   (CpG_density >= 2 & CpG_density < 4), and "high" (CpG_density >=
#'   4).
#'
#' - A column `somatic_met_level` that gives the mean methylation
#'   level of each promoter in somatic tissues.
#'
#' - A column `sperm_met_level` that gives the methylation level of
#'   each promoter in sperm.
#'
#' - A column `somatic_methylation` indicates if the promoter's mean
#'   methylation level in somatic tissues is higher than 50%.
#'
#' - A column `germline_methylation`indicates if the promoter is
#'   methylated in germline, based on the ratio with somatic tissues
#'   (FALSE if somatic_met_level is at least twice higher than
#'   germline_met_level).
#'
#' @source
#'
#' WGBS methylation data was downloaded from Encode and from GEO
#' databases. Mean methylation levels are evaluated using methylation
#' values of CpGs located in promoter region (defined as 1000 nt
#' upstream TSS and 200 nt downstream TSS) (see
#' `scripts/make_mean_methylation_in_tissues.R` for details).
#'
#' @name mean_methylation_in_tissues
#'
#' @docType data
NULL

#' Gene expression in TCGA samples
#'
#' @description
#'
#' Gene expression data in TCGA samples (tumor and peritumoral
#' samples).
#'
#' @format
#'
#' A `SummarizedExperiment` object with 24497 rows and 4141 columns
#'
#' - Rows correspond to genes (ensembl_gene_id)
#'
#' - Columns correspond to samples
#'
#' - Expression data from the assay are TPM values
#'
#' - Clinical information are stored in colData
#'
#' - Genes information are stored in rowData
#'
#' @details
#'
#' - The colData contains clinical data from TCGA as well as global
#'   hypomethylation levels obtained from paper *DNA methylation loss
#'   promotes immune evasion of tumours with high mutation and copy
#'   number load* from Jang et al., Nature Commun 2019 that were added
#'   (see `inst/scripts/make_TCGA_TPM.R` for details).
#'
#' - The rowData contains genes information and, for each gene, the
#'   percentage of tumors that are positive (TPM >= 1), and the
#'   percentage of tumors that are negative (TPM < 0.5). In column
#'   `TCGA_category`, genes are labelled as "activated" when the
#'   percentage of positive tumors is > 1, with a maximal expression higher than
#'   5 TPM, and when at least 20% of tumors are negative. Genes are labelled as
#'   "not_activated" when the percentage of positive tumors is lower than 1.
#'   Genes are labelled as "leaky" when less than 20% of tumors are negative.
#'   Genes are labelled as "lowly_expressed" when repressed (TPM <= 0.5)
#'   in at least 20%, expressed (TPM >= 1) in more than 1 % of cell lines, with
#'   a maximum expression lower than 5 TPM.
#'
#' @source
#'
#' SKCM, LUAD, LUSC, COAD, ESCA, BRCA and HNSC expression data were
#' downloaded with TCGAbiolinks (see `scripts/make_TCGA_TPM.R`
#' for details).
#'
#' @name TCGA_TPM
#'
#' @docType data
NULL

#' Methylation of all genes promoters in TCGA samples
#'
#' @description
#'
#' Methylation values of probes located within all genes
#' promoters in samples from TCGA (tumor and peritumoral samples)
#'
#' @format
#'
#' A `RangedSummarizedExperiment` object with 79445 rows and 3423
#' columns
#'
#' - Rows correspond to Infinium 450k probes
#'
#' - Columns correspond to samples
#'
#' - Methylation data from the assay are Beta values
#'
#' - Clinical information are stored in colData
#'
#' - Probe information (hg38 coordinates) are stored in rowRanges
#'
#' @source
#'
#' SKCM, LUAD, LUSC, COAD, ESCA, BRCA and HNSC methylation data were
#' downloaded with TCGAbiolinks and subsetted to select probes located
#' in CT genes promoter regions (see
#' `scripts/make_TCGA_methylation.R` for details).
#'
#' @name TCGA_methylation
#'
#' @docType data
NULL

#' CT genes description table
#'
#' @description
#'
#' Cancer-Testis (CT) genes description
#'
#' @format
#'
#' A `tibble` object with 280 rows and 47 columns.
#'
#' - Rows correspond to CT genes
#'
#' - Columns give CT genes characteristics
#'
#' @details
#'
#' When the promoter is mentionned, it has been determined as 1000 nt
#' upstream TSS and 200 nt downstream TSS.
#'
#' CT_genes characteristics column:
#'
#' - Column `CT_gene_type` indicates if the gene is a CT specific gene
#'  ("CT_gene" : testis_specific in `testis_specificity`) and activated in
#'  "TCGA_category" and "CCLE_category) or CT preferential gene ("CTP_gene" :
#'  testis_preferential in `testis_specificity`) and activated in"TCGA_category"
#'   and "CCLE_category").
#'
#' - Column `testis_specificity` gives the testis-specificity of genes
#'   assigned to each gene using `GTEX_category` and `multimapping_analysis`
#'   ("testis_specific" or "testis_preferential"). Genes were assigned
#'   "testis-preferential" if testis-specific in these categories but not testis
#'   specific in `HPA_category` or leaky in `CCLE_category` or `TCGA_category`.
#'
#' - Column `regulated_by_methylation` indicates if the gene is
#'   regulated by methylation (`TRUE`) based on DAC induction (has to
#'   be TRUE) and on promoter methylation level in normal somatic tissues
#'   (when available, has to be methylated in somatic tissues).
#'
#' - Column `X_linked` indicates if the gene is on the chromosome X
#'   (TRUE) or not (FALSE).
#'
#' - Columns `chr`, `strand` and `transcription_start_site` give the
#'   genomic location.
#'
#' - Column `GTEX_category` gives the category ("testis_specific",
#'   "testis_preferential" or "lowly_expressed") assigned to each gene
#'   using GTEx database (see `?GTEX_data` for details).
#'
#'  - Column `q75_TPM_somatic` gives the q75 expression level found
#'   in a somatic tissue (using GTEx database).
#'
#' - Column `max_TPM_somatic` gives the maximum expression level found
#'   in a somatic tissue (using GTEx database).
#'
#' - Column `ratio_testis_somatic` gives the ratio between expression in testis
#'   and the highest expression found in a somatic tissue (using GTEx database).
#'
#' - Column `TPM_testis` gives the gene expression level in testis
#'   (using GTEx database).
#'
#' - Column `lowly_expressed_in_GTEX` indicates if the gene is lowly
#'   expressed in GTEX database and thus needed to be analysed with
#'   multimapping allowed.
#'
#' - Column `multimapping_analysis` informs if the gene (flagged as
#'   "lowly_expressed" in GTEX_data) was found to be testis-specific
#'   when multi-mapped reads were counted for gene expression in
#'   normal tissues ("not_analysed" or "testis_specific") (see
#'   `?normal_tissues_multimapping_data` for details).
#'
#' - Column `HPA_RNA_single_cell_type_specific_nTPM` specifies the cell types in
#'   which genes were detected in the HPA single cell data (see
#'   `?HPA_cell_type_specificity` for details).
#'
#' - Column `max_HPA_germcell` specifies if the maximum expression value in a
#' germ cell type. (see  `?HPA_cell_type_specificity` for details).
#'
#' - Column `max_HPA_somatic` specifies if the maximum expression value in a
#' somatic cell type. (see  `?HPA_cell_type_specificity` for details).
#'
##' - Column `not_detected_in_somatic_HPA` specifies if the gene is detected or
#' not in a somatic cell type. (see  `?HPA_cell_type_specificity` for details).
#'
#' - Column `HPA_ratio_germ_som` gives the ratio between `max_HPA_germcell` and
#' `max_HPA_somatic` columns.

#'
#' - Column `percent_of_positive_CCLE_cell_lines` gives the percentage
#'   of CCLE cancer cell lines in which genes are expressed (genes
#'   were considered as expressed if TPM >= 1).
#'
#' - Column `percent_of_negative_CCLE_cell_lines` gives the percentage
#'   of CCLE cancer cell lines in which genes are repressed (TPM <=
#'   0.5).
#'
#' - Column `max_TPM_in_CCLE` gives the highest expression level of
#'   genes in CCLE cell lines.
#'
#' - Column `CCLE_category` gives the category assigned to each gene
#'   using CCLE data. "Activated" category corresponds to genes
#'   expressed in at least 1% of cell lines (TPM >= 1) and repressed in
#'   at least 20% of cell lines.
#'
#' - Column `percent_pos_tum` gives the percentage of TCGA cancer
#'   samples in which genes are expressed (genes were considered as
#'   expressed if TPM >= 1).
#'
#' - Column `percent_neg_tum` gives the percentage of TCGA cancer samples in
#'   which genes are repressed (TPM <= 0.5).
#'
#' - Column `max_TPM_in_TCGA` gives the highest expression level of
#'   genes in TCGA cancer sample.
#'
#' - Column `max_q75_in_NT` gives the maximum q75 expression in normal
#'   peritumoral tissues from TCGA.
#'
#' - Column `TCGA_category` gives the category assigned to each gene
#'   using TCGA data. "activated" category corresponds to genes
#'   expressed in at least 1% of tumors (TPM >= 1) and repressed in at
#'   least 20% of samples.  "multimapping_issue" corresponds to genes
#'   that need multi-mapping to be allowed in order to be analysed
#'   properly.
#'
#' - Columns `external_transcript_name`, `ensembl_transcript_id`, and
#'   `transcript_biotype` give the references and informations about
#'   the most biologically relevant transcript associated to each
#'   gene.
#'
#' - Column `IGV_backbone` indicates if a gene has been removed from CT genes
#'   as RNA-Seq reads were not properly aligned on exons, but were instead
#'   spread across a wide genomic region spanning the genes.
#'
#' - Column `family` gives the gene family name.
#'
#' - Column `DAC_induced` summarises the results (`TRUE` or `FALSE`)
#'   of a differential expression evaluating gene induction upon DAC
#'   treatment in a series of cell lines.
#'
#' - Column named `CpG_density`, gives the density of CpG within each
#'   promoter (number of CpG / promoter length * 100).
#'
#' - Column `CpG_promoter` classifies the promoters according to their
#'   CpG densities: "low" (CpG_density < 2), "intermediate"
#'   (CpG_density >= 2 & CpG_density < 4), and "high" (CpG_density >=
#'   4).
#'
#' - Column `somatic_met_level` that gives the mean methylation level
#'   of each promoter in somatic tissues.
#'
#' - Column `sperm_met_level` that gives the methylation level of each
#'   promoter in sperm.
#'
#' - Column `somatic_methylation` indicates if the promoter's mean
#'   methylation level in somatic tissues is higher than 50%.
#'
#' - Column `germline_methylation`indicates if the promoter is
#'   methylated in germline, based on the ratio with somatic tissues
#'   (`FALSE` if somatic_met_level is at least twice higher than
#'   `germline_met_level`).
#'
#' - Columns `oncogene` and `tumor_suppressor` informs if oncogenic
#'   and tumor-suppressor functions have been associated to genes
#'   (source: [Cancermine](http://bionlp.bcgsc.ca/cancermine/)).
#'
#' @source
#'
#' See `scripts/make_all_genes_prelim.R` and
#' `scripts/make_all_genes_and_CT_genes.R`for details on how this list of
#' curated CT genes was created.
#'
#' @name CT_genes
#'
#' @docType data
NULL

#' All genes genes description table
#'
#' @description
#'
#' All genes description, according to the analysis done for CT genes
#'
#' @format
#'
#' A `tibble` object with 24488 rows and 47 columns.
#'
#' - Rows correspond to genes
#'
#' - Columns give genes characteristics
#'
#' @details
#'
#' When the promoter is mentionned, it has been determined as 1000 nt
#' upstream TSS and 200 nt downstream TSS.
#'
#' CT_genes characteristics column:
#'
#' - Column `CT_gene_type` indicates if the gene is a CT specific gene
#'  ("CT_gene" : testis_specific in `testis_specificity`) and activated in
#'  "TCGA_category" and "CCLE_category) or CT preferential gene ("CTP_gene" :
#'  testis_preferential in `testis_specificity`) and activated in"TCGA_category"
#'   and "CCLE_category").
#'
#' - Column `testis_specificity` gives the testis-specificity of genes
#'   assigned to each gene using `GTEX_category` and `multimapping_analysis`
#'   ("testis_specific" or "testis_preferential"). Genes were assigned
#'   "testis-preferential" if testis-specific in these categories but not testis
#'   specific in `HPA_category` or leaky in `CCLE_category` or `TCGA_category`.
#'
#' - Column `regulated_by_methylation` indicates if the gene is
#'   regulated by methylation (`TRUE`) based on DAC induction (has to
#'   be TRUE) and on promoter methylation level in normal somatic tissues
#'   (when available, has to be methylated in somatic tissues).
#'
#' - Column `X_linked` indicates if the gene is on the chromosome X
#'   (TRUE) or not (FALSE).
#'
#' - Columns `chr`, `strand` and `transcription_start_site` give the
#'   genomic location.
#'
#' - Column `GTEX_category` gives the category ("testis_specific",
#'   "testis_preferential" or "lowly_expressed") assigned to each gene
#'   using GTEx database (see `?GTEX_data` for details).
#'
#'  - Column `q75_TPM_somatic` gives the q75 expression level found
#'   in a somatic tissue (using GTEx database).
#'
#' - Column `max_TPM_somatic` gives the maximum expression level found
#'   in a somatic tissue (using GTEx database).
#'
#' - Column `ratio_testis_somatic` gives the ratio between expression in testis
#'   and the highest expression found in a somatic tissue (using GTEx database).
#'
#' - Column `TPM_testis` gives the gene expression level in testis
#'   (using GTEx database).
#'
#' - Column `lowly_expressed_in_GTEX` indicates if the gene is lowly
#'   expressed in GTEX database and thus needed to be analysed with
#'   multimapping allowed.
#'
#' - Column `multimapping_analysis` informs if the gene (flagged as
#'   "lowly_expressed" in GTEX_data) was found to be testis-specific
#'   when multi-mapped reads were counted for gene expression in
#'   normal tissues ("not_analysed" or "testis_specific") (see
#'   `?normal_tissues_multimapping_data` for details).
#'
#' - Column `HPA_RNA_single_cell_type_specific_nTPM` specifies the cell types in
#'   which genes were detected in the HPA single cell data (see
#'   `?HPA_cell_type_specificity` for details).
#'
#' - Column `max_HPA_germcell` specifies if the maximum expression value in a
#' germ cell type. (see  `?HPA_cell_type_specificity` for details).
#'
#' - Column `max_HPA_somatic` specifies if the maximum expression value in a
#' somatic cell type. (see  `?HPA_cell_type_specificity` for details).
#'
##' - Column `not_detected_in_somatic_HPA` specifies if the gene is detected or
#' not in a somatic cell type. (see  `?HPA_cell_type_specificity` for details).
#'
#' - Column `HPA_ratio_germ_som` gives the ratio between `max_HPA_germcell` and
#' `max_HPA_somatic` columns.

#'
#' - Column `percent_of_positive_CCLE_cell_lines` gives the percentage
#'   of CCLE cancer cell lines in which genes are expressed (genes
#'   were considered as expressed if TPM >= 1).
#'
#' - Column `percent_of_negative_CCLE_cell_lines` gives the percentage
#'   of CCLE cancer cell lines in which genes are repressed (TPM <=
#'   0.5).
#'
#' - Column `max_TPM_in_CCLE` gives the highest expression level of
#'   genes in CCLE cell lines.
#'
#' - Column `CCLE_category` gives the category assigned to each gene
#'   using CCLE data. "Activated" category corresponds to genes
#'   expressed in at least 1% of cell lines (TPM >= 1) and repressed in
#'   at least 20% of cell lines.
#'
#' - Column `percent_pos_tum` gives the percentage of TCGA cancer
#'   samples in which genes are expressed (genes were considered as
#'   expressed if TPM >= 1).
#'
#' - Column `percent_neg_tum` gives the percentage of TCGA cancer samples in
#'   which genes are repressed (TPM <= 0.5).
#'
#' - Column `max_TPM_in_TCGA` gives the highest expression level of
#'   genes in TCGA cancer sample.
#'
#' - Column `max_q75_in_NT` gives the maximum q75 expression in normal
#'   peritumoral tissues from TCGA.
#'
#' - Column `TCGA_category` gives the category assigned to each gene
#'   using TCGA data. "activated" category corresponds to genes
#'   expressed in at least 1% of tumors (TPM >= 1) and repressed in at
#'   least 20% of samples.  "multimapping_issue" corresponds to genes
#'   that need multi-mapping to be allowed in order to be analysed
#'   properly.
#'
#' - Columns `external_transcript_name`, `ensembl_transcript_id`, and
#'   `transcript_biotype` give the references and informations about
#'   the most biologically relevant transcript associated to each
#'   gene.
#'
#' - Column `IGV_backbone` indicates if a gene has been removed from CT genes
#'   as RNA-Seq reads were not properly aligned on exons, but were instead
#'   spread across a wide genomic region spanning the genes.
#'
#' - Column `family` gives the gene family name.
#'
#' - Column `DAC_induced` summarises the results (`TRUE` or `FALSE`)
#'   of a differential expression evaluating gene induction upon DAC
#'   treatment in a series of cell lines.
#'
#' - Column named `CpG_density`, gives the density of CpG within each
#'   promoter (number of CpG / promoter length * 100).
#'
#' - Column `CpG_promoter` classifies the promoters according to their
#'   CpG densities: "low" (CpG_density < 2), "intermediate"
#'   (CpG_density >= 2 & CpG_density < 4), and "high" (CpG_density >=
#'   4).
#'
#' - Column `somatic_met_level` that gives the mean methylation level
#'   of each promoter in somatic tissues.
#'
#' - Column `sperm_met_level` that gives the methylation level of each
#'   promoter in sperm.
#'
#' - Column `somatic_methylation` indicates if the promoter's mean
#'   methylation level in somatic tissues is higher than 50%.
#'
#' - Column `germline_methylation`indicates if the promoter is
#'   methylated in germline, based on the ratio with somatic tissues
#'   (`FALSE` if somatic_met_level is at least twice higher than
#'   `germline_met_level`).
#'
#' - Columns `oncogene` and `tumor_suppressor` informs if oncogenic
#'   and tumor-suppressor functions have been associated to genes
#'   (source: [Cancermine](http://bionlp.bcgsc.ca/cancermine/)).
#' @source
#'
#' See `scripts/make_all_genes_prelim.R` and
#' `scripts/make_all_genes_and_CT_genes.R` for details on how this list of genes
#'  was created.
#'
#' @name all_genes
#'
#' @docType data
NULL

#' Gene correlations in CCLE cancer cell lines
#'
#' @description
#'
#' Correlation coefficients between Cancer-Testis genes and all genes
#' found on the CCLE database.
#'
#' @format
#'
#' A `matrix` object with 238 rows and 24483 columns.
#'
#' - Rows correspond to CT genes
#'
#' - Columns correspond to all genes from CCLE database
#'
#' @details
#'
#' Correlation coefficients (Pearson) between CT genes and all other
#' genes are given in the matrix. These correlation coefficients were
#' calculated using log transformed expression values from `CCLE_data`
#' (all cell lines).
#'
#' @source
#'
#' See `scripts/make_CCLE_correlation_matrix.R` for details.
#'
#' @name CCLE_correlation_matrix
#'
#' @docType data
NULL

#' Fetal gonad scRNAseq
#'
#' @description
#'
#' Human fetal gonad single cell RNAseq data from
#' `Single-cell roadmap of human gonadal development`
#' (Garcia-Alonso, Nature 2022)
#'
#' @format
#'
#' A `SingleCellExperiment` object with 22489 rows and 10850 columns
#'
#' - Rows correspond to genes (gene names as rownames)
#' - Columns correspond to cells
#'
#' @details
#'
#' Description of the colData:
#'
#' - Column `type` gives the gender and the cell type.
#'
#' - Column `stage` specifies if the cell type is "pre-meiotic" or "meiotic".
#'
#' - Column `germcell` is set to TRUE when the cell type is a germ cell.
#'
#' @source
#'
#' ee58527e-e1e4-465d-8dc8-800ee40f14f2.rds file dowloaded from
#' https://cellxgene.cziscience.com/collections/661a402a-2a5a-4c71-9b05-b346c57bc451Data.
#' The data were converted in a SingleCellExperiment
#' (see `scripts/make_FGC_sce.R` for details).
#'
#' @name FGC_sce
#'
#' @docType data
NULL

#' Oocytes scRNAseq
#'
#' @description
#'
#' Human oocytes single cell RNAseq data from
#' `Decoding dynamic epigenetic landscapes in human oocytes using single-cell
#' multi-omics sequencing` (Yan et al. Cell Stem Cell 2021)
#'
#' @format
#'
#' A `SingleCellExperiment` object with 26500 rows and 899 columns
#'
#' - Rows correspond to genes(gene names as rownames)
#'
#' - Columns correspond to cells
#'
#' @details
#'
#' Description of the colData:
#'
#' - Column `type` gives the cell type.
#'
#' - Column `stage` specifies if the cell type is "pre-meiotic" or "meiotic".
#'
#' - Column `germcell` is set to TRUE when the cell type is a germ cell.
#'
#' @source
#'
#' GSE154762_hO_scChaRM_count_matix.txt.gz` was downloaded from GEO
#' (accession: GSE154762). The data were converted in a SingleCellExperiment
#' (see `scripts/make_oocytes_sce.R` for details).
#'
#' @name oocytes_sce
#'
#' @docType data
NULL

#' Cell type specificities (from HPA)
#'
#' @description
#'
#' Cell type specificities based on scRNAseq data from the
#' Human Protein Atlas (https://www.proteinatlas.org)
#'
#' @format
#'
#' A `tibble` object with 24504 rows and 7 columns.
#'
#' - Rows correspond to genes (ensembl_gene_id)
#'
#' - Columns give genes cell type specificities
#'
#' @details
#'
#'
#' - Column `HPA_scRNAseq_celltype_specific_nTPM` gives the cell types in
#' which genes were detected (corresponds to column
#' `RNA single cell type specific nTPM` of proteinatlas.tsv file.
#'
#' - Column `max_HPA_germcell` specifies if the maximum expression value in a
#' germ cell type
#'
#' - Column `max_HPA_somatic` specifies if the maximum expression value in a
#' somatic cell type
#'
##' - Column `not_detected_in_somatic_HPA` specifies if the gene is detected or
#' not in a somatic cell type. Genes are flagged as
#' `TRUE` if the `max_HPA_somatic` value is equal to 0, and `FALSE` if
#' `max_HPA_somatic` value is > 0. `NA` is set when the original table from HPA
#' had no values for that gene.
#'
#' - Column `HPA_ratio_germ_som` gives the ratio between `max_HPA_germcell` and
#' `max_HPA_somatic` columns.

#' @source
#' `proteinatlas.tsv` was downloaded from the Human Protein Atlas
#' (https://www.proteinatlas.org)
#' See `scripts/make_HPA_cell_type_specificities.R` for details.
#'
#' @name HPA_cell_type_specificities
#'
#' @docType data
NULL

#' Embryo scRNAseq
#'
#' @description
#'
#' Human embryo single cell RNAseq data in FPKM from
#' `Single Cell DNA Methylome Sequencing of Human Preimplantation
## Embryos` (Zhu et al. 2018)
#'
#' @format
#'
#' A `SingleCellExperiment` object with 26255 rows and 50 columns
#'
#' - Rows correspond to genes (gene names as rownames)
#'
#' - Columns correspond to cells
#'
#' @details
#'
#' Description of the colData:
#'
#' - Column `embryo` gives the embryo the cell is coming from.
#'
#' - Column `stage` specifies the stage of the early embryo.
#'
#' - Column `sex` is the sex inference made using the expression of RPS4Y1. If
#' mean expression of RPS4Y1 is higher than 50 FPKM, the sample is male.
#'
#' @source
#'
#' 50 FPKM files were downloaded from GEO
#' (accession: GSE81233). The data were converted in a SingleCellExperiment
#' (see `scripts/make_embryo_sce_Zhu.R` for details).
#'
#' @name embryo_sce_Zhu
#'
#' @docType data
NULL

#' Embryo scRNAseq
#'
#' @description
#'
#' Human embryo single cell RNAseq data in RPKM from
#' `Single-Cell RNA-Seq Reveals Lineage and X Chromosome Dynamics in
# Human Preimplantation Embryos` (Petropulous et al, 2014)
#'
#' @format
#'
#' A `SingleCellExperiment` object with 26178 rows and 1481 columns
#'
#' - Rows correspond to genes (gene names as rownames)
#'
#' - Columns correspond to cells
#'
#' @details
#'
#' Description of the colData:
#'
#' - Column `individual` gives the sample the cell is coming from.
#'
#' - Column `stage` specifies the stage of the early embryo.
#'
#' - Column `sex` is the sex inference made using the expression of 11
#' Y-linked genes, made for each day individually.
#'
#' - Column `ambigous` indicates if the inference of the embryo's sex was
#' ambigous due to some cells expression of the Y-linked genes.
#'
#' @source
#'
#' RPKM and metadata files were downloaded from
#' https://www.ebi.ac.uk/biostudies/files/E-MTAB-3929/
#' The data were converted in a SingleCellExperiment
#' (see `scripts/make_embryo_sce_Petropoulos.R` for details).
#'
#' @name embryo_sce_Petropoulos
#'
#' @docType data
NULL

