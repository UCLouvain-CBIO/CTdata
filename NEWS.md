# CTdata 1.5

## CTdata 1.5.1

- Changed selection criteria of testis-specific and testis-preferential genes 
in GTEX 
- Added a slightly more stringent criteria of selection of testis-specific 
genes when analysing RNAseq data of normal tissues with multumapping 
- Minor modifications in selection of genes induced by DAC
- Set TCGA_catgeory in TCGA_TPM to "leaky" when the q75 expression in normal
peritumoral TCGA samples is higher than 0.5
- Added HPA_cell_type_specificity (a table giving cell type specificity of each
genes based on Human protein Atlas scRNAseq data)
- Added FGC_sce, scRNAseq data of human fetal gonads
- Added oocytes_sce, scRNAseq data of human oocytes
- Changed preliminary CT_list into all_genes_prelim list to add the 
characterisation to all genes and not only CT in the package
- Changed criteria for testis specificity, including HPA/TCGA/CCLE category in 
the definition
- Added CT_gene_type column to differentiate CT genes and CT preferential genes
- Included all genes in methylation objects (methylation_in_tissues,
mean_methylation_in_tissues and TCGA_methylation)
- Added all_genes containing all the characterisation and analysis for all genes
- Changed the criteria for regulation by methylation, needing DAC induction and 
methylation in somatic only



## CTdata 1.5.0

- New devel


# CTdata 1.4

## CTdata 1.4.0

- New release Bioc 3.19

# CTdata 1.1

## CTdata 1.1.5

- Assay from `CT_methylation_in_tissues()` converted from a tibble to
  a matrix.
- Suggest and load SingleCellExperiment in vignette.
- Update vignette figure.

## CTdata 1.1.4

- New data: `scRNAseq_HPA` and `testis_sce`.
- Updated data: `CT_methylation_in_tissues`,
  `CT_mean_methylation_in_tissues`, `TCGA_CT_methylation`, `CT_genes`
  and `CCLE_correlation_matrix`.

## CTdata 1.1.3

- Load SummarizedExperiment in vignette.

## CTdata 1.1.2

- Suggest SummarizedExperiment.

## CTdata 1.1.1

- Only display a small subset of `CCLE_correlation_matrix()` in the
  vignette.

## CTdata 1.1.0

- New Bioconductor devel release.

# CTdata 0.99

## CTdata 0.99.0

- Package submission to Biocoductor.
