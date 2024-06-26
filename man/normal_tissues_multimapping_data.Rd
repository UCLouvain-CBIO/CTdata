% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{normal_tissues_multimapping_data}
\alias{normal_tissues_multimapping_data}
\title{Gene expression values in normal tissues}
\format{
A \code{SummarizedExperiment} object with 24504 rows and 18 columns
\itemize{
\item Rows correspond to genes (ensembl_gene_id)
\item Columns correspond to normal tissues.
\item First assay, \code{TPM_no_multimapping}, gives TPM expression values
obtained when discarding multimapped reads.
\item Second assay, \code{TPM_with_multimapping}, gives TPM expression
values obtained by counting multimapped reads.
}
}
\source{
RNAseq fastq files were downloaded from Encode database (see
\code{scripts/make_normal_tissues_multimapping.R} for details).
}
\description{
Gene expression values (TPM) in a set of normal tissues obtained by
counting or not multi-mapped reads. Many CT genes belong to gene
families from which members have identical or nearly identical
sequences.  Some CT can only be detected in RNAseq data in which
multimapping reads are not discarded.
}
\details{
A column named \code{multimapping_analysis} has been added to the
rowData. It summarizes the testis specificity analysis of genes
flagged as "lowly_expressed" in GTEX_data. Genes are considered
"testis_specific" when, with multimapping allowed, they are
detectable in testis (TPM >= 1), their TPM value has increased
compared to without multimapping (ratio > 5), and their TPM value
is at least 10 times higher in testis than in any other somatic
tissue (where the maximum expression always has to be below 1 TPM).
Genes are considered "testis_preferential" when, with multimapping allowed,
they are detectable in testis (TPM >= 1), their TPM value has increased
compared to without multimapping (ratio > 5), and their TPM value
is at least 10 times higher in testis than in any other somatic
tissue (where the maximum expression is above 1 TPM).
}
