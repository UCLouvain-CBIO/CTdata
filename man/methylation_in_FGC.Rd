% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{methylation_in_FGC}
\alias{methylation_in_FGC}
\title{Methylation of CpGs within all genes promoters in FGC}
\format{
A \code{RangedSummarizedExperiment} object with 1915545 rows and 337 columns
\itemize{
\item Rows correspond to CpGs (located within all genes promoters (TSS +- 1000
nt))
\item Columns correspond to cells
\item Methylation values from scWGBS data
\item rowRanges correspond to CpG positions
}
}
\source{
scWGBS methylation data was downloaded from GEO database
(see \code{scripts/make_methylation_in_FGC.R} for details).
}
\description{
Methylation values of CpGs located within all genes
promoters in fetal germ cells. Data is based on hg19 reference genome ! From
\verb{Dissecting the epigenomic dynamics of human fetal germ cell development at single-cell resolution (Li et al. 2021)}
}
\details{
Description of the colData:
\itemize{
\item Column \code{type} indicates if the cell type is somatic or FGC
\item Column \code{time_week} specifies the time of the embryo when cells were removed.
\item Column \code{sex} indicates the sex of the cells.
\item Other information about the sequencing of each sample are clearly labelled
}
}
