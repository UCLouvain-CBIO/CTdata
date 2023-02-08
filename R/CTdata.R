##' @title All CTdata datasets
##'
##' @description
##'
##' This is the companion Package for `CTexploreR` containing omics
##' data to select and characterise CT genes.
##'
##' Data come from public databases and include expression and
##' methylation values of genes in normal and tumor samples as well as
##' in tumor cell lines, and expression in cells treated with a
##' demethylating agent is also available.
##'
##' The [CTdata()] function returns a `data.frame` with all the
##' annotated datasets provided in the package. For details on these
##' individual datasets, refer to their respective manual pages.
##'
##' See the vignette and the respective manuals pages for more details
##' about the package and the data themselves.
##'
##' @return A `data.frame` describing the data available in
##'     `CTdata`.
##'
##' @author Laurent Gatto
##'
##' @export
##'
##' @importFrom utils read.csv
##'
##' @examples
##'
##' MsDataHub()
CTdata <- function() {
   fl <- system.file("extdata", "metadata.csv", package = "CTdata")
   read.csv(fl, stringsAsFactors = FALSE)
}
