##' A short function that returns the default CTdata tags and, if
##' provided, additional data-specific tags.
##'
##' @param x `character()` containing specific tags. Default is an
##'     empty vector.
##'
##' @return
##'
##' A `character` containing the default tags and optional
##' data-specific tags.
##'
##' @examples
##'
##' CTdata:::makeTags() ## only default tags
##'
##' CTdata:::makeTags("myTag") ## one additional tag
##'
##' CTdata:::makeTags(c("myTag", "myOtherTag")) ## two additional tag
makeTags <- function(x = character()) {
    defaultTags <- c("ExperimentHub", "ExperimentData",
                      "ReproducibleResearch", "RepositoryData",
                      "Homo_sapiens_Data")
    paste(c(defaultTags, x), collapse = ":")
}
