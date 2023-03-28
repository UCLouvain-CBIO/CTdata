##' A short function that returns the default CTdata tags and, if
##' provided, additional data-specific tags.
##'
##' @param x An optional `character()` containing specific tags.
##'
##' @return
##'
##' A `character` containing the default tags and optional
##' data-specific tags. If `x` is missing or is of length 0, the
##' default tags are returned. Otherwise, a vector of length equal to
##' `length(x)` is returned.
##'
##' @examples
##'
##' CTdata:::makeTags() ## only default tags
##'
##' CTdata:::makeTags(character()) ## only default tags
##'
##' CTdata:::makeTags("myTag") ## one additional tag
##'
##' CTdata:::makeTags(c("myTag", "myOtherTag")) ## two additional tag
makeTags <- function(x) {
    defaultTags <- c("ExperimentHub", "ExperimentData",
                      "ReproducibleResearch", "RepositoryData",
                     "Homo_sapiens_Data")
    defaultTags <- paste(defaultTags, collapse = ":")
    if (missing(x) || !length(x))
        return(defaultTags)
    vapply(x, function(.x) paste(c(defaultTags, .x), collapse = ":"),
           USE.NAMES = FALSE, "")
}
