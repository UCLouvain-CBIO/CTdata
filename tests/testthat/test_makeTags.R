test_that("makeTags() works", {
    makeTags <- CTdata:::makeTags
    res0 <- c("ExperimentHub", "ExperimentData",
              "ReproducibleResearch", "RepositoryData",
              "Homo_sapiens_Data")
    res0 <- paste0(res0, collapse = ":")
    expect_identical(makeTags(), res0)
    res1 <- paste0(c(res0, "Tag1"), collapse = ":")
    expect_identical(makeTags("Tag1"), res1)
    res2 <- c(res1,
              paste0(c(res0, "Tag2"), collapse = ":"))
    expect_identical(makeTags(c("Tag1", "Tag2")), res2)
})
