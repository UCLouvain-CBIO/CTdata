test_that("makeTags() works", {
    res0 <- c("ExperimentHub", "ExperimentData",
              "ReproducibleResearch", "RepositoryData",
              "Homo_sapiens_Data")
    res0 <- paste0(res0, collapse = ":")
    expect_identical(CTdata:::makeTags(), res0)
    res1 <- paste0(c(res0, "Tag1"), collapse = ":")
    expect_identical(CTdata:::makeTags("Tag1"), res1)
    res2 <- paste0(c(res1, "Tag2"), collapse = ":")
    expect_identical(CTdata:::makeTags(c("Tag1", "Tag2")), res2)
})
