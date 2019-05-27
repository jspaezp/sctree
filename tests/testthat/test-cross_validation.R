context("test-cross_validation")

test_that("cross validation works in the bundled datasets", {
    data(small_5050_mix)
    data(small_9901_mix)

    expect_warning({
        expect_s3_class(
            cross_validate(small_5050_mix, small_9901_mix,
                           cluster = "ALL",
                           warn.imp.method = FALSE)[[1]],
            "party")
    })

    expect_warning({
        expect_s3_class(
            cross_validate(small_5050_mix, small_9901_mix,
                           cluster = "1",
                           warn.imp.method = FALSE)[[1]],
            "party")
    })

})


test_that("cross validation handles correctly non-standard names",{
    names_5050 <- rownames(small_5050_mix@assays$RNA@data)
    standard_names_version <- make.names(names_5050)
    non_standard_names <- names_5050[names_5050 != standard_names_version]

    mix_standard_names <- c(Seurat::VariableFeatures(small_5050_mix),
                            non_standard_names)

    expect_warning({
        expect_s3_class(
            cross_validate(small_5050_mix, small_9901_mix,
                           cluster = "ALL",
                           genes_use = mix_standard_names,
                           warn.imp.method = FALSE)[[1]], "party")
    })

    expect_silent({
        expect_s3_class(
            cross_validate(small_5050_mix, small_9901_mix,
                           cluster = "ALL",
                           genes_use = mix_standard_names,
                           warn.gene.removal = FALSE,
                           warn.imp.method = FALSE)[[1]], "party")
    })



})
