context("test-ranger_importances")

test_that("ranger importances from seurat objects work", {

    expect_warning({
        expect_s3_class(
            ranger_importances.seurat(
                Seurat::pbmc_small,
                cluster = "ALL",
                warn.imp.method = FALSE
            )[[3]],
            "data.frame"
        )

        expect_s3_class(
            ranger_importances.seurat(
                Seurat::pbmc_small,
                genes_use = rownames(Seurat::pbmc_small@data),
                cluster = "ALL")[[1]],
            "ranger"
        )
    })

    expect_warning(
        ranger_importances.seurat(Seurat::pbmc_small, cluster = "ALL")
    )

    expect_silent(
        ranger_importances.seurat(
            Seurat::pbmc_small,
            cluster = "ALL",
            warn.imp.method = FALSE
        )
    )

    test_that("ranger importances match the names in the original seurat object", {
        expect_warning({
            expect_equal(
                rownames(ranger_importances.seurat(
                    Seurat::pbmc_small,
                    genes_use = rownames(Seurat::pbmc_small@data),
                    cluster = "ALL")[[2]]),
                rownames(Seurat::pbmc_small@data))
        })
    })


})


test_that("ranger importances does not throw warning when called with the altman method", {
    expect_silent({
        altmann_output <- ranger_importances.seurat(
            Seurat::pbmc_small,
            genes_use = rownames(Seurat::pbmc_small@data),
            cluster = "ALL", imp_method = "altmann")
    })

    expect_s3_class(altmann_output[[3]],"data.frame")

    expect_s3_class(altmann_output[[1]],"ranger")

    test_that("ranger importances match the names in the original seurat object with the altmann method",{
        expect_equal(
            rownames(altmann_output[[2]]),
            rownames(Seurat::pbmc_small@data))

    })


    expect_silent({
        altmann_output <- ranger_importances.seurat(
            Seurat::pbmc_small,
            genes_use = rownames(Seurat::pbmc_small@data),
            cluster = "0", imp_method = "altmann")
    })

    expect_s3_class(altmann_output[[3]],"data.frame")

    expect_s3_class(altmann_output[[1]],"ranger")

    test_that("ranger importances match the names in the original seurat object with the altmann method",{
        expect_equal(
            rownames(altmann_output[[2]]),
            rownames(Seurat::pbmc_small@data))

    })
})




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
    names_5050 <- rownames(small_5050_mix@data)
    standard_names_version <- make.names(names_5050)
    non_standard_names <- names_5050[names_5050 != standard_names_version]

    mix_standard_names <- c(small_5050_mix@var.genes, non_standard_names)

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
