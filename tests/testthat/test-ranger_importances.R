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


test_that("ranger importances from seurat objects work", {
    expect_s3_class(FindAllMarkers_ranger.seurat(Seurat::pbmc_small),
                    class = "data.frame")
})

