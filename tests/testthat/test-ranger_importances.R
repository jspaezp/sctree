context("test-ranger_importances")

test_that("ranger importances from Seurat objects work", {

    expect_warning({
        expect_s3_class(
            ranger_importances.Seurat(
                Seurat::pbmc_small,
                cluster = "ALL",
                warn.imp.method = FALSE
            )[[3]],
            "data.frame"
        )

        expect_s3_class(
            ranger_importances.Seurat(
                Seurat::pbmc_small,
                genes_use = rownames(Seurat::pbmc_small@assays[[Seurat::pbmc_small@active.assay]]@data),
                cluster = "ALL")[[1]],
            "ranger"
        )
    })

    expect_warning(
        ranger_importances.Seurat(Seurat::pbmc_small, cluster = "ALL")
    )

    expect_silent(
        ranger_importances.Seurat(
            Seurat::pbmc_small,
            cluster = "ALL",
            warn.imp.method = FALSE
        )
    )

    test_that(
        paste0("ranger importances match the names ",
               "in the original seurat object"), {
        expect_warning({
            expect_equal(
                rownames(ranger_importances.Seurat(
                    Seurat::pbmc_small,
                    genes_use = rownames(Seurat::pbmc_small@assays[[Seurat::pbmc_small@active.assay]]@data),
                    cluster = "ALL")[[2]]),
                rownames(Seurat::pbmc_small@assays[[Seurat::pbmc_small@active.assay]]@data))
        })
    })


})


test_that(
    paste0("ranger importances does not throw warning ",
           "when called with the altman method"), {
    data_names <- rownames(
        Seurat::pbmc_small@assays[[Seurat::pbmc_small@active.assay]]@data)

    expect_silent({
        altmann_output <- ranger_importances.Seurat(
            Seurat::pbmc_small,
            genes_use = data_names,
            cluster = "ALL", imp_method = "altmann",
            num.trees = 10)
    })

    expect_s3_class(altmann_output[[3]],"data.frame")

    expect_s3_class(altmann_output[[1]],"ranger")


    test_that(
        paste0("ranger importances match the names in the",
               " original Seurat object with the altmann method"),{
        expect_equal(
            rownames(altmann_output[[2]]),
            data_names)

    })


    expect_silent({
        altmann_output <- ranger_importances.Seurat(
            Seurat::pbmc_small,
            genes_use = data_names,
            cluster = "0", imp_method = "altmann",
            num.trees = 10)
    })

    expect_s3_class(altmann_output[[3]],"data.frame")

    expect_s3_class(altmann_output[[1]],"ranger")

    test_that(
        paste0("ranger importances match the names in the ",
               "original Seurat object with the altmann method"),{
        expect_equal(
            rownames(altmann_output[[2]]),
            data_names)

    })
})


test_that("ranger importances from Seurat objects work", {
    expect_warning({
        expect_s3_class(FindAllMarkers_ranger.Seurat(Seurat::pbmc_small),
                        class = "data.frame")
    })

    expect_s3_class(
        FindAllMarkers_ranger.Seurat(
            Seurat::pbmc_small,
            imp_method = "altmann",
            num.trees = 10),
        class = "data.frame")
})

