context("test-ranger_importances")


test_that("ranger importances from Seurat objects work", {

    expect_warning({
        expect_s3_class(
            ranger_importances.Seurat(
                sctree::small_5050_mix,
                cluster = "ALL",
                warn.imp.method = FALSE
            ),
            "data.frame"
        )

        expect_s3_class(
            ranger_importances.Seurat(
                sctree::small_5050_mix,
                genes_use = rownames(sctree::small_5050_mix@assays[[sctree::small_5050_mix@active.assay]]@data),
                cluster = "ALL", return_what = "ranger_fit"),
            "ranger"
        )
    })

    expect_warning(
        ranger_importances.Seurat(sctree::small_5050_mix, cluster = "ALL")
    )

    expect_silent(
        ranger_importances.Seurat(
            sctree::small_5050_mix,
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
                    sctree::small_5050_mix,
                    genes_use = rownames(sctree::small_5050_mix@assays[[sctree::small_5050_mix@active.assay]]@data),
                    cluster = "ALL", return_what = "importances_ranger")),
                rownames(sctree::small_5050_mix@assays[[sctree::small_5050_mix@active.assay]]@data))
        })
    })


})


test_that(
    paste0("ranger importances does not throw warning ",
           "when called with the altman method"), {
    data_names <- rownames(
        sctree::small_5050_mix@assays[[sctree::small_5050_mix@active.assay]]@data)

    expect_silent({
        altmann_output <- ranger_importances.Seurat(
            sctree::small_5050_mix,
            genes_use = data_names,
            cluster = "ALL", imp_method = "altmann",
            num.trees = 10, return_what = c("ranger_fit", "importances_ranger"))
    })

    expect_s3_class(altmann_output[[2]], "data.frame")

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
            sctree::small_5050_mix,
            genes_use = data_names,
            cluster = "0", imp_method = "altmann",
            num.trees = 10, return_what = c("ranger_fit", "importances_ranger"))
    })

    expect_s3_class(altmann_output[[2]],"data.frame")

    expect_s3_class(altmann_output[[1]],"ranger")

    test_that(
        paste0("ranger importances match the names in the ",
               "original Seurat object with the altmann method"),{
        expect_equal(
            rownames(altmann_output[[2]]),
            data_names)

    })
})


test_that("FindAllMarkers_ranger from Seurat objects work", {
    expect_warning({
        expect_s3_class(FindAllMarkers_ranger.Seurat(sctree::small_5050_mix),
                        class = "data.frame")
    })

    expect_s3_class(
        FindAllMarkers_ranger.Seurat(
            sctree::small_5050_mix,
            imp_method = "altmann",
            num.trees = 10),
        class = "data.frame")
})


test_that("FindMarkers interface with RangerDE works", {
    expect_warning({
        sctree::FindMarkers(
            object = sctree::small_5050_mix, ident.1 = 0, ident.2 = 1,
            test.use = "RangerDE", verbose = FALSE)
    })

    expect_s3_class(
        sctree::FindMarkers(
            object = sctree::small_5050_mix, ident.1 = 0, ident.2 = 1,
            test.use = "RangerDE", verbose = FALSE,
            imp_method = "altmann",
            num.trees = 10),
        class = "data.frame")
})
