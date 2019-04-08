context("test-utils")

test_that("top_n returns the top and bottom values when given -1 and 1 for n", {
    expect_equal(unique(top_n(iris, -1, "Petal.Width")$Petal.Width), 0.1)
    expect_equal(unique(top_n(iris, 1, "Petal.Width")$Petal.Width), 2.5)
})


test_that("tsne plot wrapper works", {
    expect_s3_class(tsne_plot(Seurat::pbmc_small), "gg")
})


test_that("converting seurat to data frame works", {
    data("small_5050_mix")
    num_genes <- length(rownames(small_5050_mix@data))
    num_vargenes <- length(small_5050_mix@var.genes)
    expect_s3_class({
        as.data.frame(small_5050_mix, genes = rownames(small_5050_mix@data))
    }, "data.frame")

    expect_equal({
        ncol(as.data.frame(
            small_5050_mix,
            genes = small_5050_mix@var.genes))
    }, num_vargenes + 1)

    expect_equal({
        colnames(as.data.frame(
            small_5050_mix,
            genes = small_5050_mix@var.genes))
    }, c(small_5050_mix@var.genes, "ident"))

    expect_equal({
        colnames(as.data.frame(
            small_5050_mix,
            genes = small_5050_mix@var.genes,
            fix_names = TRUE))
    }, make.names(c(small_5050_mix@var.genes, "ident")))
})


test_that("filtering membrane genes is successfull", {
    expect_type({
        is_gene_membrane(c("CD4", "UPP1"), evidence_codes = c("EXP", "IDA"))
    }, "logical")

    expect_equal({
        is_gene_membrane(c("CD4", "UPP1"), evidence_codes = c("EXP", "IDA"))
    }, c(TRUE, FALSE))
})


test_that("Getting aliases is successfull", {
    expect_type({
        get_aliases(c("MAPK1", "CD4", "UPP1"))
    }, "list")

    gotten_aliases <- get_aliases(c("MAPK1", "CD4"))

    expect_type({
        gotten_aliases[[1]]
    }, "character")

    expect_type({
        gotten_aliases[[2]]
    }, "character")

    expect_equal({
        names(gotten_aliases)
    }, c("MAPK1", "CD4"))
})


test_that("Frequency matrix operations are consistent", {
    expect_warning(as.frequency.matrix(table(1:3, 1:3)/3.33))
    expect_warning(as.frequency.matrix(matrix(1:10, ncol = 2)))

    expect_silent(as.frequency.matrix(table(1:3, 1:3)))

    expect_equal(as.vector(as.frequency.matrix(table(1:3, 1:3))),
                 c(100, 0, 0, 0, 100, 0, 0, 0, 100))
})


test_that("Autoplot functions fulfill expectations", {
    expect_s3_class(autoplot(table(1:3, 1:3)/3.33), "gg")
    expect_s3_class(autoplot(matrix(1:10, ncol = 2)), "gg")

    g <- autoplot.table(table(1:3, 1:3), min_color = 50)
    expect_length(g$layers, 3)
    expect_s3_class(g$layers[[1]]$geom, "GeomTile")
    expect_s3_class(g$layers[[2]]$geom, "GeomTile")

    g <- autoplot.table(table(1:3, 1:3))
    expect_length(g$layers, 1)
    expect_s3_class(g$layers[[1]]$geom, "GeomTile")

    expect_equal(as.vector(as.frequency.matrix(table(1:3, 1:3))),
                 c(100, 0, 0, 0, 100, 0, 0, 0, 100))
})


