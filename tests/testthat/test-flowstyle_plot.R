context("test-flowstyle_plot")

test_that("works", {
    expect_s3_class({
        g <- plot_flowstyle(Seurat::pbmc_small, c("ACRBP", "TSC22D1", "VDAC3"))
        g
    }, "gg")
})
