context("test-flowstyle_plot")

test_that("works", {
    expect_s3_class({
        g <- plot_flowstyle(sctree::small_5050_mix, c("ASNS", "ITGA4", "ADA"))
        g
    }, "gg")
})
