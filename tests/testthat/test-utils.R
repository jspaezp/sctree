context("test-utils")

test_that("top_n returns the top and bottom values when given -1 and 1 for n", {
    expect_equal(unique(top_n(iris, -1, "Petal.Width")$Petal.Width), 0.1)
    expect_equal(unique(top_n(iris, 1, "Petal.Width")$Petal.Width), 2.5)
})
