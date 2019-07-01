context("test-tree_rules")

test_that("Tree rules are handled correctly on generic trees", {
    iris_tree <- partykit::ctree(Species ~ ., data = iris)
    my_rules <- get_concensus_rules(iris_tree)
    expect_output({
        print.concensus.rules(my_rules)
    })
})

test_that("Tree rules are handled correctly on generic trees", {
    expected_out <- c("> 0_node_2 \t(n = 52)",
                      "expressed below: S100A8 2.523156531922",
                      "",
                      "> 0_node_4 \t(n = 7)",
                      "expressed above: S100A8 2.523156531922",
                      "expressed below: S100A9 0",
                      "",
                      "> 1 \t(n = 21)",
                      "expressed above: S100A8 2.523156531922, S100A9 0")
    actual_out <- capture.output(
        as.garnett(fit_ctree(Seurat::pbmc_small), digits = Inf))
    expect_equal(actual_out, expected_out)
})

test_that("Tree rules are handled correctly on generic trees", {
    expected_out <- c("> 0_node_2 \t(n = 52)",
                      "expressed below: S100A8 2.52",
                      "",
                      "> 0_node_4 \t(n = 7)",
                      "expressed above: S100A8 2.52",
                      "expressed below: S100A9 0",
                      "",
                      "> 1 \t(n = 21)",
                      "expressed above: S100A8 2.52, S100A9 0")
    actual_out <- capture.output(
        as.garnett(fit_ctree(Seurat::pbmc_small),
                   digits = 2))
    expect_equal(actual_out, expected_out)
})

test_that("Tree rules filters correctly based on regexes", {
    expected_out <- c("> 0_node_2 \t(n = 52)",
                      "expressed below: S100A8 2.52",
                      "",
                      "> 0_node_4 \t(n = 7)",
                      "expressed above: S100A8 2.52",
                      "expressed below: S100A9 0")
    actual_out <- capture.output(
        as.garnett(fit_ctree(Seurat::pbmc_small),
                   digits = 2,
                   rules_keep = "^0"))
    expect_equal(actual_out, expected_out)
})
