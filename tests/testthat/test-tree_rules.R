context("test-tree_rules")

test_that("Tree rules are handled correctly on generic trees", {
    iris_tree <- partykit::ctree(Species ~ ., data = iris)
    my_rules <- get_concensus_rules(iris_tree)
    expect_output({
        print.concensus.rules(my_rules)
    })
})


expected_garnett_text <- c(
    "> 0_node_3 \t(n = 40)",
    "expressed below: ADA 2.84, ASNS 1.43",
    "",
    "> 0_node_6 \t(n = 8)",
    "expressed above: ASNS 1.43, ITGA4 0",
    "expressed below: ADA 2.84",
    "",
    "> 0_node_8 \t(n = 135)",
    "expressed above: ADA 2.84",
    "expressed below: TSC22D3 2.12",
    "",
    "> 0_node_9 \t(n = 18)",
    "expressed above: ADA 2.84, TSC22D3 2.12",
    "",
    "> 1 \t(n = 54)",
    "expressed above: ASNS 1.43",
    "expressed below: ADA 2.84, ITGA4 0"
)




test_that("Tree rules are handled correctly on generic trees", {

    actual_out <- capture.output(
        as.garnett(fit_ctree(sctree::small_5050_mix), digits = 2))
    expect_equal(actual_out, expected_garnett_text)
})

test_that("Tree rules are handled correctly on generic trees", {
    actual_out <- capture.output(
        as.garnett(fit_ctree(sctree::small_5050_mix),
                   digits = 2))
    expect_equal(actual_out, expected_garnett_text)
})

test_that("Tree rules filters correctly based on regexes", {
    expected_out <- c(
        "> 1 \t(n = 54)",
        "expressed above: ASNS 1.43",
        "expressed below: ADA 2.84, ITGA4 0")
    actual_out <- capture.output(
        as.garnett(fit_ctree(sctree::small_5050_mix),
                   digits = 2,
                   rules_keep = "^1"))
    expect_equal(actual_out, expected_out)
})

