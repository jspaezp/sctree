context("test-tree_rules")

test_that("Tree rules are handled correctly on generic trees", {
    iris_tree <- partykit::ctree(Species ~ ., data = iris)
    my_rules <- get_concensus_rules(iris_tree)
    expect_output({
        print.concensus.rules(my_rules)
    })
})
