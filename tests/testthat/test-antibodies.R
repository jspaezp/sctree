context("test-antibodies")

test_that("functions work on known markers", {
    expect_gt(nrow(query_cc_antibodies("CD11c")), 1)
    expect_gt(nrow(query_cc_antibodies("CD3")),1)

    expect_gt(nrow(head(query_sc_antibodies("CD11C"))), 1)
    expect_gt(nrow(head(query_sc_antibodies("CD3"))),1)

    expect_gt(length(query_biolegend_antibodies("CD11C")), 1)
    expect_gt(length(query_biolegend_antibodies("CD3")), 1)
})


test_that("functions return null or empty df on imposible markers", {
    expect_null(query_cc_antibodies("SUPERFAKEMARKER"))

    expect_null(query_sc_antibodies("SUPERFAKEMARKER"))

    expect_length(query_biolegend_antibodies("SUPERFAKEMARKER"), 0)
})
