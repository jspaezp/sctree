context("test-examples")


test_that("Examples run without errors", {
    expect_success({
        {
            sink(tempfile(), type = "output")
            suppressWarnings(test_examples(path = "../.."))
            sink()
        }
    })
})
