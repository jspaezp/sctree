context("test-examples")

documentation_files <- dir(
    "../../man",
    pattern = ".*.Rd",
    full.names = TRUE)

if (Sys.getenv("TRAVIS") == "true") {
    skip_tests <- "antibodies.Rd"
} else {
    skip_tests <- "NULL"
}

documentation_files <- documentation_files[
    !grepl(skip_tests, documentation_files)]

sink(tempfile(), type = "output")

for (doc in documentation_files) {
    print(doc)
    suppressWarnings(test_example(path = doc, doc))
}
sink()
