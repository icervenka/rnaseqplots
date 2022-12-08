expression_data <- load_expression_data(
    "../inst/extdata/rnaseq/expression_data.txt",
    list()
)

usethis::use_data(expression_data, overwrite = TRUE)