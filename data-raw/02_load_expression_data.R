expression_data <- load_data_expression(
    "inst/extdata/rnaseq/expression_data.txt",
    list()
)

usethis::use_data(expression_data, overwrite = TRUE)
