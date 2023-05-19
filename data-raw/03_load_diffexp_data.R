diffexp_data <- load_data_diffexp(
    list(
        ex_1 = "inst/extdata/rnaseq/diffexp_data_1.txt",
        ex_2 = "inst/extdata/rnaseq/diffexp_data_2.txt"
    ),
    list()
)

usethis::use_data(diffexp_data, overwrite = TRUE)
