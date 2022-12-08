diffexp_data <- load_diffexp_data(
    "../inst/extdata/rnaseq/diffexp_data.txt",
    list()
)

usethis::use_data(diffexp_data, overwrite = TRUE)