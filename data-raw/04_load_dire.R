dire <- load_pathways_dire(
    "inst/extdata/rnaseq/dire",
    list(dire_pattern = "", dire_sheet_name = "dire_all")
)

usethis::use_data(dire, overwrite = TRUE)
