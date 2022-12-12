cp_pathways <- load_pathways_cp(
    "inst/extdata/cp_pathways",
    list(cp_pathways_pattern = "_contrast")
)

usethis::use_data(cp_pathways, overwrite = TRUE)
