metadata <- load_metadata(
    "inst/extdata/rnaseq/metadata.tsv",
    list(
        "sample_colname" = "sample",
        "group_colname" = "group",
        "group_levels" = c("ctrl", "ne",  "cgrp",  "npy", "sp")
    )
)

usethis::use_data(metadata, overwrite = TRUE)
