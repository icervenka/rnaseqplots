path_list <- rjson::fromJSON(file = "../paths.json")

input_directory <- path_list$input_directory %>%
    dir_exists() %>%
    append_dir_slash()

output_directory <- path_list$output_directory %>%
    dir_exists() %>%
    append_dir_slash()

purrr::walk2(
    c(
        names(path_list$input_data_paths),
        names(path_list$param_paths)
    ),
    c(
        path_list$input_data_paths,
        path_list$param_paths
    ), function(x, y) {
        assign(x, y, envir = .GlobalEnv)
    }
)
