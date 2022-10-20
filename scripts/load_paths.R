path_list <- rjson::fromJSON(file = "paths.json")

output_directory <- path_list$output_directory %>%
  dir_exists() %>%
  append_dir_slash()

# load paths for input files
purrr::walk2(
  names(path_list$input_data_paths), path_list$input_data_paths,
  function(x, y) {
    assign(x, y, envir = .GlobalEnv)
  }
)

# load paths for params
purrr::walk2(
  names(path_list$params), path_list$params,
  function(x, y) {
    assign(x, y, envir = .GlobalEnv)
  }
)
