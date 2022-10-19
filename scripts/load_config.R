config_path = "config.json"

config <- rjson::fromJSON(file = config_path)[["plot_export_params"]] %>%
  purrr::map_dfr(data.frame)

group_levels <- rjson::fromJSON(file = config_path)[["group_levels"]]