config_path = "config.json"
config <- rjson::fromJSON(file = config_path) %>%
  purrr::map_dfr(data.frame)
