config <- rjson::fromJSON(file = "config.json") %>%
  purrr::map_dfr(data.frame)
