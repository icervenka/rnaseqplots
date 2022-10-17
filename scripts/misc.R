library(magrittr, include.only = "%>%")

append_dir_slash <- function(dir_string) {
    if (!base::endsWith(dir_string, "/")) {
        dir_string = paste0(dir_string, "/")
    }
    return(dir_string)
}

dir_exists <- function(dir_string) {
    if (!file.exists(dir_string)) {
        stop(paste0("Specified directory does not exist: ", dir_string))
    } else {
        return(dir_string)
    }
}

read_data <- function(filename) {
  ext <- tools::file_ext(filename)
  if (ext %in% c("tsv", "txt")) {
    separator = "\t"
  } else if(ext %in% c("csv")) {
    separator = ","
  }
  data <- readr::read_delim(filename,
  delim = separator,
  col_names = TRUE,
  escape_double = FALSE,
  trim_ws = TRUE)
  return(data)
}

gg_color_hue <- function(i) {
  hues <- seq(15, 375, length = i + 1)
  hcl(h = hues, l = 65, c = 100)[1:i]
}

rank_how <- function(descending) {
  if (descending == TRUE) {
    return(dplyr::desc)
  } else {
    return(identity)
  }
}

get_biomart_gene_mapping <- function(mart,
                                     id_attribute,
                                     mapped_attribute,
                                     ...) {
  if (!all(c(
    id_attribute %in% list(...)[["attributes"]],
    mapped_attribute %in% list(...)[["attributes"]]
  ))) {
    stop("Id attribute or mapped attribute not present in attribute list.")
  }

  mapping <- biomart::getBM(mart = mart, ...)
  all_values <- data.frame(
    id = list(...)[["values"]],
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(mapping, by = c("id" = id_attribute)) %>%
    dplyr::mutate(!!as.symbol(mapped_attribute) := dplyr::case_when(
      is.na(!!as.symbol(mapped_attribute)) ~ id,
      !!as.symbol(mapped_attribute) == "" ~ id,
      TRUE ~ !!as.symbol(mapped_attribute)
    )) %>%
    dplyr::group_by(id) %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(!!as.symbol(id_attribute) := id, .before = id) %>%
    dplyr::select(-id)

  return(all_values)
}
