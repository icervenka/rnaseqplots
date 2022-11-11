library(magrittr, include.only = "%>%")

validate_config <- function(item, config_dict) {
  # empty strings are sometimes needed for patterns for file listing
  if (item %in% names(config_dict)) { # && nchar(config[[item]]) > 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

append_dir_slash <- function(dir_string) {
  if (!base::endsWith(dir_string, "/")) {
    dir_string <- paste0(dir_string, "/")
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

modify_gene_symbols = function(gene_list, 
                               text_replacements, 
                               text_modifications) {
  gene_list = gene_list %>%
    stringr::str_replace_all(regex(text_replacements, ignore_case = TRUE)) %>%
    purrr::map_chr(purrr::compose(!!!text_modifications))
  return(gene_list)
}

mouseify_gene_symbols = function(gene_list) {
  modify_gene_symbols(gene_list,
                      text_replacements = c(
                        mito_atp_replacer,
                        c("MT-CO" = "MT-COX", "MT-" = "")
                      ),
                      text_modifications = list(tools::toTitleCase, tolower))
}
