gg_color_hue <- function(i) {
  hues <- seq(15, 375, length = i + 1)
  hcl(h = hues, l = 65, c = 100)[1:i]
}

rank_how = function(descending) {
  if(descending == T) {
    return(dplyr::desc)
  } else  {
    return(identity)
  }
}

get_biomart_gene_mapping = function(mart, id_attribute, mapped_attribute, ...) {
  if(!all(c(id_attribute %in% list(...)[["attributes"]],
            mapped_attribute %in% list(...)[["attributes"]]))) {
    stop("Id attribute or mapped attribute not present in attribute list.")
  }

  mapping = getBM(mart = mart, ...)
  all_values = data.frame(id = list(...)[["values"]], stringsAsFactors = F) %>%
    left_join(mapping, by = c("id" = id_attribute)) %>%
    mutate(!!as.symbol(mapped_attribute) := case_when(is.na(!!as.symbol(mapped_attribute)) ~ id,
                                                      !!as.symbol(mapped_attribute) == "" ~ id,
                                                      TRUE ~ !!as.symbol(mapped_attribute))) %>%
    group_by(id) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    mutate(!!as.symbol(id_attribute) := id, .before = id) %>%
    dplyr::select(-id)

  return(all_values)
}
