#' Wrapper function for venn.diagram
#'
#' Sets some reasonable defaults for looks. If plot_params$filename is specified
#' venn diagram is saved to this location, otherwise it is just drawn on
#' screen.
#'
#' @param l list of character vectors to plot in venn diagram.
#' @param plot_params list of output parameters for plot such as size, dpi,
#' units, file format, file name for export.
#' @param category_names character vector the same length as l containing
#' category names.
#' @param palette character vector the same length as l containg colors to use
#' as fills for categories.
#' @param font_size  integer, font size for labels inside the venn diagram.
#' @param ... other parameters to VennDiagram::venn.diagram function.
#'
#' @return NULL
#' @export
#'
#' @examples
plot_venn <- function(l,
                      plot_params,
                      category_names = "",
                      palette = viridis::viridis(2),
                      font_size = 0.6,
                      ...) {
  if (nchar(plot_params$filename) == 0) grid.newpage()

  out_filename <- paste0(
    output_directory,
    plot_params$filename,
    ".",
    plot_params$device
  )

  a <- VennDiagram::venn.diagram(
    l,
    category.names = category_names,
    filename = out_filename,
    output = TRUE,

    # Output features
    imagetype = plot_params$device,
    units = plot_params$units,
    height = plot_params$height,
    width = plot_params$width,
    resolution = plot_params$dpi,
    compression = "lzw",

    # Circles
    lwd = 2,
    lty = "blank",
    fill = palette,
    alpha = 0.4,

    # Numbers
    cex = font_size,
    fontface = "bold",
    fontfamily = "sans",
    print.mode = c("raw", "percent"),
    disable.logging = FALSE,
    ...
  )

  if (nchar(plot_params$filename) == 0) grid.draw(a)

  return(NULL)
}
