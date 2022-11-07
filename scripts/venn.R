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
}
