plot_venn <- function(l,
                      category_names = "",
                      palette = viridis::viridis(2),
                      width = 35,
                      height = 35,
                      units = "mm",
                      font_size = 0.6,
                      savename = NULL,
                      ...) {
  if (is.null(savename)) grid.newpage()

  a <- VennDiagram::venn.diagram(
    l,
    category.names = category_names,
    filename = savename,
    output = TRUE,

    # Output features
    imagetype = "tiff",
    units = units,
    height = height,
    width = width,
    resolution = 600,
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

    # Set names
    # cat.cex = 0.6,
    # cat.fontface = "bold",
    # cat.default.pos = "outer",
    # cat.pos = c(-25, 25),
    # cat.dist = c(0.055, 0.055),
    # cat.fontfamily = "sans",

    print.mode = c("raw", "percent"),
    disable.logging = FALSE,
    ...
  )

  if (is.null(savename)) grid.draw(a)
}
