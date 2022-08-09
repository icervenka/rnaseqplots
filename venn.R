plot_venn2 = function(x, y, names, savename) {
  VennDiagram::venn.diagram(
    x = list(tolower(x), tolower(y)),
    category.names = names,
    filename = savename,
    output=TRUE,

    # Output features
    imagetype="tiff" ,
    units = 'mm',
    height = 35 ,
    width = 35 ,
    resolution = 600,
    compression = "lzw",

    # Circles
    lwd = 2,
    lty = 'blank',
    fill = viridis(2),

    # Numbers
    cex = .6,
    fontface = "bold",
    fontfamily = "sans",

    # Set names
    cat.cex = 0.6,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-25, 25),
    cat.dist = c(0.055, 0.055),
    cat.fontfamily = "sans"
  )
}
