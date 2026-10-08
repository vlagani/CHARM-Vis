# Common settings for the figures: sizes in mm, font sizes in pt

# widths (A4 page, one row of a figure spans the full width)
full_width <- 180
half_width <- 88

# fonts
font_family <- 'Helvetica'
font_title <- 7 # axis and panel titles
font_text <- 6  # tick labels, legends, labels inside the plots
font_small <- 5 # minimum, for dense labels

# cell types: names and colors shared across panels
cell_type_colors <- c(`Glutamatergic neurons` = '#F8766D',
                      `Gabaergic neurons` = '#B79F00',
                      Astrocytes = '#00BA38',
                      Ependymal = '#01BFC4',
                      OPC = '#619CFF',
                      Other = '#F564E3')

# ggplot2 theme with the font sizes above
theme_figure <- function(){
  ggplot2::theme(text = ggplot2::element_text(family = font_family, size = font_text),
                 axis.title = ggplot2::element_text(size = font_title),
                 axis.text = ggplot2::element_text(size = font_text),
                 plot.title = ggplot2::element_text(size = font_title, face = 'bold'),
                 legend.title = ggplot2::element_text(size = font_text),
                 legend.text = ggplot2::element_text(size = font_text))
}

# saving a panel as pdf and png (600 dpi)
# file: path without extension; draw: function drawing the panel;
# width, height: size in mm; pointsize: base font size of the device (pt)
save_panel <- function(file, draw, width, height, pointsize = font_text){
  grDevices::cairo_pdf(paste0(file, '.pdf'), width = width / 25.4,
                       height = height / 25.4, family = font_family,
                       pointsize = pointsize)
  draw()
  grDevices::dev.off()
  grDevices::png(paste0(file, '.png'), width = width, height = height,
                 units = 'mm', res = 600, type = 'cairo',
                 family = font_family, pointsize = pointsize)
  draw()
  grDevices::dev.off()
}
