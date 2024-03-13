afr = image_read_pdf("AFR.G6PD.pdf")
e6a = ggplot() +
  annotation_raster(afr, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_blank() + theme_nothing()
meta = image_read_pdf("META_raw.G6PD.pdf")
e6b = ggplot() +
  annotation_raster(meta, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_blank() + theme_nothing()


extended_data_fig_8 = function(output_format = 'png') {
  output_type(output_format, paste0('extended_data_figure8.', output_format), height=5, width=10)
  print(ggarrange(e6a, e6b, ncol = 2, labels='auto'))
  dev.off()
}

extended_data_fig_8()
