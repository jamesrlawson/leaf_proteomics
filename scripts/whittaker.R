require(ggplot2)
ggplot(climate_locs, aes(x = prec, y = tavg)) + geom_point() + theme_classic() +
  scale_y_reverse(lim=c(30,-15)) + xlim(0,4500) + theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )
ggsave("output/whitaker_transparent.pdf", device = 'pdf', bg = "transparent")
