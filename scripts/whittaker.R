require(ggplot2)
ggplot(climate_locs, aes(x = prec, y = tavg)) + geom_point() + theme_classic() +
  scale_y_reverse(lim=c(30,-15)) + xlim(0,4500) + theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )
ggsave("output/whitaker_transparent.pdf", device = 'pdf', bg = "transparent")


require(ggplot2)
ggplot(data, aes(x = tavg, y = log10prec)) + geom_point(size = 2) + theme_classic() +
  scale_x_continuous(breaks=c(10,15,20,25)) +
  scale_y_continuous(breaks=c(-0.6,-0.4,-0.2,0,0.2,0.4)) +
# xlim(-15,30) + ylim(-0.8,0.6) + 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )
ggsave("output/clim_transparent.pdf", device = 'pdf', bg = "transparent")
