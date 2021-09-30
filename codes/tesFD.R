p <- ggplot(carb.df,
            aes(x = `Top taxa`,y = abundance,
                fill=Group)) +
  #geom_violin() +
  geom_jitter(width = 0.2, alpha=0.25, shape=21)+
  #geom_jitter() +
  #facet_grid(~Group, scales = "free") +
  theme_biome_utils() +
  theme(axis.title.x = ggplot2::element_blank(),
        #axis.text.x = ggplot2::element_text(angle=90),
        #axis.text.x = ggplot2::element_blank(),
        #axis.ticks.x = ggplot2::element_blank(),
        strip.background = ggplot2::element_blank(),
        strip.text = element_text(size = 10),
        #panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.spacing = unit(0, "line"),
        legend.position = "none") +
  scale_fill_manual(values = cols_grps) +
  ylab("Median Reads per kilobase") + xlab("") +
  labs(subtitle = "Carbohydrate degradation")
p <- p + coord_flip()
p


data <- otu.tb
core.which <- function(data, intTr, prevalenceTr) {
  d.bin <- data >= intTr
  prevalences <- rowSums(d.bin)
  nOTUs <- as.numeric(prevalences >= prevalenceTr)
  return(nOTUs)
}

core.which(otu.tb, intTr = 0.00001, prevalenceTr = 0.7)
