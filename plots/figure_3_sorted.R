#!/usr/bin/env Rscript

library(ggplot2)

data <- read.csv2("results/data/multiple_hypothesis_testing_fig3_sorted.csv")

ggplot(
  data,
  aes(
    x = hypothesis,
    y = e_value,
    colour = method,
    fill = method,
    linetype = bound
  )
) +
  stat_summary(
    fun = median,
    geom = "line"
  ) +
  stat_summary(
    geom = "ribbon",
    fun.min = \(x) quantile(x, 0.25),
    fun.max = \(x) quantile(x, 0.75),
    alpha = 0.25,
    colour = NA
  ) +
  scale_y_continuous(
    trans = "log10"
  ) +
  xlab("sorted e-value index")
ggsave(
  "results/plots/figure_3_sorted.png",
  width = 11,
  height = 8,
  units = "in",
  device = "png"
)