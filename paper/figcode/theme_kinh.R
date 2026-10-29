library(ggplot2); library(ragg)

default_font <- "sans"

okcol <- c(base03 = "#002d38", base02 = "#093946", base01 = "#5b7279",
           base00 = "#657377", base0 = "#98a8a8", base1 = "#8faaab",
           base2 = "#f1e9d2", base3 = "#fbf7ef", yellow = "#ac8300",
           orange = "#d56500", red = "#f23749", blue = "#2b90d8")

BASE_PT  <- 11    # axis titles, legend text
AXIS_PT  <- 10    # tick labels
STRIP_PT <- 10.5  # panel strip labels
TAG_PT   <- 12    # the A, B, C part labels
LW_PANEL <- 0.3   # panel border
LW_LINE  <- 0.55  # fitted curves
LW_RANGE <- 0.35  # observed intervals

ser_pal <- c(unname(okcol["blue"]), unname(okcol["orange"]),
             unname(okcol["base03"]))

k.theme <- theme_light(base_family = default_font, base_size = BASE_PT) +
  theme(text = element_text(colour = okcol["base00"]),
        axis.text = element_text(size = AXIS_PT, colour = okcol["base00"]),
        axis.title = element_text(size = BASE_PT, colour = okcol["base00"]),
        panel.grid = element_blank(),
        panel.spacing = unit(9, "pt"),
        panel.background = element_rect(fill = okcol["base3"], colour = NA),
        panel.border = element_rect(linewidth = LW_PANEL, colour = "gray1", fill = NA),
        strip.background = element_rect(fill = "white", colour = NA),
        strip.text = element_text(colour = "grey20", hjust = 0, size = STRIP_PT),
        legend.key = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.title = element_text(colour = okcol["yellow"]),
        plot.caption = element_text(colour = okcol["yellow"]),
        plot.title = element_blank(),            # no title inside the image
        plot.tag = element_text(size = TAG_PT, face = "bold", colour = "grey20",
                                hjust = 0, vjust = 1),
        plot.tag.position = c(0.005, 0.995))     # upper left corner of the part

# Journal widths, in inches, for ggsave.
W_ONE  <- 8.5 / 2.54
W_FULL <- 17.5 / 2.54
DPI    <- 400
