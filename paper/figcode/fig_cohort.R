# The cohort trends the model produces, with their 95 per cent interval.
#
#   Rscript fig_cohort.R               (run from paper1/figcode/)
#
# Writes ../figures/cohort_trend.png.
#
# Input, in data/:
#   cohort_band_all.csv   written by ROOT/cohort/band_cohort.R, one row per
#                         country, birth cohort and quantity
#
# The three quantities are the median age at first sex, the median age at first
# marriage, and the interval between the two medians. Each is the age at which a
# probability of the fitted chain crosses one half, so each is a function of the
# fixed effects and carries the delta-method interval of that function, the same
# linearisation the other figures of the paper use.
#
# Cohorts are restricted to those a survey observed to age 30, since the cohort
# term is a linear shift of the peak and a younger cohort's median age at first
# marriage would be read off an extrapolation of that shift rather than off
# observed marriages. The restriction is stated in the caption.
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr) })
source("theme_kinh.R")

cc_lab  <- c(RW = "Rwanda", MW = "Malawi", SL = "Sierra Leone", ST = "Sao Tome")
DAT_DIR <- Sys.getenv("DAT_DIR", "data")
OUT_DIR <- Sys.getenv("OUT_DIR", "../figures")
AAI_MIN <- as.integer(Sys.getenv("AAI_MIN", "30"))

f <- file.path(DAT_DIR, "cohort_band_all.csv")
if (!file.exists(f))
  stop("missing input ", f, "\nRegenerate with, from the msm tree: ",
       "Rscript cohort/band_cohort.R <CC> for CC in RW MW SL ST",
       call. = FALSE)

d <- read.csv(f) %>%
  filter(max_aai >= AAI_MIN) %>%
  mutate(country = factor(cc_lab[cc], levels = cc_lab),
         row = factor(quantity, levels = c("afs", "afm", "gap"),
                      labels = c("Median age (years)", "Median age (years)",
                                 "Interval (years)")),
         series = factor(quantity, levels = c("afs", "afm", "gap"),
                         labels = c("First sex", "First marriage",
                                    "Debut to marriage")))

ser_col <- c("First sex" = unname(okcol["blue"]),
             "First marriage" = unname(okcol["orange"]),
             "Debut to marriage" = "gray1")

## part labels A to H, reading across the rows
span <- diff(range(d$yob))
tag <- d %>% distinct(country, row) %>% arrange(row, country) %>%
  mutate(lab = LETTERS[seq_len(n())])
rng <- d %>% group_by(row) %>%
  summarise(ymin = min(lo), ymax = max(hi), .groups = "drop")
tag <- left_join(tag, rng, by = "row") %>%
  mutate(x = min(d$yob) + 0.015 * span, y = ymax + 0.055 * (ymax - ymin))

p <- ggplot(d, aes(yob, fit, colour = series, fill = series)) +
  geom_blank(data = tag, aes(x, y), inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.35, colour = NA) +
  geom_line(linewidth = LW_LINE) +
  geom_text(data = tag, aes(x, y, label = lab), inherit.aes = FALSE,
            family = default_font, size = TAG_PT / .pt, fontface = "bold",
            hjust = 0, vjust = 1, colour = "grey20") +
  facet_grid(row ~ country, scales = "free_y", switch = "y") +
  scale_colour_manual(values = ser_col, name = NULL) +
  scale_fill_manual(values = ser_col, name = NULL) +
  ## interior decades only: a break at either end of a panel collides with the
  ## break at the end of the panel beside it at this width
  scale_x_continuous(breaks = seq(1960, 1980, 10),
                     expand = expansion(mult = 0.02)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = "Year of birth", y = NULL) +
  k.theme +
  theme(strip.placement = "outside",
        strip.text.y.left = element_text(angle = 90, hjust = 0.5),
        legend.position = "bottom",
        legend.margin = margin(t = -2, b = 0),
        legend.text = element_text(size = BASE_PT))

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(OUT_DIR, "cohort_trend.png"), plot = p, device = ragg::agg_png,
       width = W_FULL, height = 4.1, dpi = DPI, bg = "white")
cat("wrote", file.path(OUT_DIR, "cohort_trend.png"), "\n")

## the values the appendix quotes, oldest and youngest cohort kept per country
q <- d %>% group_by(country, quantity) %>%
  filter(yob == min(yob) | yob == max(yob)) %>%
  arrange(country, quantity, yob) %>%
  transmute(yob, fit = round(fit, 2), lo = round(lo, 2), hi = round(hi, 2))
print(as.data.frame(q), row.names = FALSE)
