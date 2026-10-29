# Observed and fitted proportion divorced, widowed or in a later union on the
# three post-marriage margins, one figure each.
#
#   Rscript fig_post_band.R            (run from paper1/figcode/)
#
# Writes ../figures/post_band_dur.png, post_band_yom.png and post_band_afm.png,
# three outcomes by row and four countries by column.
#
# Input, one per country, in data/:
#   post_band_<CC>.csv   the adopted additive surface, columns
#     by     the margin, "dur" completed years since first marriage in two-year
#            cells, "yom" year of first marriage in five-year cohorts, "afm" age
#            at first marriage in single years
#     key    the bin value on that margin
#     state  M, D1, W1, R, D2, W2 and the grouped outcomes divorced, widowed and
#            remarried, which are the three rows drawn here
#     N      weighted number in the cell, rows the unweighted count
#     obs    observed proportion
#     fit    fitted proportion, lo and hi its 95 per cent interval over ndraw
#            draws of the posterior
#
# The interval on an observed point is the Jeffreys one, the 2.5 and 97.5 per
# cent points of Beta(N q + 1/2, N (1 - q) + 1/2) at the observed proportion q
# and the weighted number N in the cell. It is not stored in the input. On the
# year-of-marriage margin it reproduces the interval of the earlier
# data/fit_vs_obs_by_marriage_year.csv to five decimal places, so the points
# of these figures carry the interval the published ones carried.
#
# A cell holding under 2 per cent of its country's weight on its own margin is
# omitted, which is the rule the manuscript captions state. On the duration
# margin it keeps 59 cells, the 177 the Results report over three outcomes.
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
source("theme_kinh.R")

cc_lab <- c(RW = "Rwanda", MW = "Malawi", SL = "Sierra Leone", ST = "Sao Tome")
DAT_DIR <- Sys.getenv("DAT_DIR", "data")
OUT_DIR <- Sys.getenv("OUT_DIR", "../figures")
SHARE_MIN <- 0.02

need <- sprintf("post_band_%s.csv", names(cc_lab))
miss <- need[!file.exists(file.path(DAT_DIR, need))]
if (length(miss))
  stop("missing input, expected in ", DAT_DIR, "/:\n  ",
       paste(miss, collapse = "\n  "), call. = FALSE)

d <- bind_rows(lapply(names(cc_lab), function(k) {
  x <- read.csv(file.path(DAT_DIR, sprintf("post_band_%s.csv", k)))
  # the weight of a cell is its share of the country's weight on its own margin,
  # read off the distinct cells so that the nine states are not counted nine times
  tot <- x %>% distinct(by, key, N) %>% group_by(by) %>%
    summarise(Ntot = sum(N), .groups = "drop")
  x %>% left_join(tot, by = "by") %>% mutate(cc = k, share = N / Ntot)
})) %>%
  filter(state %in% c("divorced", "widowed", "remarried"), share >= SHARE_MIN) %>%
  mutate(country = factor(cc_lab[cc], levels = cc_lab),
         row = factor(state, levels = c("divorced", "widowed", "remarried"),
                      labels = c("Divorced", "Widowed", "In a later union")),
         obs_lo = qbeta(0.025, N * obs + 0.5, N * (1 - obs) + 0.5),
         obs_hi = qbeta(0.975, N * obs + 0.5, N * (1 - obs) + 0.5))

margin <- function(dat, file, xlab, xbreaks = waiver()) {
  dat <- dat %>% group_by(country, row) %>% mutate(w = N / max(N)) %>%
    ungroup() %>%
    group_by(row) %>%
    mutate(ymax = 1.08 * max(c(fit, hi, obs), na.rm = TRUE)) %>% ungroup() %>%
    mutate(hi_cut = pmin(obs_hi, ymax), lo_cut = pmin(obs_lo, ymax),
           over = obs > ymax, obs_cut = pmin(obs, ymax),
           band_hi_cut = pmin(hi, ymax))
  blank <- dat %>% distinct(country, row, ymax) %>%
    mutate(x = min(dat$key), y = ymax)
  span <- diff(range(dat$key))
  tag <- dat %>% distinct(country, row, ymax) %>% arrange(row, country) %>%
    mutate(lab = LETTERS[seq_len(n())], x = min(dat$key) + 0.015 * span,
           y = 0.985 * ymax)
  p <- ggplot() +
    geom_blank(data = blank, aes(x, y)) +
    geom_ribbon(data = dat, aes(key, ymin = lo, ymax = band_hi_cut,
                                fill = row), alpha = 0.3) +
    geom_linerange(data = dat, aes(key, ymin = lo_cut, ymax = hi_cut,
                                   colour = row), linewidth = LW_RANGE) +
    geom_point(data = filter(dat, !over), aes(key, obs, size = w, colour = row),
               shape = 16) +
    geom_point(data = filter(dat, over), aes(key, obs_cut, size = w,
                                             colour = row),
               shape = 2, stroke = 0.4) +
    geom_line(data = dat, aes(key, fit, colour = row), linewidth = LW_LINE) +
    geom_text(data = tag, aes(x, y, label = lab), family = default_font,
              size = TAG_PT / .pt, fontface = "bold", hjust = 0, vjust = 1,
              colour = "grey20") +
    scale_size_continuous(range = c(0.3, 2.6), guide = "none") +
    # one colour per outcome, the same on the three margins, so the divorced row
    # of this figure and of the two appendix figures read alike; the strip label
    # names the row, so no legend is drawn
    scale_colour_manual(values = setNames(ser_pal, levels(dat$row)),
                        guide = "none") +
    scale_fill_manual(values = setNames(ser_pal, levels(dat$row)),
                      guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
    scale_x_continuous(breaks = xbreaks, expand = expansion(mult = 0.03)) +
    facet_grid(row ~ country, scales = "free_y") +
    labs(x = xlab, y = "Share of women") + k.theme
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(OUT_DIR, file), plot = p, device = ragg::agg_png,
         width = W_FULL, height = 5.6, dpi = DPI, bg = "white")
  cat("wrote", file.path(OUT_DIR, file), "\n")
}

# Full journal width: four countries side by side leave each panel under 2 cm at
# one column. The y scale is common to the four countries of a row and free
# between the three rows, since the later-union proportion reaches three times
# the widowed one and one scale over twelve panels flattens the widowed row.
margin(filter(d, by == "dur"), "post_band_dur.png", "Years since first marriage")
# two breaks on the year of marriage: four-digit labels at 10 pt collide at
# three or more across a panel 3.5 cm wide
margin(filter(d, by == "yom"), "post_band_yom.png", "Year of first marriage",
       c(1980, 2000))
margin(filter(d, by == "afm"), "post_band_afm.png", "Age at first marriage")
