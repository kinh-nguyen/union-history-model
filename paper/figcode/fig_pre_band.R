# The adopted pre-marriage model with its 95 per cent interval, two figures.
#
#   Rscript fig_pre_band.R             (run from paper1/figcode/)
#
# Writes ../figures/pre_band_age.png and ../figures/pre_band_dur.png.
#
# Inputs, one pair per country, in data/:
#   fit5_prob_<CC>_ph_c2.csv   written by ROOT/code/fit_ph.R   with STEM=ph_c2
#   fit5_band_<CC>_ph_c2.csv   written by ROOT/code/band_ph.R  with STEM=ph_c2
# The eight files were staged into data/ from the cluster on 2026-10-27. The
# script stops with the list of the missing ones if any of them is removed.
#
# One model only, since the comparison belongs in the appendix. The interval of
# the fitted curve is read from the band file, whose method column records it as
# the delta-method one. Point area and the bars on the panel floor both give the
# number of women at risk in the cell, the bars on a square-root scale; the y
# limit of each row stops just above the largest fitted value, interval and
# well-populated observed value of that row; an observed value above that limit
# is drawn as an open triangle on the edge.
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
source("theme_kinh.R")

cc_lab <- c(RW = "Rwanda", MW = "Malawi", SL = "Sierra Leone", ST = "Sao Tome")
STEM    <- Sys.getenv("STEM", "ph_c2")
DAT_DIR <- Sys.getenv("DAT_DIR", "data")        # optional override
OUT_DIR <- Sys.getenv("OUT_DIR", "../figures")  # optional override
AGE_MAX <- 40; DUR_MAX <- 30

need <- c(sprintf("fit5_prob_%s_%s.csv", names(cc_lab), STEM),
          sprintf("fit5_band_%s_%s.csv", names(cc_lab), STEM))
miss <- need[!file.exists(file.path(DAT_DIR, need))]
if (length(miss))
  stop("missing input, expected in ", DAT_DIR, "/:\n  ",
       paste(miss, collapse = "\n  "),
       "\nRegenerate with, from ROOT: STEM=", STEM,
       " Rscript code/fit_ph.R <CC> && STEM=", STEM,
       " Rscript code/band_ph.R <CC>", call. = FALSE)

d <- bind_rows(lapply(names(cc_lab), function(k) {
  p <- read.csv(file.path(DAT_DIR, sprintf("fit5_prob_%s_%s.csv", k, STEM)))
  b <- read.csv(file.path(DAT_DIR, sprintf("fit5_band_%s_%s.csv", k, STEM)))
  left_join(p[p$N > 0, ], b[, c("x", "transition", "axis", "band_lo", "band_hi")],
            by = c("x", "transition", "axis"))
})) %>%
  mutate(country = factor(cc_lab[cc], levels = cc_lab),
         row = factor(transition, levels = c("V->X", "V->M", "X->M"),
                      labels = c("Sexual debut\n(V to X)",
                                 "Marriage in\nthe debut year",
                                 "Marriage after\ndebut (X to M)"))) %>%
  filter((axis == "age" & x >= 12 & x <= AGE_MAX) |
         (axis == "dur" & x >= 0  & x <= DUR_MAX))

pane <- function(dat, by_row, file, width, height, xlab) {
  grp <- if (by_row) "row" else "axis"
  dat <- dat %>% group_by(country, row) %>% mutate(w = N / max(N)) %>% ungroup() %>%
    mutate(g = .data[[grp]])
  ylim <- bind_rows(transmute(dat, g, v = fitted_prob),
                    transmute(dat, g, v = band_hi),
                    transmute(filter(dat, w >= 0.01), g, v = obs)) %>%
    group_by(g) %>% summarise(ymax = 1.08 * max(v, na.rm = TRUE), .groups = "drop")
  dat <- left_join(dat, ylim, by = "g") %>%
    mutate(hi_cut = pmin(hi, ymax), lo_cut = pmin(lo, ymax),
           over = obs > ymax, obs_cut = pmin(obs, ymax),
           bar = 0.085 * ymax * sqrt(w), band_hi_cut = pmin(band_hi, ymax))
  blank <- dat %>% distinct(country, row, ymax) %>% mutate(x = min(dat$x), y = ymax)
  # part labels A, B, C ... in the upper left of every panel, reading across the
  # rows, with the note on the bars on the line under the first of them
  span <- diff(range(dat$x))
  tag <- dat %>% distinct(country, row, ymax) %>%
    arrange(row, country) %>%
    mutate(lab = LETTERS[seq_len(n())], x = min(dat$x) + 0.015 * span,
           y = 0.985 * ymax)
  note <- tag %>% slice(1) %>%                 # second line of the first panel
    mutate(y = 0.86 * ymax, lab = "bars: women at risk")
  p <- ggplot() +
    geom_blank(data = blank, aes(x, y)) +
    geom_col(data = dat, aes(x, bar), width = 0.9, fill = okcol["base1"],
             alpha = 0.45, colour = NA) +
    geom_ribbon(data = dat, aes(x, ymin = band_lo, ymax = band_hi_cut,
                                fill = row), alpha = 0.3) +
    geom_linerange(data = dat, aes(x, ymin = lo_cut, ymax = hi_cut,
                                   colour = row), linewidth = LW_RANGE) +
    geom_point(data = filter(dat, !over), aes(x, obs, size = w, colour = row),
               shape = 16) +
    geom_point(data = filter(dat, over), aes(x, obs_cut, size = w, colour = row),
               shape = 2, stroke = 0.4) +
    geom_line(data = dat, aes(x, fitted_prob, colour = row),
              linewidth = LW_LINE) +
    geom_text(data = note, aes(x, y, label = lab), family = default_font,
              size = AXIS_PT / .pt, hjust = 0, vjust = 1,
              colour = okcol["base00"]) +
    geom_text(data = tag, aes(x, y, label = lab), family = default_font,
              size = TAG_PT / .pt, fontface = "bold", hjust = 0, vjust = 1,
              colour = "grey20") +
    scale_size_continuous(range = c(0.3, 2.6), guide = "none") +
    # one colour per transition, held across the two figures, so the marriage
    # intensity after debut is the same colour on the age axis and on the
    # duration axis; the strip label names the row, so no legend is drawn
    scale_colour_manual(values = setNames(ser_pal, levels(dat$row)),
                        guide = "none") +
    scale_fill_manual(values = setNames(ser_pal, levels(dat$row)),
                      guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
    scale_x_continuous(expand = expansion(mult = 0.02)) +
    labs(x = xlab, y = "Annual probability") + k.theme
  p <- if (by_row) p + facet_grid(row ~ country, scales = "free_y") else
                   p + facet_wrap(~ country, nrow = 1)
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(OUT_DIR, file), plot = p, device = ragg::agg_png,
         width = width, height = height, dpi = DPI, bg = "white")
  cat("wrote", file.path(OUT_DIR, file), "\n")
}

# Full journal width for both. Four countries side by side at 8.5 cm would leave
# each panel under 2 cm, which no amount of type size rescues.
# The y scale is common to the four countries of a row and free between rows:
# the three transitions differ in level by a factor of two, and one scale over
# all twelve panels flattens the first row.
pane(filter(d, axis == "age"), TRUE, "pre_band_age.png", W_FULL, 5.6,
     "Age (years)")
pane(filter(d, axis == "dur"), FALSE, "pre_band_dur.png", W_FULL, 2.6,
     "Years since sexual debut")
