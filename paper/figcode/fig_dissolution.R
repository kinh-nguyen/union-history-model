# Dissolution of a first marriage against the share of women out of a union at
# the same duration, with the 95 per cent interval of each.
#
#   Rscript fig_dissolution.R          (run from paper1/figcode/)
#
# Writes ../figures/dissolution.png.
#
# Input, one per country, in data/:
#   post_bandx_<CC>.csv   the adopted additive surface, the columns of
#                         post_band_<CC>.csv and two further states,
#     dissolved     the share whose first marriage has ended, one minus the
#                   share still in it, so it counts the divorced, the widowed
#                   and those who have gone on to a later union
#     unpartnered   the share out of a union at that duration, the divorced and
#                   the widowed from either union
# Both are summed inside each draw before its quantile is taken, which a sum of
# the stored quantiles of the separate states would not reproduce.
#
# The two together are the statement the introduction cites from the
# literature, that dissolution of a first union is frequent while repartnering
# keeps the number out of a union at any one time small: the gap between the two
# curves is the share who have left a first marriage and are in a later union.
#
# The interval on an observed point is the Jeffreys one and a cell holding under
# 2 per cent of its country's weight on the duration margin is omitted, as in
# fig_post_band.R.
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
source("theme_kinh.R")

cc_lab <- c(RW = "Rwanda", MW = "Malawi", SL = "Sierra Leone", ST = "Sao Tome")
DAT_DIR <- Sys.getenv("DAT_DIR", "data")
OUT_DIR <- Sys.getenv("OUT_DIR", "../figures")
SHARE_MIN <- 0.02

need <- sprintf("post_bandx_%s.csv", names(cc_lab))
miss <- need[!file.exists(file.path(DAT_DIR, need))]
if (length(miss))
  stop("missing input, expected in ", DAT_DIR, "/:\n  ",
       paste(miss, collapse = "\n  "),
       "\nRegenerate with, from msm_post/tmb_par: ",
       "Rscript predict_bandx_cc.R <CC> \"<adopted options>\" 400",
       call. = FALSE)

d <- bind_rows(lapply(names(cc_lab), function(k) {
  x <- read.csv(file.path(DAT_DIR, sprintf("post_bandx_%s.csv", k)))
  tot <- x %>% distinct(by, key, N) %>% group_by(by) %>%
    summarise(Ntot = sum(N), .groups = "drop")
  x %>% left_join(tot, by = "by") %>% mutate(cc = k, share = N / Ntot)
})) %>%
  filter(by == "dur", state %in% c("dissolved", "unpartnered"),
         share >= SHARE_MIN) %>%
  mutate(country = factor(cc_lab[cc], levels = cc_lab),
         series = factor(state, levels = c("dissolved", "unpartnered"),
                         labels = c("First marriage has ended",
                                    "Out of a union")),
         obs_lo = qbeta(0.025, N * obs + 0.5, N * (1 - obs) + 0.5),
         obs_hi = qbeta(0.975, N * obs + 0.5, N * (1 - obs) + 0.5))

ser_col <- c("First marriage has ended" = unname(okcol["orange"]),
             "Out of a union" = unname(okcol["blue"]))

ymax <- 1.08 * max(c(d$fit, d$hi, d$obs), na.rm = TRUE)
d <- d %>% mutate(w = N / max(N),
                  hi_cut = pmin(obs_hi, ymax), lo_cut = pmin(obs_lo, ymax),
                  band_hi_cut = pmin(hi, ymax))
span <- diff(range(d$key))
tag <- d %>% distinct(country) %>% arrange(country) %>%
  mutate(lab = LETTERS[seq_len(n())], x = min(d$key) + 0.02 * span,
         y = 0.985 * ymax)

p <- ggplot(d, aes(key, colour = series, fill = series)) +
  geom_blank(data = tag, aes(x, y), inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = lo, ymax = band_hi_cut), alpha = 0.3, colour = NA) +
  geom_linerange(aes(ymin = lo_cut, ymax = hi_cut), linewidth = LW_RANGE,
                 alpha = 0.75) +
  geom_point(aes(y = obs, size = w), shape = 16) +
  geom_line(aes(y = fit), linewidth = LW_LINE) +
  geom_text(data = tag, aes(x, y, label = lab), inherit.aes = FALSE,
            family = default_font, size = TAG_PT / .pt, fontface = "bold",
            hjust = 0, vjust = 1, colour = "grey20") +
  facet_wrap(~ country, nrow = 1) +
  scale_colour_manual(values = ser_col, name = NULL) +
  scale_fill_manual(values = ser_col, name = NULL) +
  scale_size_continuous(range = c(0.3, 2.2), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
  scale_x_continuous(expand = expansion(mult = 0.03)) +
  labs(x = "Years since first marriage", y = "Share of women") +
  k.theme +
  theme(legend.position = "bottom", legend.margin = margin(t = -2, b = 0),
        legend.text = element_text(size = BASE_PT))

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(OUT_DIR, "dissolution.png"), plot = p, device = ragg::agg_png,
       width = W_FULL, height = 3.2, dpi = DPI, bg = "white")
cat("wrote", file.path(OUT_DIR, "dissolution.png"), "\n")

## the values the appendix quotes, the longest duration kept per country
q <- d %>% group_by(country, state) %>% filter(key == max(key)) %>%
  transmute(key, fit = round(fit, 3), lo = round(lo, 3), hi = round(hi, 3))
print(as.data.frame(q), row.names = FALSE)
