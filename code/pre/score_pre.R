# Scores for the pre-marriage comparison table.
#
#   Rscript code/score_pre.R
#
# For every fit of the phase template, the unpenalised log likelihood is
# recomputed with the priors switched off, and the cell deviance is summed on
# each axis of the figure over cells holding at least a hundredth of the women at
# risk of their panel. Writes fitted/summary/pre_table.csv.
suppressPackageStartupMessages({ library(TMB); library(here); library(dplyr) })
DLL <- "msm_ph"
dyn.load(TMB::dynlib(here("code", DLL)))
invisible(TMB::config(tape.parallel = 0, DLL = DLL))
TMB::openmp(1)
CCS <- c("RW", "MW", "SL", "ST")
STEMS <- strsplit(Sys.getenv("STEMS", "ph_off,ph_c2"), ",")[[1]]

cc_all <- vroom::vroom(here("data/cc8.csv.bz2"), show_col_types = FALSE)
out <- list()
for (CC in CCS) {
  raw <- cc_all %>% filter(cc == CC, sex == "f") %>%
    mutate(afs = if_else(is.na(afs), 0, afs), afm = if_else(is.na(afm), 0, afm),
           coh = (yob - 1975) / 10)
  wom <- raw %>% filter(start == 0) %>%
    transmute(afs = as.integer(afs), afm = as.integer(afm),
              aai = as.integer(aai), coh, n)
  cohs <- sort(unique(wom$coh)); A_MAX <- 60L
  mk <- function(nopen) list(afs = wom$afs, afm = wom$afm,
    aai = pmin(wom$aai, A_MAX), ic = match(wom$coh, cohs) - 1L, coh_u = cohs,
    a_max = A_MAX, n = wom$n, nopen = as.integer(nopen))
  for (ST in STEMS) {
    f <- here("fit", sprintf("fit5_%s_%s.rds", CC, ST))
    if (!file.exists(f)) { cat("missing", CC, ST, "\n"); next }
    o <- readRDS(f)
    mp <- if (length(o$map)) o$map else NULL
    obj_u <- TMB::MakeADFun(mk(1L), o$pl, map = mp, DLL = DLL, silent = TRUE)
    nll_u <- obj_u$fn(o$par)
    npar <- length(o$par)
    p <- read.csv(here("fitted", "prob", sprintf("fit5_prob_%s_%s.csv", CC, ST)))
    p <- p %>% group_by(transition, axis) %>% mutate(w = N / max(N)) %>% ungroup() %>%
      filter(w >= 0.01, N > 0)
    dev1 <- function(d) {
      y <- d$events; n <- d$N; pr <- pmin(pmax(d$fitted_prob, 1e-12), 1 - 1e-12)
      2 * sum(ifelse(y > 0, y * log(y / (n * pr)), 0) +
              ifelse(n - y > 0, (n - y) * log((n - y) / (n * (1 - pr))), 0))
    }
    out[[paste(CC, ST)]] <- data.frame(cc = CC, spec = ST, npar = npar,
      nll_unpen = nll_u, aic_unpen = 2 * nll_u + 2 * npar,
      dev_age = dev1(filter(p, axis == "age")),
      dev_dur = dev1(filter(p, axis == "dur")))
    cat(sprintf("%s %-7s npar %d AIC %.1f dev age %.1f dur %.1f\n", CC, ST, npar,
                2 * nll_u + 2 * npar, out[[paste(CC, ST)]]$dev_age,
                out[[paste(CC, ST)]]$dev_dur))
  }
}
tab <- bind_rows(out)
write.csv(tab, here("fitted", "summary", "pre_table.csv"), row.names = FALSE)
agg <- tab %>% group_by(spec) %>%
  summarise(npar = sum(npar), aic = sum(aic_unpen), dev_age = sum(dev_age),
            dev_dur = sum(dev_dur), .groups = "drop")
print(as.data.frame(agg))
write.csv(agg, here("fitted", "summary", "pre_table_total.csv"), row.names = FALSE)
cat("done\n")
