# Delta-method interval for the fitted curves of the phase model.
#
#   Rscript code/band_delta.R <CC>
#
# The reweighted normal proposal was screened on 2026-10-22 and rejected again
# (effective sample size 5.9 of 1500 proposals in Sao Tome), so the interval here
# is the linearisation at the mode, stated as such: the Jacobian of each one-year
# cell probability with respect to the fixed effects is taken by central
# differences, the variance is g' V g with V the covariance from sdreport, and the
# interval is formed on the log-odds scale so it stays inside zero and one.
#
# Writes fitted/prob/fit5_band_<CC>_<STEM>.csv.
suppressPackageStartupMessages({ library(TMB); library(here); library(dplyr) })
source(here("code", "ph_cells.R"))

CC <- commandArgs(trailingOnly = TRUE)[1]; if (is.na(CC)) CC <- "MW"
STEM <- Sys.getenv("STEM", "ph_c2")
DLL <- "msm_ph"
cat("country", CC, "stem", STEM, "\n")

cc <- vroom::vroom(here("data/cc8.csv.bz2"), show_col_types = FALSE)
raw <- cc %>% filter(cc == CC, sex == "f") %>%
  mutate(afs = if_else(is.na(afs), 0, afs), afm = if_else(is.na(afm), 0, afm),
         coh = (yob - 1975) / 10)
wom <- raw %>% filter(start == 0) %>%
  transmute(afs = as.integer(afs), afm = as.integer(afm), aai = as.integer(aai),
            coh, n)
cohs <- sort(unique(wom$coh)); A_MAX <- 60L; AGES <- 1:A_MAX
dat <- list(afs = wom$afs, afm = wom$afm, aai = pmin(wom$aai, A_MAX),
            ic = match(wom$coh, cohs) - 1L, coh_u = cohs, a_max = A_MAX,
            n = wom$n, nopen = 0L)

o <- readRDS(here("fit", sprintf("fit5_%s_%s.rds", CC, STEM)))
pl <- o$pl
dyn.load(TMB::dynlib(here("code", DLL)))
invisible(TMB::config(tape.parallel = 0, DLL = DLL))
TMB::openmp(1)
obj <- TMB::MakeADFun(dat, pl, DLL = DLL, silent = TRUE)
stopifnot(max(abs(obj$par - o$par)) < 1e-8)
sdr <- TMB::sdreport(obj, par.fixed = o$par)
V <- sdr$cov.fixed; V <- (V + t(V)) / 2
cat("parameters", length(o$par), "largest sd", round(max(sqrt(diag(V))), 3), "\n")

ex <- ph_exposure(raw, cohs, AGES)
unflat <- function(p, pl) {
  for (nm in unique(names(p))) pl[[nm]] <- unname(p[names(p) == nm])
  pl
}
g <- function(p) ph_cells(unflat(p, pl), ex, cohs, AGES, with_obs = FALSE)$fitted_prob
base <- ph_cells(pl, ex, cohs, AGES, with_obs = TRUE)
p0 <- base$fitted_prob
cat("cells", length(p0), "\n")

## central differences, one column per fixed effect
th <- o$par
J <- matrix(NA_real_, length(p0), length(th))
t0 <- proc.time()[["elapsed"]]
for (k in seq_along(th)) {
  h <- 1e-4 * max(1, abs(th[k]))
  tp <- th; tp[k] <- th[k] + h
  tm <- th; tm[k] <- th[k] - h
  J[, k] <- (g(tp) - g(tm)) / (2 * h)
  if (k %% 5 == 0) cat("column", k, "of", length(th), " ",
                       round(proc.time()[["elapsed"]] - t0), "s\n")
}
## on the log-odds scale, so the interval cannot leave (0, 1)
Jl <- J / (p0 * (1 - p0))
se <- sqrt(pmax(rowSums((Jl %*% V) * Jl), 0))
lo <- plogis(qlogis(p0) - 1.96 * se)
hi <- plogis(qlogis(p0) + 1.96 * se)
cat("median relative half width", round(median((hi - lo) / (2 * p0)), 3), "\n")

out <- base[, c("x", "transition", "axis")]
out$band_lo <- lo; out$band_hi <- hi
out$cc <- CC; out$spec <- STEM; out$method <- "delta"
write.csv(out, here("fitted", "prob", sprintf("fit5_band_%s_%s.csv", CC, STEM)),
          row.names = FALSE)
cat("done", CC, nrow(out), "cells\n")
