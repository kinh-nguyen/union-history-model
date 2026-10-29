# Pre-marriage block with a phase after debut (code/msm_ph.cpp), one country.
#
#   Rscript code/fit_ph.R <CC>
#
# The likelihood is per woman, so the data are one row per woman-group, taken
# from the first episode of each group in data/cc8.csv.bz2. Shared parameters are
# started from the fixed-sharpness age-clock fit; the sharpness stays at one half
# and the two new quantities are the log multiplier of the marriage intensity in
# the first phase and the log rate of leaving it.
#
# Writes fitted/prob/fit5_prob_<CC>_ph.csv in the shape code/diag5.R writes, so
# code/fig_prefit_rp.R can draw it, and fitted/summary/fit5_sum_<CC>_ph.csv.
suppressPackageStartupMessages({ library(TMB); library(here); library(dplyr) })
for (.d in c("haz", "par", "prob", "summary", "other"))
  dir.create(here("fitted", .d), recursive = TRUE, showWarnings = FALSE)

CC <- commandArgs(trailingOnly = TRUE)[1]; if (is.na(CC)) CC <- "MW"
N_THREADS <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
state_order <- c("V", "X", "M", "D1", "W1", "R", "D2", "W2")
DLL <- "msm_ph"
# STEM separates the outputs of two runs of the same model, so a login-node run
# and a queued run never write the same file
STEM <- Sys.getenv("STEM", "ph")
cat("country", CC, "threads", N_THREADS, "\n")

## ---- data, one row per woman-group ------------------------------------------
cc <- vroom::vroom(here("data/cc8.csv.bz2"), show_col_types = FALSE)
raw <- cc %>% filter(cc == CC, sex == "f") %>%
  mutate(afs = if_else(is.na(afs), 0, afs), afm = if_else(is.na(afm), 0, afm),
         coh = (yob - 1975) / 10)
wom <- raw %>% filter(start == 0) %>%
  transmute(afs = as.integer(afs), afm = as.integer(afm), aai = as.integer(aai),
            coh, n)
cat("women-groups", nrow(wom), "weight", round(sum(wom$n)), "\n")
cat("path types",
    sum(wom$afs == 0 & wom$afm == 0), sum(wom$afs > 0 & wom$afm == 0),
    sum(wom$afm > 0 & wom$afs >= wom$afm),
    sum(wom$afm > 0 & wom$afs > 0 & wom$afs < wom$afm), "\n")

cohs0 <- sort(unique(wom$coh))
A_MAX <- 60L
dat <- list(afs = wom$afs, afm = wom$afm, aai = pmin(wom$aai, A_MAX),
            ic = match(wom$coh, cohs0) - 1L, coh_u = cohs0, a_max = A_MAX,
            n = wom$n, nopen = 0L)

## ---- starting values from the fixed-sharpness age-clock fit ------------------
# tauf_c0 is the free-sharpness fit on the age axis, which is the right start for
# all three transitions here, since the phase model puts every bowl on age
f0 <- here("fit", sprintf("fit5_%s_tauf_c0.rds", CC))
stopifnot(file.exists(f0))
pl0 <- readRDS(f0)$pl
# the peak curvature that reproduces the warm start's sharpness exactly
c2_0 <- 2 * exp(pl0$log_gam) * exp(pl0$log_lh_tau)^2 /
        (cosh(pl0$lh_del)^2 * exp(pl0$mu)^2 * exp(pl0$log_sigma)^2)
init <- list(mu = pl0$mu, log_sigma = pl0$log_sigma, lh_del = pl0$lh_del,
             log_c2 = log(c2_0), log_lh_tau = pl0$log_lh_tau,
             b_coh = pl0$b_coh, itc = pl0$itc[1:3], d1 = 0, log_lambda = 0)
map <- list()
# NOPHASE freezes the first phase out of the model, which makes X1 and X2 the same
# state and recovers the age-clock log hazard inside this likelihood, so the two
# can be compared on one scale
if (nzchar(Sys.getenv("NOPHASE"))) { init$d1 <- 0; map$d1 <- factor(NA) }

dyn.load(TMB::dynlib(here("code", DLL)))
# one tape, since the annual matrices are built once per age and cohort and the
# work per woman is a few matrix products; a parallel tape per thread was what
# exhausted memory on the login node
invisible(TMB::config(tape.parallel = 0, DLL = DLL))
TMB::openmp(1)
obj <- TMB::MakeADFun(dat, init, map = if (length(map)) map else NULL,
                      DLL = DLL, silent = TRUE)
cat("free parameters", length(obj$par), " start nll", obj$fn(obj$par), "\n")
ck0 <- here("fit", sprintf("fit5_%s_%s_pass.rds", CC, STEM))
if (nzchar(Sys.getenv("RESUME")) && file.exists(ck0)) {
  p0 <- readRDS(ck0)$par
  if (length(p0) == length(obj$par)) {
    obj$par <- p0
    cat("resumed from the saved pass, nll", obj$fn(p0), "\n")
  }
}

# the login node kills a process at three CPU-hours, so each pass is bounded and
# the fit is written after every pass; trace prints movement to the log
ctl <- list(trace = 10, eval.max = as.integer(Sys.getenv('EVALMAX', '120')),
            iter.max = as.integer(Sys.getenv('ITERMAX', '60')), abs.tol = 1e-12)
bnd <- function(par) {
  lo <- rep(-12, length(par)); hi <- rep(12, length(par))
  lo[names(par) == "mu"] <- -3;         hi[names(par) == "mu"] <- 5
  lo[names(par) == "log_lh_tau"] <- -3; hi[names(par) == "log_lh_tau"] <- 1.5
  lo[names(par) == "log_lambda"] <- -4; hi[names(par) == "log_lambda"] <- 3
  lo[names(par) == "d1"] <- -3;         hi[names(par) == "d1"] <- 6
  lo[names(par) == "log_c2"] <- -12;    hi[names(par) == "log_c2"] <- 6
  list(lo = lo, hi = hi)
}
b <- bnd(obj$par)
t0 <- proc.time()[["elapsed"]]
fit <- nlminb(pmin(pmax(obj$par, b$lo + 1e-6), b$hi - 1e-6), obj$fn, obj$gr,
              lower = b$lo, upper = b$hi, control = ctl)
ck <- here("fit", sprintf("fit5_%s_%s_pass.rds", CC, STEM))
saveRDS(fit, ck)
cat(sprintf("pass 1 nll %.4f max|grad| %.3g %.0fs\n", fit$objective,
            max(abs(obj$gr(fit$par))), proc.time()[["elapsed"]] - t0))
for (k in 2:12) {
  f2 <- nlminb(fit$par, obj$fn, obj$gr, lower = b$lo, upper = b$hi, control = ctl)
  gain <- fit$objective - f2$objective; fit <- f2
  saveRDS(fit, ck)
  cat(sprintf("pass %d nll %.4f max|grad| %.3g %.0fs\n", k, fit$objective,
              max(abs(obj$gr(fit$par))), proc.time()[["elapsed"]] - t0))
  if (gain < 1e-6 * max(1, abs(fit$objective)) ||
      max(abs(obj$gr(fit$par))) < 1e-3) break
}
t_fit <- proc.time()[["elapsed"]] - t0
cat(sprintf("nll %.4f max|grad| %.3g %.0fs\n", fit$objective,
            max(abs(obj$gr(fit$par))), t_fit))
pl <- init
for (nm in unique(names(fit$par))) pl[[nm]] <- unname(fit$par[names(fit$par) == nm])
cat("lambda", round(exp(pl$log_lambda), 4), " X1 multiplier",
    round(exp(pl$d1), 3), "\n")
saveRDS(list(par = fit$par, pl = pl, map = map, nll = fit$objective),
        here("fit", sprintf("fit5_%s_%s.rds", CC, STEM)))

## ---- fitted one-year cell probabilities, the definition of code/diag5.R -------
sasl <- function(a, mu_pk, sigma, del, gam, tau) {
  la <- log(a / mu_pk); z <- (exp(la) - 1) / sigma
  u <- z + sinh(del); r <- sqrt(u * u + 1)
  w <- ifelse(u > 0, u + r, 1 / (r - u)); wt <- exp(tau * log(w))
  s <- 0.5 * (wt - 1/wt) * cosh(tau * del) - 0.5 * (wt + 1/wt) * sinh(tau * del)
  -gam * s * s
}
gam_pl <- exp(pl$log_c2) * cosh(pl$lh_del)^2 * exp(2 * pl$mu) *
          exp(pl$log_sigma)^2 / (2 * exp(pl$log_lh_tau)^2)
rates_at <- function(a, coh_i) sapply(1:3, function(z)
  exp(pl$itc[z] + sasl(a, exp(pl$mu[z] + pl$b_coh[z] * coh_i),
                       exp(pl$log_sigma[z]), pl$lh_del[z],
                       gam_pl[z], exp(pl$log_lh_tau[z]))))
lam <- exp(pl$log_lambda); mult <- exp(pl$d1)
AGES <- 1:60
cohs <- sort(unique(wom$coh))
# annual transition matrix per age and cohort
Pm <- array(0, c(length(AGES), length(cohs), 4, 4))
for (ia in seq_along(AGES)) for (ic in seq_along(cohs)) {
  q <- rates_at(AGES[ia], cohs[ic])
  Q <- matrix(0, 4, 4)
  Q[1, 2] <- q[1]; Q[1, 4] <- q[2]
  Q[2, 3] <- lam;  Q[2, 4] <- q[3] * mult
  Q[3, 4] <- q[3]
  diag(Q) <- -rowSums(Q)
  Pm[ia, ic, , ] <- as.matrix(Matrix::expm(Q))
}
# phase mix among women still unmarried, by debut age and cohort
phase <- array(NA_real_, c(length(AGES), length(cohs), length(AGES), 2))
for (ic in seq_along(cohs)) for (s in AGES) {
  if (s >= max(AGES)) next
  v <- c(Pm[s, ic, 1, 2], Pm[s, ic, 1, 3])         # entering X during the debut year
  for (a in (s + 1):max(AGES)) {
    phase[s, ic, a, ] <- v / sum(v)
    P <- Pm[a, ic, , ]
    v <- c(v[1] * P[2, 2], v[1] * P[2, 3] + v[2] * P[3, 3])
  }
}

## person-years, as code/diag5.R builds them, from the episode file
e0 <- raw %>% filter(end - start >= 1, A %in% c("V", "X")) %>% mutate(dur = end - start)
ex <- e0[rep(seq_len(nrow(e0)), e0$dur), ]
ex$age <- unlist(lapply(seq_len(nrow(e0)), function(i) e0$start[i] + seq_len(e0$dur[i]) - 1))
ex$last <- ex$age == (ex$start + ex$dur - 1)
ex <- ex %>% filter(age >= 1, age <= max(AGES))
ex$ic <- match(ex$coh, cohs)
ex$ia <- match(ex$age, AGES)

inV <- ex$A == "V"
p_vx <- numeric(nrow(ex)); p_vm <- numeric(nrow(ex)); p_xm <- numeric(nrow(ex))
p_vx[inV] <- Pm[cbind(ex$ia[inV], ex$ic[inV], 1, 2)] + Pm[cbind(ex$ia[inV], ex$ic[inV], 1, 3)]
p_vm[inV] <- Pm[cbind(ex$ia[inV], ex$ic[inV], 1, 4)]
inX <- !inV
sX <- pmax(as.integer(ex$afs[inX]), 1L)
ph <- cbind(phase[cbind(sX, ex$ic[inX], ex$ia[inX], 1)],
            phase[cbind(sX, ex$ic[inX], ex$ia[inX], 2)])
ph[is.na(ph)] <- 0.5
p_xm[inX] <- ph[, 1] * Pm[cbind(ex$ia[inX], ex$ic[inX], 2, 4)] +
             ph[, 2] * Pm[cbind(ex$ia[inX], ex$ic[inX], 3, 4)]

e2 <- ex %>% mutate(p_vx = p_vx, p_vm = p_vm, p_xm = p_xm)
mk <- function(d, ev, pr, tr, ax, xv) d %>% mutate(x = xv, .ev = ev, .pr = pr) %>%
  group_by(x) %>%
  summarise(N = sum(n), events = sum(n * .ev), fitted_rate = NA_real_,
            fitted_prob = sum(n * .pr) / sum(n), .groups = "drop") %>%
  mutate(transition = tr, axis = ax)
tab <- bind_rows(
  mk(filter(e2, A == "V"), with(filter(e2, A == "V"), last & Z == "X"),
     filter(e2, A == "V")$p_vx, "V->X", "age", filter(e2, A == "V")$age),
  mk(filter(e2, A == "V"), with(filter(e2, A == "V"), last & Z == "M"),
     filter(e2, A == "V")$p_vm, "V->M", "age", filter(e2, A == "V")$age),
  mk(filter(e2, A == "X"), with(filter(e2, A == "X"), last & Z == "M"),
     filter(e2, A == "X")$p_xm, "X->M", "age", filter(e2, A == "X")$age),
  mk(filter(e2, A == "X", afs < 50), with(filter(e2, A == "X", afs < 50), last & Z == "M"),
     filter(e2, A == "X", afs < 50)$p_xm, "X->M", "dur",
     filter(e2, A == "X", afs < 50)$age - filter(e2, A == "X", afs < 50)$afs)) %>%
  mutate(obs = events / N,
         lo = qbeta(0.025, events + 0.5, N - events + 0.5),
         hi = qbeta(0.975, events + 0.5, N - events + 0.5),
         cc = CC, spec = STEM)
write.csv(tab, here("fitted", "prob", sprintf("fit5_prob_%s_%s.csv", CC, STEM)),
          row.names = FALSE)
write.csv(data.frame(cc = CC, spec = STEM, npar = length(fit$par),
                     nll = fit$objective, lambda = exp(pl$log_lambda),
                     x1_multiplier = exp(pl$d1),
                     maxgrad = max(abs(obj$gr(fit$par))), fit_s = round(t_fit, 1)),
          here("fitted", "summary", sprintf("fit5_sum_%s_%s.csv", CC, STEM)),
          row.names = FALSE)
cat("done", CC, nrow(tab), "rows\n")
