# Model-versus-data diagnostics for the fit5 specifications, on the right scale.
#
#   Rscript diag5.R <CC>
#
# fit4.R and fit5.R report the fitted hazard as an exposure-weighted mean RATE
# and compare it against the observed cell frequency events/N, which is an annual
# PROBABILITY. The two agree only while the rate is small: at 0.2/yr the gap is
# 10 per cent, at 1/yr it is 46 per cent and at 3/yr it is a factor of three. The
# incumbent's rates are small enough that the convention did not matter, but a
# hazard indexed on duration since debut is of order 1/yr immediately after
# debut, so on the rate convention the new specifications appear to overshoot
# the data badly while their likelihood improves by hundreds of units. This
# script recomputes the fitted side as the exact one-year cell probability the
# likelihood uses, so the comparison is like for like. Both columns are written.
#
# The three-state chain V -> X -> M with rates a (V->X), b (V->M), c (X->M) and
# an outflow r from M has a closed-form annual transition matrix, which is used
# here instead of a numerical matrix exponential per episode-year (there are up
# to a million of them) and checked against Matrix::expm on a random sample.
#   V(1)    = exp(-s),  s = a + b
#   P[V,X]  = a D(s, c)
#   P[X,M]  = c D(c, r)
#   P[V,M]  = b D(s, r) + a c (D(c, r) - D(s, r)) / (s - c)
# with D(x, y) = (exp(-y) - exp(-x))/(x - y) = exp(-x) expm1(x - y)/(x - y),
# written that way because the ratio is 0/0 at x = y and the limit is exp(-x).
suppressPackageStartupMessages({ library(TMB); library(here); library(dplyr); library(Matrix) })
for (.d in c("haz", "par", "prob", "summary", "other")) dir.create(here("fitted", .d), recursive = TRUE, showWarnings = FALSE)


CC <- commandArgs(trailingOnly = TRUE)[1]; if (is.na(CC)) CC <- "ST"
SUFFIX <- Sys.getenv("SUFFIX", "")   # read the fits of a suffixed re-run
SPECS  <- Sys.getenv("SPECS", "")
state_order <- c("V", "X", "M", "D1", "W1", "R", "D2", "W2")

DATAFILE <- Sys.getenv("DATAFILE", "data/cc8.csv.bz2")
cc <- vroom::vroom(here(DATAFILE), show_col_types = FALSE)
raw <- cc %>% filter(cc == CC, sex == "f") %>%
  mutate(start = if_else(start == 0, 1, start)) %>%
  mutate(A = match(A, state_order) - 1, Z = match(Z, state_order) - 1,
         dur = end - start, afs = if_else(afs == 0 | is.na(afs), 100, afs),
         coh = (yob - 1975) / 10) %>%
  filter(dur >= 1, A <= 1) %>% select(A, Z, start, dur, afs, coh, n)
ex <- raw[rep(seq_len(nrow(raw)), raw$dur), ]
ex$age <- unlist(lapply(seq_len(nrow(raw)), function(i)
  raw$start[i] + seq_len(raw$dur[i]) - 1))
ex$last <- ex$age == (ex$start + ex$dur - 1)
cat("country", CC, "episode-years", nrow(ex), "\n")

Dfun <- function(x, y) {            # (exp(-y) - exp(-x))/(x - y), stable at x = y
  d <- x - y
  ifelse(abs(d) < 1e-9, exp(-x), exp(-x) * expm1(d) / d)
}
sp <- function(x) log1p(exp(x))

mk_haz <- function(pl, clock, lh = 1L) {
  omega <- 1 / (1 + exp(-pl$logit_omega))
  # additive background intensity on X->M (fit_floor.R); absent in older fits
  h0 <- exp(if (is.null(pl$log_h0)) -30 else pl$log_h0)
  if (lh == 0L) return(function(a, z, afs_i, coh_i) {
    clk <- if (z == 3 && clock %in% c(1, 2, 4, 6)) sp(a - afs_i)
           else if (z == 3 && clock %in% c(3, 5)) sp(a - omega * afs_i) else a
    mu  <- pl$mu[z] + pl$b_coh[z] * coh_i +
           (if (z == 3 && clock == 6) pl$b_mafs * (pmin(afs_i, 30) - 16) else 0)
    sg  <- exp(pl$log_sigma[z]); nu <- exp(pl$log_nu[z]); ta <- exp(pl$log_tau[z])
    x   <- log(clk); zz <- (x - mu) / (sg * ta)
    t   <- ta * log(zz + sqrt(zz * zz + 1)) - nu
    lr  <- -log(sg) - 0.5 * log(2 * pi) - 0.5 * log(1 + zz * zz) +
            log(cosh(t)) - 0.5 * sinh(t)^2 - x
    e   <- pl$itc[z] + lr - log(1 - pnorm(sinh(t)) + .Machine$double.eps)
    if (z == 3) {
      if (clock %in% c(0, 3)) e <- e + pl$b_tx * sp(a - afs_i) +
                                   pl$b_ld * log1p(sp(a - afs_i)) -
                                   log1p(exp(pl$log_xi) * sp(a - afs_i))
      else if (clock == 1)    e <- e + pl$b_ax * (a - 20)
      else if (clock %in% c(4, 5, 6)) e <- e + pl$b_afs * (pmin(afs_i, 30) - 16)
    }
    exp(e) + if (z == 3) h0 else 0
  })
  function(a, z, afs_i, coh_i) {
    pk  <- exp(pl$mu[z] + pl$b_coh[z] * coh_i +
                 (if (z == 3 && clock == 6) pl$b_mafs * (pmin(afs_i, 30) - 16) else 0))
    clk <- if (z == 3 && clock %in% c(1, 2, 4, 6)) sp(a - afs_i)
           else if (z == 3 && clock %in% c(3, 5)) sp(a - omega * afs_i) else a
    pw  <- exp(pl$log_lh_p[z]); tau <- exp(pl$log_lh_tau[z])
    zz  <- ((clk / pk)^pw - 1) / (pw * exp(pl$log_sigma[z]))
    u   <- zz + sinh(pl$lh_del[z]); r <- sqrt(u * u + 1)
    w   <- ifelse(u > 0, u + r, 1 / (r - u)); wt <- exp(tau * log(w))
    s   <- 0.5 * (wt - 1/wt) * cosh(tau * pl$lh_del[z]) -
           0.5 * (wt + 1/wt) * sinh(tau * pl$lh_del[z])
    e   <- pl$itc[z] - exp(pl$log_gam[z]) * s * s
    if (z == 3) {
      if (clock %in% c(0, 3)) e <- e + pl$b_tx * sp(a - afs_i) +
                                   pl$b_ld * log1p(sp(a - afs_i)) -
                                   log1p(exp(pl$log_xi) * sp(a - afs_i))
      else if (clock == 1)    e <- e + pl$b_ax * (a - 20)
      else if (clock %in% c(4, 5, 6)) e <- e + pl$b_afs * (pmin(afs_i, 30) - 16)
    }
    exp(e) + if (z == 3) h0 else 0
  }
}

specs <- if (nzchar(SPECS)) strsplit(SPECS, ",")[[1]] else
  c("shash_leak", "shash", "lh1_leak", "lh1", "base_leak", "base", "tauf_c0",
    "tauf_c3", "tauf_c4", "tauf_c5", "taupf_c5", "tauf_c6", "ld_base", "ld_tauf", "fr_base")
out <- list()
for (sp_nm in specs) {
  f <- here("fit", sprintf("fit5_%s_%s%s.rds", CC, sp_nm, SUFFIX))
  if (!file.exists(f)) { cat("missing", sp_nm, "\n"); next }
  o <- readRDS(f); pl <- o$pl; clock <- o$spec$clock
  r_leak <- if (isTRUE(o$spec$leak)) 2 else 0     # M outflow: exp(0) + exp(0)
  hz <- mk_haz(pl, clock, if (is.null(o$spec$lh)) 1L else o$spec$lh)
  a_r <- hz(ex$age, 1, ex$afs, ex$coh)
  b_r <- hz(ex$age, 2, ex$afs, ex$coh)
  c_r <- hz(ex$age, 3, ex$afs, ex$coh)
  s_r <- a_r + b_r
  p_vx <- a_r * Dfun(s_r, c_r)
  p_xm <- c_r * Dfun(c_r, r_leak)
  p_vm_dir <- b_r * Dfun(s_r, r_leak)     # married without entering X
  p_vm_rte <- a_r * c_r * (Dfun(c_r, r_leak) - Dfun(s_r, r_leak)) / (s_r - c_r)
  p_vm <- p_vm_dir + p_vm_rte             # the route share is pi
  # the closed form against a numerical matrix exponential, on a random sample
  set.seed(1); k <- sample(nrow(ex), min(300, nrow(ex)))
  err <- max(sapply(k, function(i) {
    Q <- matrix(0, 3, 3); Q[1, 2] <- a_r[i]; Q[1, 3] <- b_r[i]; Q[2, 3] <- c_r[i]
    diag(Q) <- c(-s_r[i], -c_r[i], -r_leak)
    P <- as.matrix(Matrix::expm(Q))
    max(abs(c(P[1, 2] - p_vx[i], P[1, 3] - p_vm[i], P[2, 3] - p_xm[i])))
  }))
  stopifnot(err < 1e-9)
  e2 <- ex %>% mutate(rate = if_else(A == 0, a_r, c_r),
                      rate2 = b_r, prob = if_else(A == 0, p_vx, p_xm),
                      prob2 = p_vm, prob2_rte = p_vm_rte)
  tab <- bind_rows(
    e2 %>% filter(A == 0) %>% group_by(x = age) %>%
      summarise(N = sum(n), events = sum(n * (last & Z == 1)),
                fitted_rate = sum(n * rate) / sum(n),
                fitted_prob = sum(n * prob) / sum(n), .groups = "drop") %>%
      mutate(transition = "V->X", axis = "age"),
    e2 %>% filter(A == 0) %>% group_by(x = age) %>%
      summarise(N = sum(n), events = sum(n * (last & Z == 2)),
                fitted_rate = sum(n * rate2) / sum(n),
                fitted_prob = sum(n * prob2) / sum(n),
                fitted_prob_route = sum(n * prob2_rte) / sum(n), .groups = "drop") %>%
      mutate(transition = "V->M", axis = "age"),
    e2 %>% filter(A == 1) %>% group_by(x = age) %>%
      summarise(N = sum(n), events = sum(n * (last & Z == 2)),
                fitted_rate = sum(n * rate) / sum(n),
                fitted_prob = sum(n * prob) / sum(n), .groups = "drop") %>%
      mutate(transition = "X->M", axis = "age"),
    e2 %>% filter(A == 1, afs < 50) %>% group_by(x = age - afs) %>%
      summarise(N = sum(n), events = sum(n * (last & Z == 2)),
                fitted_rate = sum(n * rate) / sum(n),
                fitted_prob = sum(n * prob) / sum(n), .groups = "drop") %>%
      mutate(transition = "X->M", axis = "dur")) %>%
    mutate(obs = events / N,
           lo = qbeta(0.025, events + 0.5, N - events + 0.5),
           hi = qbeta(0.975, events + 0.5, N - events + 0.5),
           cc = CC, spec = paste0(sp_nm, SUFFIX))
  cat(sprintf("%-10s clock %d leak %d  expm check %.2e  rows %d\n",
              sp_nm, clock, r_leak, err, nrow(tab)))
  out[[sp_nm]] <- tab
}
write.csv(bind_rows(out), here("fitted", "prob", sprintf("fit5_prob_%s%s.csv", CC, SUFFIX)), row.names = FALSE)
cat("done", CC, "\n")
