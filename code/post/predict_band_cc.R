# Posterior bands for binned state shares, for fit_post_crisis.R fits.
# Draws from N(p_hat, H^+) restricted to directions with curvature above CUT of the largest (use 1e-6).
# Predictions use the transposed-generator report, which is the correct one in double mode.
suppressPackageStartupMessages({ library(TMB); library(jsonlite); library(splines) })
set.seed(1)
src <- readLines("fit_post_crisis.R"); cut <- grep("^obj <- MakeADFun", src)
ARGV <- commandArgs(TRUE); CCX <- ARGV[1]; nm0 <- ARGV[2]; ND <- if (length(ARGV) > 2) as.integer(ARGV[3]) else 400
CUT <- if (length(ARGV) > 3) as.numeric(ARGV[4]) else 1e-6
a <- c(CCX, nm0, "-", "4"); commandArgs <- function(...) a
eval(parse(text = src[1:(cut - 1)]))
tag <- gsub("[^A-Za-z0-9]+", "_", nm0)
x <- readRDS(sprintf("fit2d_%s_%s_par.rds", CCX, tag)); p <- x$par
obj0 <- MakeADFun(data, pars, map = map, DLL = DLLN, silent = TRUE)
stopifnot(length(obj0$par) == length(p))
cat(sprintf("REBUILD objective %.4f stored %.4f diff %.2e\n", obj0$fn(p), x$objective, obj0$fn(p) - x$objective)); flush.console()
H <- optimHess(p, obj0$fn, obj0$gr); H <- (H + t(H)) / 2
e <- eigen(H, symmetric = TRUE); kp <- e$values > CUT * max(e$values)
L <- e$vectors[, kp] %*% diag(1 / sqrt(e$values[kp]))
cat("sampling directions", sum(kp), "of", length(p), "\n"); flush.console()
d1 <- data; d1$tq <- 1L
obj1 <- MakeADFun(d1, pars, map = map, DLL = DLLN, silent = TRUE)
PRhat <- obj1$report(p)$PR
cat("rowsum", range(rowSums(PRhat)), "\n"); flush.console()
yom <- gr$gy0[rw$g + 1] + 1930; afm <- gr$start[rw$g + 1]
keys <- list(dur = pmin((rw$dur - 1) %/% 2, 16) * 2 + 1, yom = pmin(pmax((yom %/% 5) * 5, 1960), 2015), afm = pmin(pmax(afm, 14), 36))
# the three outcome groups of the figure are summed inside each draw, so their
# quantiles come from the grouped quantity and not from a sum of quantiles
grp <- function(P) cbind(P, divorced = P[, 2] + P[, 5], widowed = P[, 3] + P[, 6],
                        remarried = P[, 4])
bin <- function(P0, k) { P <- grp(P0); o <- lapply(1:9, function(s) tapply(rw$n * P[, s], k, sum)); N <- tapply(rw$n, k, sum)
  m <- do.call(cbind, o) / as.numeric(N); rownames(m) <- names(N); m }
S <- c("M", "D1", "W1", "R", "D2", "W2", "divorced", "widowed", "remarried")
A <- lapply(keys, function(k) array(NA_real_, c(length(unique(k)), 9, ND)))
ok <- 0; tried <- 0
while (ok < ND && tried < 5 * ND) {
  tried <- tried + 1
  P <- obj1$report(as.numeric(p + L %*% rnorm(sum(kp))))$PR
  if (!all(is.finite(P)) || max(abs(rowSums(P) - 1)) > 1e-3) next
  ok <- ok + 1
  for (bn in names(keys)) A[[bn]][, , ok] <- bin(P, keys[[bn]])
}
cat("accepted", ok, "of", tried, "draws\n")
out <- list()
for (bn in names(keys)) {
  k <- keys[[bn]]; N <- tapply(rw$n, k, sum); nr_ <- tapply(rep(1, nrow(rw)), k, sum)
  obsm <- bin(model.matrix(~ factor(rw$Z, levels = 0:5) - 1), k); hat <- bin(PRhat, k)
  Ab <- A[[bn]][, , seq_len(ok), drop = FALSE]
  lo <- apply(Ab, 1:2, quantile, 0.025); hi <- apply(Ab, 1:2, quantile, 0.975)
  for (s in 1:9) out[[length(out) + 1]] <- data.frame(by = bn, key = as.numeric(rownames(hat)), state = S[s],
      N = as.numeric(N), rows = as.numeric(nr_), obs = obsm[, s], fit = hat[, s], lo = lo[, s], hi = hi[, s], ndraw = ok)
}
write.csv(do.call(rbind, out), sprintf("band_%s_%s.csv", CCX, tag), row.names = FALSE)
cat("written\n")
