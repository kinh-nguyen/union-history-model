# Annual probability of divorce, widowhood and remarriage by single year of age,
# post-marriage block, exposure weighted like the pre-marriage figures.
#
#   Rscript post_age.R <CC> <ops> [NDRAW]
#
# The surveys record no transition times after marriage, so this margin has no
# observed counterpart; it is the model's own schedule with its 95 per cent
# interval. The template msm_post_rep.cpp reports the three intensities a row is
# exposed to in its own year, and the probabilities are
#   divorce     q_D / (q_D + q_W) * (1 - exp(-(q_D + q_W)))
#   widowhood   q_W / (q_D + q_W) * (1 - exp(-(q_D + q_W)))
#   remarriage  1 - exp(-q_R)
# Each is averaged over women at risk of it, that is women in a first union for
# the first two and women whose first union has ended for the third, with the
# model's own occupancies as the weights.
suppressPackageStartupMessages({ library(TMB); library(jsonlite); library(splines) })
set.seed(2)
src <- readLines("fit_post_crisis.R"); cut <- grep("^obj <- MakeADFun", src)
ARGV <- commandArgs(TRUE); CCX <- ARGV[1]; nm0 <- ARGV[2]
ND <- if (length(ARGV) > 2) as.integer(ARGV[3]) else 200
CUT <- 1e-6
a <- c(CCX, nm0, "-", "4"); commandArgs <- function(...) a
eval(parse(text = src[1:(cut - 1)]))
tag <- gsub("[^A-Za-z0-9]+", "_", nm0)
x <- readRDS(sprintf("fit2d_%s_%s_par.rds", CCX, tag)); p <- x$par

dyn.load(dynlib("msm_post_rep"))
obj0 <- MakeADFun(data, pars, map = map, DLL = "msm_post_rep", silent = TRUE)
stopifnot(length(obj0$par) == length(p))
cat(sprintf("REBUILD %.4f stored %.4f diff %.2e\n", obj0$fn(p), x$objective,
            obj0$fn(p) - x$objective)); flush.console()
H <- optimHess(p, obj0$fn, obj0$gr); H <- (H + t(H)) / 2
e <- eigen(H, symmetric = TRUE); kp <- e$values > CUT * max(e$values)
L <- e$vectors[, kp] %*% diag(1 / sqrt(e$values[kp]))
cat("sampling directions", sum(kp), "of", length(p), "\n"); flush.console()
# the report is read from the transposed generator, which is the correct
# orientation in double mode; without it the occupancy stays in the first state
d1 <- data; d1$tq <- 1L
obj1 <- MakeADFun(d1, pars, map = map, DLL = "msm_post_rep", silent = TRUE)
cat("rowsum", range(rowSums(obj1$report(p)$PR)), "\n"); flush.console()

age <- gr$start[rw$g + 1] + rw$dur - 1
key <- pmin(pmax(age, 15), 60)
curves <- function(par) {
  rp <- obj1$report(par)
  Q <- rp$QH; PR <- rp$PR
  s <- Q[, 1] + Q[, 2]; pe <- 1 - exp(-s)
  pd <- ifelse(s > 0, Q[, 1] / s, 0) * pe
  pw <- ifelse(s > 0, Q[, 2] / s, 0) * pe
  pr <- 1 - exp(-Q[, 3])
  wM <- rw$n * PR[, 1]                      # in a first union
  wE <- rw$n * (PR[, 2] + PR[, 3])          # first union ended, not remarried
  cbind(divorce = tapply(wM * pd, key, sum) / tapply(wM, key, sum),
        widowhood = tapply(wM * pw, key, sum) / tapply(wM, key, sum),
        remarriage = tapply(wE * pr, key, sum) / tapply(wE, key, sum))
}
hat <- curves(p)
A <- array(NA_real_, c(nrow(hat), 3, ND))
ok <- 0; tried <- 0
while (ok < ND && tried < 5 * ND) {
  tried <- tried + 1
  cv <- try(curves(as.numeric(p + L %*% rnorm(sum(kp)))), silent = TRUE)
  if (inherits(cv, "try-error") || !all(is.finite(cv))) next
  ok <- ok + 1; A[, , ok] <- cv
}
cat("accepted", ok, "of", tried, "draws\n")
Ab <- A[, , seq_len(ok), drop = FALSE]
lo <- apply(Ab, 1:2, quantile, 0.025); hi <- apply(Ab, 1:2, quantile, 0.975)
wM <- tapply(rw$n, key, sum)
out <- do.call(rbind, lapply(1:3, function(s) data.frame(
  cc = CCX, rate = colnames(hat)[s], age = as.numeric(rownames(hat)),
  fit = hat[, s], lo = lo[, s], hi = hi[, s], N = as.numeric(wM), ndraw = ok)))
write.csv(out, sprintf("postage_%s_%s.csv", CCX, tag), row.names = FALSE)
cat("written", nrow(out), "rows\n")
