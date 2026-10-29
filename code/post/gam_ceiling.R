# multinomial ceiling for the post-marriage block (RW): nested binary GAMs, same rows and weights as TMB
suppressPackageStartupMessages(library(mgcv))
d <- read.csv("rowP_RW.csv")
br <- list(b1 = list(sub = rep(TRUE, nrow(d)), y = d$Z == 0),     # in first union vs ended
           b2 = list(sub = d$Z != 0,          y = d$Z %in% 3:5),  # later union (R,D2,W2) vs D1/W1
           b3 = list(sub = d$Z %in% 1:2,      y = d$Z == 2),      # W1 vs D1
           b4 = list(sub = d$Z %in% 3:5,      y = d$Z == 3),      # R vs D2/W2
           b5 = list(sub = d$Z %in% 4:5,      y = d$Z == 5))      # W2 vs D2
fm <- list(G1 = y ~ te(dur, yom, k = c(8, 8)) + s(afm, k = 6),
           G2 = y ~ te(dur, yom, afm, k = c(6, 6, 4)))
out <- list(); PP <- list()
for (m in names(fm)) {
  p <- matrix(NA, nrow(d), 5)
  for (b in names(br)) {
    i <- br[[b]]$sub; dd <- d[i, ]; dd$y <- as.numeric(br[[b]]$y[i])
    t0 <- Sys.time()
    f <- suppressWarnings(bam(fm[[m]], family = binomial, weights = n, data = dd, method = "fREML", discrete = FALSE))
    p[, match(b, names(br))] <- predict(f, newdata = d, type = "response")
    pi <- fitted(f); ll <- sum(dd$n * (dd$y * log(pi) + (1 - dd$y) * log(1 - pi)))
    out[[length(out) + 1]] <- data.frame(model = m, branch = b, rows = nrow(dd), N = sum(dd$n), ll = ll,
                                         edf = sum(f$edf), secs = as.numeric(Sys.time() - t0, units = "secs"))
    cat(m, b, round(ll, 2), round(sum(f$edf), 2), "\n")
  }
  q <- with(list(p1 = p[,1], p2 = p[,2], p3 = p[,3], p4 = p[,4], p5 = p[,5]),
            cbind(p1, (1-p1)*(1-p2)*(1-p3), (1-p1)*(1-p2)*p3, (1-p1)*p2*p4, (1-p1)*p2*(1-p4)*(1-p5), (1-p1)*p2*(1-p4)*p5))
  colnames(q) <- paste0(m, 0:5); PP[[m]] <- q
}
res <- do.call(rbind, out); write.csv(res, "gam_ceiling_RW.csv", row.names = FALSE)
write.csv(cbind(d[, c("g", "dur", "Z", "n", "afm", "yom")], PP$G1, PP$G2), "gam_rowP_RW.csv", row.names = FALSE)
print(aggregate(cbind(ll, edf) ~ model, res, sum))
