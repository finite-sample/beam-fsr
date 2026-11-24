# Each check compares against code that shares nothing with beam_path, so
# agreement is evidence rather than a restatement.
suppressMessages({
  library(bestsubset)
  source("R/methods.R")
})

set.seed(11)
for (cfg in list(c(100, 12, 0), c(100, 12, 0.7))) {
  d <- sim.xy(cfg[1], cfg[2], 1, rho = cfg[3], snr = 1, beta.type = 2)
  x <- d$x
  colnames(x) <- paste0("X", seq_len(ncol(x)))
  p <- ncol(x)

  # Width 1 is forward stepwise: match leaps' forward path exactly.
  fwd <- regsubsets(x, d$y, method = "forward", nvmax = p)
  b1 <- beam_path(x, d$y, 1)$best
  for (k in seq_len(p)) {
    cf <- coef(fwd, k)
    want <- numeric(p + 1)
    want[c(1, match(names(cf)[-1], colnames(x)) + 1)] <- cf
    stopifnot(isTRUE(all.equal(b1[, k + 1], want, tolerance = 1e-8)))
  }

  # A beam wide enough to hold every subset is best subset: its best member's
  # RSS must equal leaps' exhaustive optimum at every size.
  wide <- beam_path(x, d$y, choose(p, p %/% 2))
  ex <- best_subset_path(x, d$y)
  stopifnot(isTRUE(all.equal(wide$rss, ex$rss, tolerance = 1e-10)))

  # The averages are the stated averages of the members' OLS fits, refit
  # here with lm() on the member subsets.
  bp <- beam_path(x, d$y, 5)
  for (k in c(2, 5, 9)) {
    fits <- sapply(bp$beams[[k]], function(S) {
      b <- numeric(p + 1)
      b[c(1, S + 1)] <- coef(lm(d$y ~ x[, S, drop = FALSE]))
      b
    })
    rss <- sapply(bp$beams[[k]], function(S) {
      sum(resid(lm(d$y ~ x[, S, drop = FALSE]))^2)
    })
    w <- rss^(-nrow(x) / 2)
    stopifnot(isTRUE(all.equal(bp$avg[, k + 1], rowMeans(fits), tolerance = 1e-8)))
    stopifnot(isTRUE(all.equal(bp$bic[, k + 1], drop(fits %*% (w / sum(w))),
                               tolerance = 1e-8)))
    stopifnot(isTRUE(all.equal(bp$rss[k + 1], min(rss), tolerance = 1e-10)))
  }
}

# Complete subset regression with every subset enumerated equals a brute
# force average over combn().
d <- sim.xy(60, 6, 1, snr = 1)
cs <- csr_path(d$x, d$y, kmax = 3)
brute <- rowMeans(apply(combn(6, 2), 2, function(S) {
  b <- numeric(7)
  b[c(1, S + 1)] <- coef(lm(d$y ~ d$x[, S]))
  b
}))
stopifnot(isTRUE(all.equal(cs[, 3], brute, tolerance = 1e-8)))

# Highly correlated columns: an AR(1) design with rho = 0.99, and one with
# near-duplicate column pairs. beam_path updates its projections one column
# at a time; the reference below scores every candidate with a fresh lm.fit,
# applies the same collinearity rule and drops duplicate subsets the same way.
beam_ref <- function(x, y, width, tol = 1e-8) {
  rss_of <- function(S, v) sum(lm.fit(cbind(1, x[, S, drop = FALSE]), v)$residuals^2)
  ss0 <- colSums(sweep(x, 2, colMeans(x))^2)
  beam <- list(integer(0))
  out <- list(beams = list(), rss = sum((y - mean(y))^2))
  for (k in seq_len(ncol(x))) {
    cand <- list()
    for (S in beam) for (j in setdiff(seq_len(ncol(x)), S)) {
      if (rss_of(S, x[, j]) > tol * ss0[j]) cand <- c(cand, list(sort(c(S, j))))
    }
    if (!length(cand)) break
    cand <- cand[!duplicated(vapply(cand, paste, "", collapse = ","))]
    score <- vapply(cand, rss_of, 0, v = y)
    o <- order(score)[seq_len(min(width, length(cand)))]
    beam <- cand[o]
    out$beams[[k]] <- beam
    out$rss[k + 1] <- score[o[1]]
  }
  out
}
set.seed(12)
d <- sim.xy(100, 10, 1, rho = 0.99, snr = 6, beta.type = 2)
near <- sim.xy(100, 10, 1, rho = 0.5, snr = 6, beta.type = 2)
near$x[, 2] <- near$x[, 1] + 1e-3 * rnorm(100)
near$x[, 6] <- near$x[, 5] + 1e-5 * rnorm(100)
for (dd in list(d, near)) for (w in c(1, 4)) {
  got <- beam_path(dd$x, dd$y, w)
  want <- beam_ref(dd$x, dd$y, w)
  stopifnot(identical(got$beams, want$beams))
  stopifnot(isTRUE(all.equal(got$rss[!is.na(got$rss)], want$rss, tolerance = 1e-10)))
}

# Complete subset regression fits most subsets by Cholesky and falls back to
# QR when the subset is ill conditioned. With a pair of columns 1e-2 apart
# (Cholesky path) and 1e-7 apart (QR path, where Cholesky alone was off by
# about 1e-3), every size must match lm.fit over the same enumerated subsets.
for (eps in c(1e-2, 1e-7)) {
  dd <- sim.xy(100, 8, 1, rho = 0.5, snr = 6, beta.type = 2)
  dd$x[, 2] <- dd$x[, 1] + eps * rnorm(100)
  cs <- csr_path(dd$x, dd$y)
  for (k in 1:6) {
    brute <- rowMeans(apply(combn(8, k), 2, function(S) ols_coef(dd$x, dd$y, S, 8)))
    stopifnot(isTRUE(all.equal(cs[, k + 1], brute, tolerance = 1e-9)))
  }
}

# Bootstrap designs with p close to n are rank deficient: no error, and the
# path stops before the residual degrees of freedom run out.
d <- sim.xy(60, 50, 1, snr = 1)
i <- sample.int(60, 60, replace = TRUE)
path <- beam_path(d$x[i, ], d$y[i], 3)$best
stopifnot(max(colSums(path[-1, ] != 0)) <= length(unique(i)) - 1)

cat("all method tests passed\n")
