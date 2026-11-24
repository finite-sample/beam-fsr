library(glmnet)
library(leaps)
library(BMA)

# Every method returns a coefficient matrix with the intercept in row 1 and
# one column per model size k = 0, 1, ..., kmax (or per tuning value), so a
# size or tuning value can be picked on validation data with pick().

ols_coef <- function(x, y, vars, p) {
  b <- numeric(p + 1)
  b[c(1, vars + 1)] <- lm.fit(cbind(1, x[, vars, drop = FALSE]), y)$coefficients
  b
}

# Beam search over subsets for least squares. Every member of the beam is
# extended by each unused variable, candidates are scored by training RSS,
# duplicate subsets are dropped, and the `width` best become the next beam.
# width = 1 is forward stepwise; a width at least choose(p, k) at every k is
# best subset. Hand-rolled because no maintained R package exposes the beam
# itself (leaps returns only its nbest table, OKRidge is Python and ridge).
#
# For each size the beam is turned into one coefficient vector three ways:
#   best  the lowest-RSS member (what beam search usually reports)
#   avg   equal-weight average of all members' OLS fits
#   bic   BIC-weighted average; members share a size, so the weights are
#         proportional to RSS^(-n/2)
beam_path <- function(x, y, width, kmax = min(ncol(x), nrow(x) - 2),
                      tol = 1e-8) {
  n <- nrow(x)
  p <- ncol(x)
  xm <- colMeans(x)
  # Each beam member carries the centered [x y] with its own columns
  # projected out, plus the rows of its R factor, both built by modified
  # Gram-Schmidt one column at a time. A child then costs one projection of
  # its parent's matrix instead of a fresh QR, and its OLS fit is a
  # triangular solve. MGS on the augmented [x y] is backward stable for least
  # squares (Bjorck 1967), and nothing here forms x'x, so correlated columns
  # lose no more accuracy than the QR they replace.
  z0 <- cbind(sweep(x, 2, xm), y - mean(y))
  ss0 <- colSums(z0[, seq_len(p), drop = FALSE]^2)
  out <- lapply(c(best = 1, avg = 1, bic = 1), function(i) {
    m <- matrix(0, p + 1, kmax + 1)
    m[1, 1] <- mean(y)
    m
  })
  beams <- vector("list", kmax)
  rss <- c(sum(z0[, p + 1]^2), rep(NA_real_, kmax))
  beam <- list(list(S = integer(0), z = z0, R = matrix(0, 0, p + 1)))
  last <- 0
  for (k in seq_len(kmax)) {
    from <- integer(0)
    add <- integer(0)
    score <- numeric(0)
    for (m in seq_along(beam)) {
      z <- beam[[m]]$z
      ss <- colSums(z[, seq_len(p), drop = FALSE]^2)
      # A column already explained by S (collinear, or a duplicated bootstrap
      # row pattern) adds nothing and would make the fit singular.
      ok <- setdiff(which(ss > tol * ss0), beam[[m]]$S)
      if (!length(ok)) next
      r <- z[, p + 1]
      gain <- drop(crossprod(z[, ok, drop = FALSE], r))^2 / ss[ok]
      from <- c(from, rep(m, length(ok)))
      add <- c(add, ok)
      score <- c(score, sum(r^2) - gain)
    }
    if (!length(score)) break
    # Walk the candidates best first, dropping subsets already reached from
    # another member, until the beam is full; only the survivors are built.
    keys <- character(0)
    pick_c <- integer(0)
    for (cc in order(score)) {
      key <- paste(sort(c(beam[[from[cc]]]$S, add[cc])), collapse = ",")
      if (key %in% keys) next
      keys <- c(keys, key)
      pick_c <- c(pick_c, cc)
      if (length(pick_c) == width) break
    }
    beam <- lapply(pick_c, function(cc) {
      par <- beam[[from[cc]]]
      j <- add[cc]
      q <- par$z[, j] / sqrt(sum(par$z[, j]^2))
      row <- drop(crossprod(q, par$z))
      list(S = c(par$S, j), z = par$z - tcrossprod(q, row), R = rbind(par$R, row))
    })
    beams[[k]] <- lapply(beam, function(b) sort(b$S))
    sc <- score[pick_c]
    rss[k + 1] <- sc[1]
    fits <- vapply(beam, function(b) {
      cf <- backsolve(b$R[, b$S, drop = FALSE], b$R[, p + 1])
      f <- numeric(p + 1)
      f[c(1, b$S + 1)] <- c(mean(y) - sum(xm[b$S] * cf), cf)
      f
    }, numeric(p + 1))
    # A saturated fit has RSS ~ 0, and rounding can push it below zero.
    lw <- -n / 2 * log(pmax(sc, .Machine$double.eps * rss[1]))
    w <- exp(lw - max(lw))
    out$best[, k + 1] <- fits[, 1]
    out$avg[, k + 1] <- rowMeans(fits)
    out$bic[, k + 1] <- drop(fits %*% (w / sum(w)))
    last <- k
  }
  # Sizes the search could not reach repeat the largest model, so paths from
  # different bootstrap samples can be averaged column by column.
  if (last < kmax) {
    for (m in names(out)) out[[m]][, (last + 2):(kmax + 1)] <- out[[m]][, last + 1]
  }
  c(out, list(rss = rss, beams = beams[seq_len(last)]))
}

# Breiman (1994, TR 421 sec. 5): forward stepwise refit on bootstrap samples,
# coefficient paths averaged size by size.
bag_fsr <- function(x, y, B, kmax = min(ncol(x), nrow(x) - 2)) {
  n <- nrow(x)
  acc <- 0
  for (b in seq_len(B)) {
    i <- sample.int(n, n, replace = TRUE)
    acc <- acc + beam_path(x[i, ], y[i], 1, kmax)$best
  }
  acc / B
}

# Exhaustive best subset (leaps branch and bound); the ceiling for any beam
# whose members are ranked by training RSS.
best_subset_path <- function(x, y, kmax = ncol(x)) {
  p <- ncol(x)
  s <- summary(regsubsets(x, y, nvmax = kmax, method = "exhaustive",
                          really.big = TRUE))
  path <- matrix(0, p + 1, kmax + 1)
  path[1, 1] <- mean(y)
  for (k in seq_len(kmax)) path[, k + 1] <- ols_coef(x, y, which(s$which[k, -1]), p)
  list(path = path, rss = c(sum((y - mean(y))^2), s$rss))
}

# Complete subset regressions (Elliott, Gargano & Timmermann 2013): the
# equal-weight average of OLS fits over all k-variable subsets. When
# choose(p, k) exceeds max_subsets, the average is over a uniform random
# sample of max_subsets subsets instead; that shortcut is ours, because the
# exact average is infeasible at moderate k once p is much above 20.
#
# Each subset is fit from the centered cross-products by a Cholesky solve,
# which squares the condition number of x[, S]. A subset whose Cholesky
# factor has reciprocal condition below 1e-3 (so kappa(x_S'x_S) above about
# 1e6, where the solve could lose more than ~1e-10 relative accuracy) is
# refit by QR instead, as is one whose cross-product is not positive definite.
csr_path <- function(x, y, kmax = min(ncol(x), nrow(x) - 2),
                     max_subsets = 500) {
  p <- ncol(x)
  xm <- colMeans(x)
  xc <- sweep(x, 2, xm)
  G <- crossprod(xc)
  cy <- drop(crossprod(xc, y - mean(y)))
  fit <- function(S) {
    R <- tryCatch(chol(G[S, S, drop = FALSE]), error = function(e) NULL)
    if (is.null(R) || rcond(R, triangular = TRUE) < 1e-3) return(ols_coef(x, y, S, p))
    cf <- backsolve(R, backsolve(R, cy[S], transpose = TRUE))
    b <- numeric(p + 1)
    b[c(1, S + 1)] <- c(mean(y) - sum(xm[S] * cf), cf)
    b
  }
  path <- matrix(0, p + 1, kmax + 1)
  path[1, 1] <- mean(y)
  for (k in seq_len(kmax)) {
    subs <- if (choose(p, k) <= max_subsets) {
      combn(p, k, simplify = FALSE)
    } else {
      replicate(max_subsets, sort.int(sample.int(p, k)), simplify = FALSE)
    }
    path[, k + 1] <- rowMeans(vapply(subs, fit, numeric(p + 1)))
  }
  path
}

# Bayesian model averaging over Occam's window (Raftery 1995; BMA::bicreg,
# defaults). leaps finds the candidate models, so p is capped at maxCol = 31.
bicreg_coef <- function(x, y) {
  fit <- bicreg(x, y, strict = FALSE, OR = 20)
  b <- numeric(ncol(x) + 1)
  b[1] <- fit$postmean[1]
  b[match(fit$namesx, colnames(x)) + 1] <- fit$postmean[-1]
  b
}

relaxed_lasso_path <- function(x, y, gammas = c(0, 0.25, 0.5, 0.75, 1)) {
  fit <- glmnet(x, y, relax = TRUE)
  do.call(cbind, lapply(gammas, function(g) as.matrix(coef(fit, gamma = g))))
}

val_err <- function(path, xval, yval) {
  colMeans((yval - cbind(1, xval) %*% path)^2)
}

pick <- function(path, xval, yval) {
  path[, which.min(val_err(path, xval, yval)), drop = TRUE]
}
