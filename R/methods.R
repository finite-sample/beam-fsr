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
  xc <- sweep(x, 2, colMeans(x))
  yc <- y - mean(y)
  ss0 <- colSums(xc^2)
  out <- lapply(c(best = 1, avg = 1, bic = 1), function(i) {
    m <- matrix(0, p + 1, kmax + 1)
    m[1, 1] <- mean(y)
    m
  })
  beams <- vector("list", kmax)
  rss <- c(sum(yc^2), rep(NA_real_, kmax))
  beam <- list(integer(0))
  last <- 0
  for (k in seq_len(kmax)) {
    cand <- list()
    score <- numeric(0)
    for (S in beam) {
      if (length(S)) {
        Q <- qr.Q(qr(xc[, S, drop = FALSE]))
        r <- yc - Q %*% crossprod(Q, yc)
        xp <- xc - Q %*% crossprod(Q, xc)
      } else {
        r <- yc
        xp <- xc
      }
      ss <- colSums(xp^2)
      # A column already explained by S (collinear, or a duplicated bootstrap
      # row pattern) adds nothing and would make the fit singular.
      ok <- setdiff(which(ss > tol * ss0), S)
      if (!length(ok)) next
      gain <- drop(crossprod(xp[, ok, drop = FALSE], r))^2 / ss[ok]
      cand <- c(cand, lapply(ok, function(j) sort(c(S, j))))
      score <- c(score, sum(r^2) - gain)
    }
    if (!length(cand)) break
    keep <- !duplicated(vapply(cand, paste, "", collapse = ","))
    cand <- cand[keep]
    score <- score[keep]
    o <- order(score)[seq_len(min(width, length(score)))]
    beam <- cand[o]
    beams[[k]] <- beam
    rss[k + 1] <- score[o[1]]
    fits <- vapply(beam, function(S) ols_coef(x, y, S, p), numeric(p + 1))
    # A saturated fit has RSS ~ 0, and rounding can push it below zero.
    lw <- -n / 2 * log(pmax(score[o], .Machine$double.eps * rss[1]))
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
csr_path <- function(x, y, kmax = min(ncol(x), nrow(x) - 2),
                     max_subsets = 500) {
  p <- ncol(x)
  path <- matrix(0, p + 1, kmax + 1)
  path[1, 1] <- mean(y)
  for (k in seq_len(kmax)) {
    subs <- if (choose(p, k) <= max_subsets) {
      combn(p, k, simplify = FALSE)
    } else {
      replicate(max_subsets, sort(sample.int(p, k)), simplify = FALSE)
    }
    path[, k + 1] <- rowMeans(vapply(subs, function(S) ols_coef(x, y, S, p),
                                     numeric(p + 1)))
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
