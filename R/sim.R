library(bestsubset)
source("R/methods.R")

# Relative test error, exact for a fresh draw from the simulation design:
# E(y0 - x0'b)^2 / sigma^2. The null model scores 1 + SNR, the oracle 1.
rte <- function(b, d) {
  e <- b[-1] - d$beta
  as.numeric(1 + t(e) %*% d$Sigma %*% e / d$sigma^2)
}

# One replication of the Hastie, Tibshirani & Tibshirani (2020) design:
# training and validation sets of size n; every path is tuned on validation.
one_rep <- function(n, p, rho, snr, beta_type, width = 10, B = 50,
                    kmax = min(p, 50)) {
  d <- sim.xy(n, p, nval = n, rho = rho, s = 5, beta.type = beta_type,
              snr = snr)
  colnames(d$x) <- colnames(d$xval) <- paste0("X", seq_len(p))
  tuned <- function(path) rte(pick(path, d$xval, d$yval), d)

  fs <- beam_path(d$x, d$y, 1, kmax)
  bp <- beam_path(d$x, d$y, width, kmax)
  out <- c(
    fsr = tuned(fs$best),
    beam_best = tuned(bp$best),
    beam_avg = tuned(bp$avg),
    beam_bic = tuned(bp$bic),
    bag_fsr = tuned(bag_fsr(d$x, d$y, B, kmax)),
    csr = tuned(csr_path(d$x, d$y, kmax)),
    relaxed_lasso = tuned(relaxed_lasso_path(d$x, d$y))
  )
  # Both are exhaustive over subsets, so only run where that is feasible.
  if (p <= 30) {
    out <- c(out,
             best_subset = tuned(best_subset_path(d$x, d$y)$path),
             bicreg = rte(bicreg_coef(d$x, d$y), d))
  }
  out
}
