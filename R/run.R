# Usage: Rscript R/run.R <reps> <outdir> [p]
# Each grid cell is saved on its own and finished cells are skipped, so an
# interrupted run loses only the cells in flight.
suppressMessages(source("R/sim.R"))
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
reps <- as.integer(args[1])
outdir <- args[2]
p <- if (length(args) > 2) as.integer(args[3]) else 20
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

grid <- expand.grid(n = 100, p = p, rho = c(0, 0.35, 0.7),
                    snr = c(0.05, 0.25, 1, 6), beta_type = c(1, 2))

RNGkind("L'Ecuyer-CMRG")
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  f <- file.path(outdir, sprintf("cell_p%d_rho%s_snr%s_beta%d.rds",
                                 g$p, g$rho, g$snr, g$beta_type))
  if (file.exists(f)) next
  set.seed(1000 * p + i)
  res <- mclapply(seq_len(reps), function(r) {
    suppressWarnings(one_rep(g$n, g$p, g$rho, g$snr, g$beta_type))
  }, mc.cores = detectCores())
  saveRDS(cbind(g, rep = seq_len(reps), do.call(rbind, res)), f)
  message(sprintf("%d/%d %s", i, nrow(grid), basename(f)))
}
