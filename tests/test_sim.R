suppressMessages(source("R/sim.R"))

# A replication's data must not depend on how many cores ran the cell.
g <- data.frame(n = 60, p = 8, rho = 0.35, snr = 1, beta_type = 2)
runs <- lapply(c(1, 2, 4), function(cores) {
  run_cell(g, reps = 5, seed = 7, cores = cores, width = 3, B = 5)
})
stopifnot(identical(runs[[1]], runs[[2]]), identical(runs[[1]], runs[[3]]))
# and the replications are not copies of one another.
stopifnot(!anyDuplicated(runs[[1]]$fsr))

cat("all simulation tests passed\n")
