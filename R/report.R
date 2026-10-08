# Usage: Rscript R/report.R
# Reads results/raw, writes results/summary.csv and figures/rte.png, and
# rewrites the README block between the results markers, so every number in
# the README comes from the saved simulations.
suppressMessages(library(ggplot2))

# Cells at p > 30 have no best_subset or bicreg columns; fill them with NA.
cells <- lapply(list.files("results/raw", full.names = TRUE), readRDS)
cols <- unique(unlist(lapply(cells, names)))
res <- do.call(rbind, lapply(cells, function(d) {
  d[setdiff(cols, names(d))] <- NA
  d[cols]
}))
methods <- c(beam_best = "Beam, best member", beam_avg = "Beam, average",
             beam_bic = "Beam, BIC-weighted average",
             bag_fsr = "Bagged forward stepwise", csr = "Complete subset regression",
             relaxed_lasso = "Relaxed lasso", best_subset = "Best subset",
             bicreg = "BMA (bicreg)")
long <- do.call(rbind, lapply(intersect(names(methods), names(res)), function(m) {
  data.frame(p = res$p, snr = res$snr, rho = res$rho, beta_type = res$beta_type,
             method = m, diff = res[[m]] - res$fsr)
}))
long <- long[!is.na(long$diff), ]

summ <- aggregate(diff ~ p + snr + method, long, function(z) {
  c(mean = mean(z), half = 1.96 * sd(z) / sqrt(length(z)), n = length(z))
})
summ <- cbind(summ[1:3], as.data.frame(summ$diff))
dir.create("figures", showWarnings = FALSE)
write.csv(summ, "results/summary.csv", row.names = FALSE)

fmt_table <- function(pp) {
  s <- summ[summ$p == pp, ]
  snrs <- sort(unique(s$snr))
  head <- paste0("| Method | ", paste0("SNR ", snrs, collapse = " | "), " |")
  rule <- paste0("|---|", strrep("---|", length(snrs)))
  rows <- vapply(intersect(names(methods), unique(s$method)), function(m) {
    cells <- vapply(snrs, function(v) {
      r <- s[s$method == m & s$snr == v, ]
      sprintf("%+.3f ± %.3f", r$mean, r$half)
    }, "")
    paste0("| ", methods[[m]], " | ", paste(cells, collapse = " | "), " |")
  }, "")
  r <- res[res$p == pp, ]
  title <- sprintf(paste("**p = %d** (n = 100; %d replications per SNR, pooled over",
                         "correlations %s and coefficient pattern%s %s)"),
                   pp, max(s$n), paste(sort(unique(r$rho)), collapse = ", "),
                   if (length(unique(r$beta_type)) > 1) "s" else "",
                   paste(sort(unique(r$beta_type)), collapse = " and "))
  c(title, "", head, rule, rows)
}

block <- c("<!-- results:start -->",
           "Relative test error minus forward stepwise's, paired by replication",
           "(mean ± 95% CI half-width; negative = better than forward stepwise).",
           "")
for (pp in sort(unique(summ$p))) block <- c(block, fmt_table(pp), "")
block <- c(block, "![Relative test error vs forward stepwise](figures/rte.png)",
           "<!-- results:end -->")

readme <- readLines("README.md")
a <- grep("<!-- results:start -->", readme, fixed = TRUE)
b <- grep("<!-- results:end -->", readme, fixed = TRUE)
writeLines(c(readme[seq_len(a - 1)], block, readme[-seq_len(b)]), "README.md")

summ$label <- factor(methods[summ$method], levels = methods)
g <- ggplot(summ, aes(factor(snr), mean, colour = label, group = label)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_pointrange(aes(ymin = mean - half, ymax = mean + half),
                  position = position_dodge(width = 0.6), size = 0.25) +
  facet_wrap(~ p, labeller = label_both) +
  # Complete subset regression reaches +1.9 at high SNR and would flatten
  # every other method; the table carries its full values.
  coord_cartesian(ylim = c(-0.08, 0.12)) +
  labs(x = "Signal-to-noise ratio", colour = NULL,
       y = "Relative test error minus forward stepwise's",
       caption = "Complete subset regression is off scale at SNR 1 and 6; see the table.") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 3))
ggsave("figures/rte.png", g, width = 9, height = 5.5, dpi = 150, bg = "white")
cat("wrote results/summary.csv, figures/rte.png, README results block\n")
