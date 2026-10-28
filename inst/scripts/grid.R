# Monte Carlo grid (ADEMP). Resumable: one rds per condition.
# Jalankan di RStudio: source("grid.R"). Kondisi yang sudah ada di folder out dilewati.
# ---- atur di sini ----
out <- "hasil_grid"; nrep <- 200; cores <- 1
# ----------------------
suppressMessages({library(irtDFI); library(mirt)})
dir.create(out, showWarnings = FALSE)
conds <- expand.grid(N = c(250, 500, 1000), J = c(8, 12), K = c(3, 5))
for (i in seq_len(nrow(conds))) {
  f <- file.path(out, sprintf("cond_%02d.rds", i)); if (file.exists(f)) next
  N <- conds$N[i]; J <- conds$J[i]; K <- conds$K[i]
  set.seed(5000 + i)
  a <- seq(0.8, 2.2, length.out = J)
  b <- t(sapply(seq_len(J), function(j) seq(-1.5, 1.5, length.out = K - 1) + rnorm(K - 1, 0, 0.15)))
  b <- t(apply(b, 1, sort)); d <- -a * b
  big <- as.data.frame(simdata(a = matrix(a, ncol = 1), d = d, N = 20000, itemtype = "graded"))
  pop <- mirt(big, 1, itemtype = "graded", verbose = FALSE)
  t0 <- Sys.time()
  r <- dfi_irt(pop, nrep = nrep, N = N, L = 3, ld = 0.7, cores = cores, seed = i)
  saveRDS(list(cond = conds[i, ], cutoffs = r$cutoffs, fixed = r$fixed, draws = r$draws,
               secs = as.numeric(difftime(Sys.time(), t0, units = "secs"))), f)
  cat("done", i, format(Sys.time()), "\n")
}
