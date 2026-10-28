# Plug-in comparison: IRT-based dynamic cutoffs (irtDFI, GPCM + C2 indices)
# vs. an ordinal-CFA DFI procedure (own lavaan/WLSMV re-implementation of the
# logic of McNeish, 2023; NOT the authors' code). Truth = GPCM.
# Jalankan di RStudio: source("compare_cfa.R"); atur R, nrep, cores, out di bawah
# ---- atur di sini ----
R <- 50; nrep <- 50; cores <- 1; out <- "hasil_perbandingan_cfa.rds"
# ----------------------
suppressMessages({library(irtDFI); library(mirt); library(lavaan)})
J <- 10; K <- 4; N <- 500
a <- seq(0.8, 2.2, length.out = J)
set.seed(99)
b <- t(apply(t(sapply(1:J, function(j) seq(-1.5, 1.5, length.out = K - 1) + rnorm(K - 1, 0, .15))), 1, sort))
gen <- function(N, ld = 0, pairs = NULL) {
  z <- rnorm(N); psi <- matrix(0, N, J)
  if (ld > 0) for (p in seq_len(nrow(pairs))) { s <- rnorm(N); psi[, pairs[p, 1]] <- s; psi[, pairs[p, 2]] <- s }
  o <- matrix(0L, N, J)
  for (j in 1:J) {
    th <- z + ld * psi[, j]
    num <- cbind(0, t(apply(outer(th, b[j, ], function(t, bb) a[j] * (t - bb)), 1, cumsum)))
    P <- exp(num); P <- P / rowSums(P)
    o[, j] <- rowSums(runif(N) > t(apply(P, 1, cumsum))) |> pmin(K - 1L)
  }
  d <- as.data.frame(o); names(d) <- paste0("x", 1:J); d
}
mod_syn <- paste("f =~", paste0("x", 1:J, collapse = " + "))
cfa_idx <- function(d) {
  f <- try(suppressWarnings(cfa(mod_syn, d, ordered = names(d), estimator = "WLSMV")), silent = TRUE)
  if (inherits(f, "try-error") || !lavInspect(f, "converged")) return(NULL)
  fm <- try(fitMeasures(f, c("rmsea.scaled", "srmr", "cfi.scaled")), silent = TRUE)
  if (inherits(fm, "try-error")) NULL else list(idx = unname(fm), fit = f)
}
cfa_cut <- function(fit, nrep, N) {
  lam <- standardizedSolution(fit); lam <- lam[lam$op == "=~", "est.std"]
  pe <- parameterEstimates(fit); th <- pe[pe$op == "|", ]
  taus <- lapply(paste0("x", 1:J), function(v) sort(th$est[th$lhs == v]))
  M <- NULL
  for (i in seq_len(nrep)) {
    eta <- rnorm(N)
    d <- as.data.frame(sapply(1:J, function(j) findInterval(lam[j] * eta + sqrt(1 - lam[j]^2) * rnorm(N), taus[[j]])))
    names(d) <- paste0("x", 1:J)
    if (any(sapply(d, function(v) length(unique(v)) < K))) next
    r <- cfa_idx(d); if (!is.null(r)) M <- rbind(M, r$idx)
  }
  c(rmsea = unname(quantile(M[, 1], .95)), srmr = unname(quantile(M[, 2], .95)),
    cfi = unname(quantile(M[, 3], .05)))
}
one <- function(i, ld) {
  set.seed(31000 * (ld > 0) + i)
  d <- gen(N, ld, rbind(c(1, 2), c(3, 4)))
  if (any(sapply(d, function(v) length(unique(v)) < K))) return(NULL)
  m <- suppressMessages(mirt(d, 1, "gpcm", verbose = FALSE))
  r <- dfi_irt(m, nrep = nrep, L = 0, seed = i)
  c0 <- r$cutoffs[r$cutoffs$level == 0, ]; obs <- r$observed
  ci <- cfa_idx(d); if (is.null(ci)) return(NULL)
  cc <- cfa_cut(ci$fit, nrep, N)
  data.frame(rep = i, ld = ld,
    irt_rmsea = obs["RMSEA"] > c0$cutoff[c0$index == "RMSEA"],
    irt_srmsr = obs["SRMSR"] > c0$cutoff[c0$index == "SRMSR"],
    irt_cfi = obs["CFI"] < c0$cutoff[c0$index == "CFI"],
    cfa_rmsea = ci$idx[1] > cc["rmsea"], cfa_srmr = ci$idx[2] > cc["srmr"],
    cfa_cfi = ci$idx[3] < cc["cfi"],
    fix_rmsea = obs["RMSEA"] > .05, fix_srmsr = obs["SRMSR"] > .05, fix_cfi = obs["CFI"] < .95,
    row.names = NULL)
}
res <- list()
for (ld in c(0, 0.7)) {
  rr <- parallel::mclapply(seq_len(R), function(i) one(i, ld), mc.cores = cores)
  res[[as.character(ld)]] <- do.call(rbind, Filter(is.data.frame, rr))
  saveRDS(res, out); cat("done ld", ld, format(Sys.time()), "\n")
}
