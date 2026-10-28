# ---- internals -----------------------------------------------------------

.fit_it <- function(data, itemtype) {
  suppressMessages(suppressWarnings(
    mirt::mirt(data, 1, itemtype = itemtype, verbose = FALSE)))
}

.item_type <- function(mod) {
  it <- unique(mirt::extract.mirt(mod, "itemtype"))
  if (length(it) != 1L || !it %in% c("graded", "gpcm", "Rasch"))
    stop("The model must be unidimensional with all items 'graded', 'gpcm' or 'Rasch'.",
         call. = FALSE)
  it
}

# Draw responses from the fitted model. Local dependence enters through a
# nuisance factor psi shared by an item pair: effective trait
# sigma * (z + ld * psi) for the pair members.
.gen <- function(mod, N, level, ld, pairs) {
  cf <- mirt::coef(mod, simplify = TRUE)
  sdv <- sqrt(as.numeric(cf$cov[1, 1]))
  J <- length(mirt::extract.mirt(mod, "K"))
  z <- stats::rnorm(N)
  psi <- matrix(0, N, J)
  if (level > 0) for (l in seq_len(level)) {
    p <- pairs[l, ]; s <- stats::rnorm(N)
    psi[, p[1]] <- s; psi[, p[2]] <- s
  }
  out <- matrix(0L, N, J)
  for (j in seq_len(J)) {
    th <- sdv * (z + ld * psi[, j])
    P <- mirt::probtrace(mirt::extract.item(mod, j), matrix(th, ncol = 1))
    u <- stats::runif(N)
    out[, j] <- rowSums(u > t(apply(P, 1, cumsum))) |> pmin(ncol(P) - 1L)
  }
  colnames(out) <- colnames(mirt::extract.mirt(mod, "data"))
  as.data.frame(out)
}

.indices <- function(mod, type = "C2") {
  r <- try(suppressWarnings(mirt::M2(mod, type = type)), silent = TRUE)
  if (inherits(r, "try-error")) return(rep(NA_real_, 4))
  c(RMSEA = r$RMSEA, SRMSR = r$SRMSR, CFI = r$CFI, TLI = r$TLI)
}

.one_rep <- function(mod, itemtype, N, level, ld, pairs, type, categories) {
  d <- .gen(mod, N, level, ld, pairs)
  # keep the original category structure; drop unusable replicates
  if (any(vapply(seq_along(d), function(j) length(unique(d[[j]])) < categories[j], NA)))
    return(rep(NA_real_, 4))
  m <- try(.fit_it(d, itemtype), silent = TRUE)
  if (inherits(m, "try-error")) return(rep(NA_real_, 4))
  .indices(m, type)
}

#' Dynamic fit index cutoffs for polytomous IRT models
#'
#' Derives tailored cutoffs for RMSEA, SRMSR, CFI and TLI (C2 based) from
#' data generated under the fitted model and under increasing local
#' dependence.
#'
#' @param model A fitted unidimensional \code{mirt} model (\code{graded},
#'   \code{gpcm} or \code{Rasch}).
#' @param nrep Replications per level.
#' @param N Sample size of the simulated data sets (default: the sample size of the model's data).
#' @param L Number of misspecification levels (\code{L} item pairs,
#'   one added per level); needs at least \code{2 * L} items.
#' @param ld Strength of the local dependence: the shared nuisance factor
#'   has standard deviation \code{ld} on the standardized latent metric.
#' @param specific Quantile of the Level 0 distribution used as the cutoff
#'   for a correctly specified model.
#' @param sensitiv Share of misspecified replications that should be flagged.
#' @param type \code{"C2"} (default) or \code{"M2*"}.
#' @param seed Integer seed (reproducible, also with \code{cores > 1}).
#' @param cores Number of cores (forked, not on Windows).
#' @return An object of class \code{dfi_irt}.
#' @export
dfi_irt <- function(model, nrep = 200, N = NULL, L = 3, ld = 0.7, specific = 0.95,
                    sensitiv = 0.95, type = "C2", seed = 1, cores = 1L) {
  itemtype <- .item_type(model)
  if (is.null(N)) N <- nrow(mirt::extract.mirt(model, "data"))
  J <- length(mirt::extract.mirt(model, "K"))
  if (J < 2 * L) stop("At least 2 * L items are needed.", call. = FALSE)
  cats <- vapply(seq_len(J), function(j) length(unique(mirt::extract.mirt(model, "data")[, j])), 1L)
  pairs <- if (L > 0) cbind(seq(1, 2 * L, by = 2), seq(2, 2 * L, by = 2)) else matrix(0L, 0, 2)
  obs <- .indices(model, type)
  jobs <- expand.grid(level = 0:L, rep = seq_len(nrep))
  run <- function(i) {
    set.seed(seed * 1000003L + i)
    .one_rep(model, itemtype, N, jobs$level[i], ld, pairs, type, cats)
  }
  res <- if (cores > 1L && .Platform$OS.type != "windows")
    parallel::mclapply(seq_len(nrow(jobs)), run, mc.cores = cores) else lapply(seq_len(nrow(jobs)), run)
  M <- do.call(rbind, lapply(res, function(r) if (inherits(r, "try-error") || is.null(r)) rep(NA_real_, 4) else r))
  colnames(M) <- c("RMSEA", "SRMSR", "CFI", "TLI")
  draws <- data.frame(level = jobs$level, M)
  draws <- draws[stats::complete.cases(draws), ]
  upper <- c(RMSEA = TRUE, SRMSR = TRUE, CFI = FALSE, TLI = FALSE)
  cut <- do.call(rbind, lapply(names(upper), function(ix) {
    v0 <- draws[[ix]][draws$level == 0]
    c0 <- stats::quantile(v0, if (upper[ix]) specific else 1 - specific, names = FALSE)
    rows <- data.frame(index = ix, level = 0, cutoff = c0, detection = NA_real_)
    for (l in seq_len(L)) {
      vl <- draws[[ix]][draws$level == l]
      cl <- stats::quantile(vl, if (upper[ix]) 1 - sensitiv else sensitiv, names = FALSE)
      det <- if (upper[ix]) mean(vl > c0) else mean(vl < c0)
      rows <- rbind(rows, data.frame(index = ix, level = l, cutoff = cl, detection = det))
    }
    rows
  }))
  fixed <- c(RMSEA = 0.05, SRMSR = 0.05, CFI = 0.95, TLI = 0.95)
  flag <- function(ix, v) if (upper[ix]) v > fixed[ix] else v < fixed[ix]
  fx <- do.call(rbind, lapply(names(upper), function(ix) {
    data.frame(index = ix, fixed_rule = fixed[ix],
               reject_true = mean(flag(ix, draws[[ix]][draws$level == 0])),
               reject_level1 = mean(flag(ix, draws[[ix]][draws$level == 1])),
               reject_levelL = mean(flag(ix, draws[[ix]][draws$level == L])),
               row.names = NULL)
  }))
  structure(list(observed = obs, cutoffs = cut, fixed = fx, draws = draws,
                 call = list(N = N, J = J, nrep = nrep, L = L, ld = ld,
                             itemtype = itemtype, n_usable = nrow(draws), type = type)),
            class = "dfi_irt")
}

#' @export
print.dfi_irt <- function(x, digits = 3, ...) {
  cat(sprintf("Dynamic fit cutoffs (%s), N = %d, J = %d, %d usable of %d replications\n",
              x$call$itemtype, x$call$N, x$call$J, x$call$n_usable,
              (x$call$L + 1) * x$call$nrep))
  cat("\nObserved fit:\n"); print(round(x$observed, digits))
  cat("\nCutoffs (Level 0 = correct model; detection = share of misspecified draws flagged by the Level 0 cutoff):\n")
  d <- x$cutoffs; d$cutoff <- round(d$cutoff, digits); d$detection <- round(d$detection, 2)
  print(d, row.names = FALSE)
  cat("\nFixed rules of thumb (share of replications rejected):\n")
  f <- x$fixed; f[-1] <- round(f[-1], 2); print(f, row.names = FALSE)
  invisible(x)
}

#' Plot the simulated fit distributions
#' @param x A \code{dfi_irt} object.
#' @param index One of RMSEA, SRMSR, CFI, TLI.
#' @return A ggplot object.
#' @export
dfi_plot <- function(x, index = c("RMSEA", "SRMSR", "CFI", "TLI")) {
  index <- match.arg(index)
  d <- x$draws; d$Level <- factor(d$level)
  ggplot2::ggplot(d, ggplot2::aes(.data[[index]], fill = Level)) +
    ggplot2::geom_density(alpha = 0.4) +
    ggplot2::geom_vline(xintercept = x$observed[[index]], linetype = 2) +
    ggplot2::labs(x = index, y = "Density") + ggplot2::theme_bw()
}
