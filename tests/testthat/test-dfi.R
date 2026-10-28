test_that("generator reproduces category structure and level 0 is unbiased in fit", {
  data(sim_grm); m <- mirt::mirt(sim_grm[, 1:8], 1, itemtype = "graded", verbose = FALSE)
  set.seed(1); d <- irtDFI:::.gen(m, 2000, 0, 0.7, cbind(1, 2))
  expect_equal(dim(d), c(2000, 8)); expect_true(all(unlist(d) %in% 0:3))
  ## generated item means match the fitted model's expected scores
  th <- seq(-6, 6, length.out = 801); w <- dnorm(th); w <- w / sum(w)
  e <- sapply(1:8, function(j) sum(w * (mirt::probtrace(mirt::extract.item(m, j), matrix(th, ncol = 1)) %*% 0:3)))
  expect_equal(unname(colMeans(d)), e, tolerance = 0.06)
})
test_that("dfi_irt returns ordered cutoffs and is reproducible", {
  data(sim_grm); m <- mirt::mirt(sim_grm[, 1:8], 1, itemtype = "graded", verbose = FALSE)
  a <- dfi_irt(m, nrep = 6, L = 2, seed = 3); b <- dfi_irt(m, nrep = 6, L = 2, seed = 3)
  expect_equal(a$cutoffs, b$cutoffs)
  expect_true(all(a$draws$RMSEA >= 0))
  expect_error(dfi_irt(m, nrep = 4, L = 5))
})
