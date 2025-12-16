test_that("lcmc works", {
  # generate example data
  n <- 10
  x <- matrix(1:(n * 2), nrow = n, ncol = 2)
  y <- matrix(c(1:3, 8:4, 6:10, 11:8, 20:18), nrow = n, ncol = 2)

  # evaluate LCMC for perfect matching
  ks <- 1:5
  lcmc_out <- lcmc(x, x, ks = ks)
  expect_true(is.data.frame(lcmc_out))
  expect_equal(colnames(lcmc_out), c("k", "lcmc", "adjusted_lcmc"))
  expect_equal(lcmc_out$k, ks)
  expect_equal(lcmc_out$lcmc, rep(1, length(ks)))
  expect_equal(lcmc_out$adjusted_lcmc, 1 - (ks / (n - 1)))

  # evaluate LCMC for random matching
  ks <- 1:5
  lcmc_out <- lcmc(x, y, ks = ks)
  expect_equal(lcmc_out$k, ks)
  expect_equal(lcmc_out$lcmc, c(0.7, 0.8, 0.6, 0.6, 0.7))
  expect_equal(lcmc_out$adjusted_lcmc, lcmc_out$lcmc - (ks / (n - 1)))
})
