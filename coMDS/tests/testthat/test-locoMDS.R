test_that("locoMDS works", {
  set.seed(123)

  # set up data
  tau <- 0.1
  taus <- c(0.1, 0.5)
  percentile <- 0.5
  percentiles <- c(0.3, 0.5, 0.9)
  itmax <- 500
  data(iris)
  X <- iris[, 1:4] |>
    dplyr::distinct()
  n <- nrow(X)
  pca_scores <- prcomp(X, center = TRUE, scale = TRUE)$x
  tsne_scores <- Rtsne::Rtsne(X, dims = 2, perplexity = 30, verbose = FALSE)$Y
  umap_scores <- umap::umap(X, n_components = 2, verbose = FALSE)$layout
  embed_list <- list(
    pca = pca_scores,
    tsne = tsne_scores,
    umap = umap_scores
  )
  D_list <- list(
    pca = dist(pca_scores),
    tsne = dist(tsne_scores),
    umap = dist(umap_scores)
  )

  # run locoMDS with single hyperparameter
  locomds_out <- locoMDS(
    D_list, tau = tau, percentile = percentile, ndim = 2, itmax = itmax
  )
  # check output
  expected_names <- c(
    "gspace", "cweights", "stress", "rss", "spp", "spps", "sps", "sps_norm",
    "ndim", "model", "tau", "percentile", "niter", "nobj", "type", "constraint",
    "call"
  )
  expect_true(is.list(locomds_out))
  expect_s3_class(locomds_out, "locoMDS")
  expect_equal(names(locomds_out), expected_names)
  # expect_snapshot(locomds_out)

  # run locoMDS with single hyperparameter and embedding list
  locomds_out_embed <- locoMDS(
    embed_list = embed_list, tau = tau, percentile = percentile, ndim = 2,
    itmax = itmax
  )
  # check output
  expect_equal(locomds_out, locomds_out_embed)

  # run locoMDS with multiple hyperparameters
  locomds_tuned_out <- locoMDS(
    D_list, tau = taus, percentile = percentiles, ndim = 2, itmax = itmax
  )
  # check output
  expect_length(locomds_tuned_out, length(taus) * length(percentiles))
  expect_equal(
    locomds_out,
    extract_locoMDS(locomds_tuned_out, tau = tau, percentile = percentile)
  )

  # check tune_locoMDS
  ks <- c(1, 2, 5, round(10^seq(1, log10(nrow(X) * 0.7), length.out = 7)))
  expect_warning(tune_locoMDS(locomds_out, data = X, ks = ks))
  tuned_out <- suppressWarnings(
    tune_locoMDS(locomds_tuned_out, data = X, ks = ks)
  )
  expect_equal(names(tuned_out), c("lcmc", "plot"))
  expect_equal(
    nrow(tuned_out$lcmc),
    length(ks) * length(taus) * length(percentiles)
  )
  expect_equal(
    colnames(tuned_out$lcmc),
    c("k", "tau", "percentile", "lcmc", "adjusted_lcmc")
  )
  # expect_snapshot(tuned_out$lcmc)
  expect_s3_class(tuned_out$plot, "gg")
})
