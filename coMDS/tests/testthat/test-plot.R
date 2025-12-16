test_that("plot_coMDS works", {
  skip()
  skip_on_cran()
  skip_on_ci()

  set.seed(123)

  # set up data
  tau <- 0.1
  taus <- c(0.1, 0.5)
  percentile <- 0.5
  percentiles <- c(0.3, 0.5, 0.9)
  itmax <- 500
  data(iris)
  iris <- dplyr::distinct(iris)
  X <- iris[, 1:4]
  y <- iris$Species
  n <- nrow(X)
  pca_scores <- prcomp(X, center = TRUE, scale = TRUE)$x
  tsne_scores <- Rtsne::Rtsne(X, dims = 2, perplexity = 30, verbose = FALSE)$Y
  umap_scores <- umap::umap(X, n_components = 2, verbose = FALSE)$layout
  D_list <- list(
    pca = dist(pca_scores),
    tsne = dist(tsne_scores),
    umap = dist(umap_scores)
  )

  # run coMDS
  comds_out <- coMDS(D_list, ndim = 2, itmax = itmax)

  # run locoMDS with single hyperparameter
  locomds_out <- locoMDS(
    D_list, tau = tau, percentile = percentile, ndim = 2, itmax = itmax
  )

  # run locoMDS with multiple hyperparameters
  locomds_tuned_out <- locoMDS(
    D_list, tau = taus, percentile = percentiles, ndim = 2, itmax = itmax
  )

  fit_ls <- list(
    "CoMDS" = comds_out,
    "locoMDS_single" = locomds_out,
    "locoMDS_multi" = locomds_tuned_out
  )
  for (fit_name in names(fit_ls)) {
    fit <- fit_ls[[fit_name]]
    for (type in c("scores", "relative_errors")) {
      vdiffr::expect_doppelganger(
        sprintf("plot-%s-%s", fit_name, type),
        plot_coMDS(fit, color = y, type = type)
      )
    }
  }
})
