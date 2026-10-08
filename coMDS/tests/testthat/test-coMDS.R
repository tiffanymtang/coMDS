test_that("coMDS works", {
  # set up data
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
  D_list_norm <- normalize_distances(D_list)

  # run coMDS
  comds_out <- coMDS(D_list, ndim = 2)
  # check output
  expected_names <- c(
    "gspace", "cweights", "stress", "rss", "spp", "spps", "sps", "sps_norm",
    "ndim", "model", "niter", "nobj", "type", "constraint", "call"
  )
  expect_true(is.list(comds_out))
  expect_s3_class(comds_out, "coMDS")
  expect_equal(names(comds_out), expected_names)
  truth <- smacof::smacofIndDiff(
    D_list_norm, ndim = 2,
    itmax = 100, eps = 1e-6 * length(D_list_norm) * n * (n - 1) / 2
  )
  for (name in setdiff(expected_names, c("model", "call", "sps_norm"))) {
    if (name == "cweights") {
      expect_equal(
        purrr::map(comds_out[[name]], ~ unname(.x)),
        purrr::map(truth[[name]], ~ unname(.x))
      )
    } else {
      expect_equal(
        unname(comds_out[[name]]),
        unname(truth[[name]])
      )
    }
  }

  # run coMDS with embedding list
  comds_out_embed <- coMDS(embed_list = embed_list, ndim = 2)
  # check output
  expect_equal(comds_out, comds_out_embed)

  # check coMDS with ndim = 1
  expect_error(coMDS(D_list, ndim = 1), NA)
  expect_error(coMDS(D_list, ndim = 1, constraint = "identity"), NA)
  expect_error(coMDS(D_list, ndim = 1, constraint = "indscal"), NA)

  # check with weights
  comds_out_wt1 <- coMDS(D_list, ndim = 2, weights = c(1, 2, 20), itmax = 500)
  comds_out_wt2 <- coMDS(D_list, ndim = 2, weights = c(1, 1, 1))
  expect_true(comds_out_wt1$sps["umap"] < comds_out_wt2$sps["umap"])

  # check with adding duplicate distance matrix
  comds_out_wt3 <- coMDS(D_list, ndim = 2, weights = c(2, 1, 1))
  comds_out_wt4 <- coMDS(
    c(list(pca_copy = dist(pca_scores)), D_list),
    ndim = 2, weights = c(1, 1, 1, 1)
  )
  expect_equal(
    comds_out_wt3$sps[[1]],
    sum(comds_out_wt4$sps[1:2]),
    tolerance = 1e-1
  )
  expect_equal(
    abs(cor(comds_out_wt3$gspace[, 1], comds_out_wt4$gspace[, 1])),
    1,
    tolerance = 1e-4
  )
  expect_equal(
    abs(cor(comds_out_wt3$gspace[, 2], comds_out_wt4$gspace[, 2])),
    1,
    tolerance = 1e-4
  )
})

test_that("coMDS source weights are invariant to a common scale", {
  set.seed(123)
  X <- matrix(rnorm(30), nrow = 10)
  D_list <- list(
    all_dimensions = dist(X),
    selected_dimensions = dist(X[, c(1, 3)]),
    squared_values = dist(X^2)
  )
  init <- cmdscale(D_list[[1]], k = 2)
  weights <- c(1, 0.5, 0.5)

  fit <- coMDS(
    D_list, ndim = 2, constraint = "identity", weights = weights,
    init = init, eps = 1e-6, itmax = 1000
  )
  scaled_fit <- coMDS(
    D_list, ndim = 2, constraint = "identity", weights = 10 * weights,
    init = init, eps = 1e-6, itmax = 1000
  )

  expect_equal(scaled_fit$stress, fit$stress, tolerance = 1e-10)
  expect_equal(scaled_fit$niter, fit$niter)
  expect_equal(scaled_fit$gspace, fit$gspace, tolerance = 1e-10)
})

test_that("coMDS half-weight duplicates equal one full-weight source", {
  set.seed(321)
  X <- matrix(rnorm(30), nrow = 10)
  A <- dist(X)
  B <- dist(cbind(X[, 1] + X[, 2], X[, 3]^2))

  fit <- coMDS(
    list(A = A, B = B), weights = c(1, 1),
    eps = 1e-8, itmax = 1000
  )
  duplicate_fit <- coMDS(
    list(A = A, B1 = B, B2 = B),
    weights = c(1, 0.5, 0.5), eps = 1e-8, itmax = 1000
  )

  expect_equal(
    as.vector(dist(duplicate_fit$gspace)),
    as.vector(dist(fit$gspace)), tolerance = 1e-6
  )
  expect_equal(duplicate_fit$stress, fit$stress, tolerance = 1e-8)
  expect_equal(duplicate_fit$rss, fit$rss, tolerance = 1e-6)
  expect_equal(duplicate_fit$spp, fit$spp, tolerance = 1e-6)
  expect_equal(
    unname(c(duplicate_fit$sps[1], sum(duplicate_fit$sps[2:3]))),
    unname(fit$sps), tolerance = 1e-6
  )
  expect_equal(
    unname(c(duplicate_fit$sps_norm[1], sum(duplicate_fit$sps_norm[2:3]))),
    unname(fit$sps_norm), tolerance = 1e-6
  )
})
