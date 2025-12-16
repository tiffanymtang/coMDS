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
