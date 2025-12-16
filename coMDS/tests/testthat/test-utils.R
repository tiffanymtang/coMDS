test_that("normalize_distances works", {
  set.seed(123)
  data(iris)
  # remove duplicates so that tSNE can run
  X <- dplyr::distinct(iris[, 1:4])

  # fit various dimension reduction methods
  pca_scores <- prcomp(X, center = TRUE, scale = TRUE)$x
  tsne_scores <- Rtsne::Rtsne(X, dims = 2, perplexity = 30, verbose = FALSE)$Y
  umap_scores <- umap::umap(X, n_components = 2, verbose = FALSE)$layout
  dr_list <- list(
    pca = pca_scores,
    tsne = tsne_scores,
    umap = umap_scores
  )

  # convert embeddings to distance matrices
  dist_list <- purrr::map(dr_list, ~ dist(.x, method = "euclidean"))
  # normalize distance matrices
  dist_list_norm <- normalize_distances(dist_list)
  # check that output is a list of dist objects
  expect_true(all(
    purrr::map_lgl(
      dist_list_norm,
      function(D) {
        inherits(D, "dist")
      }
    )
  ))
  # check that top singular value of each normalized distance matrix is 1
  expect_true(all(
    purrr::map_lgl(
      dist_list_norm,
      function(D) {
        D_mat <- as.matrix(D)
        if (any(is.na(D_mat))) {
          diag(D_mat) <- NA
          na_rows <- which(apply(D_mat, 1, function(x) all(is.na(x))))
          D_mat <- D_mat[-na_rows, -na_rows]
          diag(D_mat) <- 0
        }
        D_svd <- svd(D_mat, nu = 1, nv = 1)
        return(abs(D_svd$d[1] - 1) < 1e-6)
      }
    )
  ))
})
