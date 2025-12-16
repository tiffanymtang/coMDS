#' Consensus Multidimensional Scaling (CoMDS)
#'
#' @param D_list List of distance matrices.
#' @param embed_list Optional list of embeddings to compute distances from. If
#'   provided, the pairwise Euclidean distances will be computed between the
#'   embedded samples and used in lieu of `D_list`.
#' @param ndim Number of dimensions for the output embedding. Default is 2.
#' @param constraint Type of weight constraint to use in the CoMDS algorithm.
#'   Options are "indscal" (default), "idioscal", or "identity"
#' @param ... Additional arguments passed to the [smacofIndDiff()] function.
#'
#' @returns A list containing the following components:
#' - `gspace`: Consensus embedding matrix of size `n x ndim`, where `n` is the
#'   number of samples.
#' - `cweights`: A list of learned weight matrices for each input embedding/
#'   distance matrix.
#' - `stress`: The final stress value of the CoMDS solution.
#' - `spp`: Stress per point (in percent).
#' - `rss`: Residual sum of squares.
#' - `spps`: Stress per point per input embedding/distance matrix (in percent,
#'   conditional on the input embedding/distance matrix).
#' - `sps`: Stress per input embedding/distance matrix (in percent).
#' - `sps_norm`: Stress per input embedding/distance matrix (in percent),
#'   normalized by sum of source weights.
#' - `ndim`: Number of dimensions in consensus embedding.
#' - `model`: Type of smacof model.
#' - `niter`: Number of iterations run.
#' - `nobj`: Number of points/samples.
#' - `type`: MDS type used (e.g., "interval", "ratio", "ordinal", or "mspline").
#' - `constraint`: Type of weight constraint used in the CoMDS algorithm.
#' - `call`: The matched call to the function.
#'
#' @examples
#' data(iris)
#' # remove duplicates so that tSNE can run
#' X <- dplyr::distinct(iris[, 1:4])
#'
#' # fit various dimension reduction methods
#' pca_scores <- prcomp(X, center = TRUE, scale = TRUE)$x
#' tsne_scores <- Rtsne::Rtsne(X, dims = 2, perplexity = 30, verbose = FALSE)$Y
#' umap_scores <- umap::umap(X, n_components = 2, verbose = FALSE)$layout
#' dr_list <- list(
#'   pca = pca_scores,
#'   tsne = tsne_scores,
#'   umap = umap_scores
#' )
#'
#' # fit coMDS using dimension reduction embeddings directly as input
#' comds_out <- coMDS(embed_list = dr_list, ndim = 2)
#'
#' # fit coMDS using distance matrices as input
#' dist_list <- purrr::map(dr_list, ~ dist(.x, method = "euclidean"))
#' comds_out <- coMDS(D_list = dist_list, ndim = 2)
#'
#' # can re-run with larger number of iterations or smaller tolerance if desired
#' comds_out <- coMDS(D_list = dist_list, ndim = 2, itmax = 500, eps = 1e-8)
#'
#' @export
coMDS <- function(D_list, embed_list = NULL, ndim = 2,
                  constraint = c("indscal", "idioscal", "identity"),
                  ...) {
  # setup
  constraint <- match.arg(constraint)
  if (!is.null(embed_list)) {
    D_list <- purrr::map(embed_list, ~ dist(.x, method = "euclidean"))
  }

  # normalize distances for numerical stability
  D_list <- normalize_distances(D_list)

  # run coMDS
  result <- smacofIndDiff(
    D_list, ndim = ndim, constraint = constraint, ...
  )
  return(result)
}
