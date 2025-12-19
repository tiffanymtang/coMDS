# Local Consensus Multidimensional Scaling (LoCoMDS)

Local Consensus Multidimensional Scaling (LoCoMDS)

## Usage

``` r
locoMDS(
  D_list,
  embed_list = NULL,
  tau,
  percentile,
  ndim = 2,
  constraint = c("indscal", "idioscal", "identity"),
  ...
)
```

## Arguments

- D_list:

  List of distance matrices.

- embed_list:

  Optional list of embeddings to compute distances from. If provided,
  the pairwise Euclidean distances will be computed between the embedded
  samples and used in lieu of `D_list`.

- tau:

  Regularization parameter. Can be a scalar or a vector.

- percentile:

  Percentile used to define local neighborhoods. Larger percentiles
  result in larger neighborhoods. Can be a scalar or a vector.

- ndim:

  Number of dimensions for the output embedding. Default is 2.

- constraint:

  Type of weight constraint to use in the CoMDS algorithm. Options are
  "indscal" (default), "idioscal", or "identity"

- ...:

  Additional arguments passed to the
  [`smacofIndDiff()`](https://tiffanymtang.github.io/coMDS/reference/smacofIndDiff.md)
  function.

## Value

A list containing the following components:

- `gspace`: Consensus embedding matrix of size `n x ndim`, where `n` is
  the number of samples.

- `cweights`: A list of learned weight matrices for each input
  embedding/ distance matrix.

- `stress`: The final stress value of the CoMDS solution.

- `rss`: Residual sum of squares.

- `spp`: Stress per point (in percent).

- `spps`: Stress per point per input embedding/distance matrix (in
  percent, conditional on the input embedding/distance matrix).

- `sps`: Stress per input embedding/distance matrix (in percent).

- `sps_norm`: Stress per input embedding/distance matrix (in percent),
  normalized by sum of source weights.

- `ndim`: Number of dimensions in consensus embedding.

- `model`: Type of smacof model.

- `tau`: The regularization parameter used in the LoCoMDS algorithm.

- `percentile`: The percentile used to define local neighborhoods.

- `niter`: Number of iterations run.

- `nobj`: Number of points/samples.

- `type`: MDS type used (e.g., "interval", "ratio", "ordinal", or
  "mspline").

- `constraint`: Type of weight constraint used in the CoMDS algorithm.

- `call`: The matched call to the function. If multiple hyperparameter
  combinations are provided, the output will be a list of such LoCoMDS
  objects, each corresponding to a different combination of `tau` and
  `percentile`.

## Examples

``` r
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

# fit LoCoMDS using dimension reduction embeddings directly as input
locomds_out <- locoMDS(
  embed_list = dr_list, ndim = 2, tau = 0.1, percentile = 0.5
)
#> Warning: Iteration limit reached! Increase itmax argument!

# fit LoCoMDS using distance matrices as input
dist_list <- purrr::map(dr_list, ~ dist(.x, method = "euclidean"))
locomds_out <- locoMDS(
  D_list = dist_list, ndim = 2, tau = 0.1, percentile = 0.5
)
#> Warning: Iteration limit reached! Increase itmax argument!

# can re-run with larger number of iterations or smaller tolerance if desired
locomds_out <- locoMDS(
  D_list = dist_list, ndim = 2, tau = 0.1, percentile = 0.5,
  itmax = 500, eps = 1e-8
)

# fit LoCoMDS with multiple possible hyperparameters
locomds_out <- locoMDS(
  D_list = dist_list, ndim = 2, tau = c(0.01, 0.1), percentile = c(0.5, 0.8)
)
#> Warning: Iteration limit reached! Increase itmax argument!
#> Warning: Iteration limit reached! Increase itmax argument!
```
