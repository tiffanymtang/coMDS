# Extract LoCoMDS fit for particular tau and percentile

Extract LoCoMDS fit for particular tau and percentile

## Usage

``` r
extract_locoMDS(obj, tau, percentile)
```

## Arguments

- obj:

  Output of
  [`locoMDS()`](https://tiffanymtang.github.io/coMDS/reference/locoMDS.md).

- tau:

  Regularization parameter. Can be a scalar or a vector.

- percentile:

  Percentile used to define local neighborhoods. Larger percentiles
  result in larger neighborhoods. Can be a scalar or a vector.

## Value

The locoMDS object corresponding to the provided tau and percentile

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

# fit LoCoMDS with multiple possible hyperparameters
locomds_fits_all <- locoMDS(
  embed_list = dr_list, ndim = 2, tau = c(0.01, 0.1), percentile = c(0.5, 0.8)
)

# extract LoCoMDS with specific hyperparameters
locomds_fit <- extract_locoMDS(locomds_fits_all, tau = 0.1, percentile = 0.5)
```
