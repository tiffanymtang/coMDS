# Tune hyperparameters for LoCoMDS

Tune hyperparameters for LoCoMDS

## Usage

``` r
tune_locoMDS(obj, data, ks = NULL, make_plots = TRUE)
```

## Arguments

- obj:

  Output of
  [`locoMDS()`](https://tiffanymtang.github.io/coMDS/reference/locoMDS.md).

- data:

  Original data used to evaluate LCMC.

- ks:

  Vector of neighborhood sizes to evaluate LCMC. If `NULL`, defaults to
  a range of sizes between 1 and 70% of the number of samples in `data`.

- make_plots:

  Logical indicating whether to create a plot of LCMC values.

## Value

A list with two components:

- `lcmc`: A data frame containing the LCMC values for each neighborhood
  size and hyperparameter combination.

- `plots`: A ggplot object showing the LCMC values across neighborhood
  sizes and hyperparameters.

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

# tune LoCoMDS hyperparameters
locomds_tuning <- tune_locoMDS(locomds_fits_all, data = X)
```
