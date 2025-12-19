# Check consensus embeddings for uniqueness of points

Check consensus embeddings for uniqueness of points

## Usage

``` r
check_consensus_embeddings(y, fit = NULL)
```

## Arguments

- y:

  Raw data/embeddings in new (dimension-reduced) space.

- fit:

  (Optional) Fitted `locoMDS` object corresponding to the embeddings in
  `y`. Only used for warning messages about uniqueness of points.

## Value

Original consensus embeddings if sufficiently unique; otherwise, the
rounded consensus embeddings with a warning.
