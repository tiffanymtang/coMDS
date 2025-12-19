# Evaluate Local Continuity Meta-Criterion (LCMC)

Evaluate Local Continuity Meta-Criterion (LCMC)

## Usage

``` r
lcmc(x, y, ks, fit = NULL)
```

## Arguments

- x:

  Raw data/embeddings in original space.

- y:

  Raw data/embeddings in new (dimension-reduced) space.

- ks:

  Vector of neighborhood sizes to evaluate LCMC.

- fit:

  (Optional) Fitted `locoMDS` object corresponding to the embeddings in
  `y`. Only used for warning messages about uniqueness of points.

## Value

A data frame with columns:

- `k`: Neighborhood size.

- `lcmc`: Local Continuity Meta-Criterion value for the given
  neighborhood size.

- `adjusted_lcmc`: Adjusted LCMC value, accounting for the neighborhood
  size.

## Examples

``` r
# generate example data
n <- 100
p1 <- 4
p2 <- 2
x <- matrix(rnorm(n * p1), nrow = n, ncol = p1)
y <- matrix(rnorm(n * p2), nrow = n, ncol = p2)

# evaluate LCMC for neighborhood sizes 5, 10, and 20
lcmc_out <- lcmc(x, y, ks = c(5, 10, 20))
```
