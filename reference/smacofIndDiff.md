# SMACOF for Individual Differences

Variant of
[`smacof::smacofIndDiff()`](https://rdrr.io/pkg/smacof/man/smacofIndDiff.html)
but with regularization (`weightlam`) in case of singularity in weight
matrices (e.g., due to missingness) and optional weight scaling
(`weights`) which multiples weights for each data source by some amount.

## Usage

``` r
smacofIndDiff(
  delta,
  ndim = 2,
  type = c("ratio", "interval", "ordinal", "mspline"),
  constraint = c("indscal", "idioscal", "identity"),
  weights = NULL,
  weightmat = NULL,
  weightlam = 0,
  init = "torgerson",
  ties = "primary",
  verbose = FALSE,
  modulus = 1,
  itmax = 100,
  eps = 1e-06,
  spline.degree = 2,
  spline.intKnots = 2
)
```

## Arguments

- delta:

  A list of dissimilarity matrices or a list objects of class `dist`

- ndim:

  Number of dimensions

- type:

  MDS type: `"interval"`, `"ratio"`, `"ordinal"` (nonmetric MDS), or
  `"mspline"`

- constraint:

  Either `"indscal"`, `"idioscal"`, or `"identity"` (see details)

- weights:

  Vector of length equal to the number of distance matrices. If
  provided, weights for each distance matrix are multiplied by the
  corresponding value in `weights`. Default is `NULL` which will give
  each distance matrix the same weight (i.e., scaled by 1).

- weightmat:

  Optional matrix with dissimilarity weights

- weightlam:

  Regularization parameter for weight matrices to avoid singularity
  matrices. Default is 0 (i.e., no regularization). If regularization is
  needed, set to a small value (e.g., 1e-6).

- init:

  Matrix with starting values for configurations (optional)

- ties:

  Tie specification for non-metric MDS

- verbose:

  If `TRUE`, intermediate stress is printed out

- modulus:

  Number of smacof iterations per monotone regression call

- itmax:

  Maximum number of iterations

- eps:

  Convergence criterion

- spline.degree:

  Degree of the spline for `"mspline"` MDS type

- spline.intKnots:

  Number of interior knots of the spline for `"mspline"` MDS type
