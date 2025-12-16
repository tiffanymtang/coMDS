#' Generate data from mixture of 3 Gaussians
#'
#' @returns A list with elements:
#' - X: A data frame containing the generated covariates data
#' - labels: A vector containing the cluster labels (1, 2, or 3)
mixture_of_gaussians_dgp <- function() {
  k <- 3
  cluster_sizes <- c(150, 200, 250)

  # cluster means
  mus <- list(c(-3, -2, 0), c(2, -4, 1), c(0, 6, 6))

  # cluster covariances
  sig <- matrix(c(1, 0.3, 0.2, 0.3, 1, 0.4, 0.2, 0.4, 1), ncol = k)
  sd <- 0.5
  Z1 <- matrix(rnorm(k * k, 0, sd), nrow = k)
  Z2 <- matrix(rnorm(k * k, 0, sd), nrow = k)
  Z3 <- matrix(rnorm(k * k, 0, sd), nrow = k)
  sigs <- list(
    sig + Z1 %*% t(Z1),
    sig + Z2 %*% t(Z2),
    sig + Z3 %*% t(Z3)
  )

  # generate data
  X <- purrr::pmap(
    list(cluster_size = cluster_sizes, mu = mus, sig = sigs),
    function(cluster_size, mu, sig) {
      as.data.frame(MASS::mvrnorm(n = cluster_size, mu = mu, Sigma = sig))
    }
  ) |>
    purrr::list_rbind() |>
    scale()
  labels <- as.factor(rep(1:k, times = cluster_sizes))
  return(list(X = X, labels = labels))
}


#' Generate data from swiss roll
#'
#' @param n_samples Number of samples to generate
#' @param gap Gap between the turns of the swiss roll
#' @param width Width of the swiss roll
#' @param t_min Minimum value of the parameter t
#' @param t_max Maximum value of the parameter t
#' @param noise Standard deviation of Gaussian noise to add to the data
#'
#' @returns A list with elements:
#' - X: A data frame containing the covariates data
#' - labels: A vector containing the parameter t used to generate the swiss roll
swiss_roll_dgp <- function(n_samples = 1000, gap = 1, width = 1,
                           t_min = 1 * pi, t_max = 3.5 * pi, noise = 0.1) {
  # generate swiss roll data
  t <- runif(n_samples, t_min, t_max)
  z <- runif(n_samples, -width, width)
  x <- gap * t * cos(t)
  y <- gap * t * sin(t)

  # add Gaussian noise
  X <- cbind(x, y, z) + matrix(rnorm(n_samples * 3, sd = noise), ncol = 3)
  X <- scale(X)

  return(list(X = X, labels = t))
}


#' Load real data from a given path
#'
#' @param path Path to the .RData file containing the data
#'
#' @returns A list with elements:
#' - X: A data frame containing the covariates data
#' - labels: A vector containing the cluster labels
real_data_dgp <- function(path) {
  load(path)
  return(list(X = data, labels = info))
}
