#' Faster version of `smacof::myGenInv`
#' @keywords internal
myGenInv <- function(x) {
  n <- dim(x)[1]
  nn <- 1 / n
  # return(solve(x + nn) - nn)
  return(chol2inv(chol(x + nn)) - nn)
}


#' Normalize dissimilarities for CoMDS, allowing different sources to have
#' different overall weights
#' @keywords internal
normDissN <- function(diss, wghts, m) {
  N <- length(diss) * m
  dissnorm <- diss / sqrt(sum(wghts * diss^2, na.rm = TRUE)) * sqrt(N)
  return(dissnorm)
}


#' Regularized version of `smacof::vmat()` to avoid singularity issues
#' @keywords internal
vmat <- function(wgths, weightlam = 0) {
  v <- as.matrix(wgths)
  r <- rowSums(v) # row margins of weight matrix + regularization
  return(diag(r) - v + diag(weightlam, nrow(v)))
}


#' Copy of `smacof:::strucprep` since not exported from `smacof`
#' @keywords internal
strucprep <- function(x) {
  distvec <- as.vector(x[lower.tri(x)])
  n <- dim(x)[1]
  dissim <- structure(
    distvec,
    Size = n, call = quote(as.dist.default(m = b)), class = "dist",
    Diag = FALSE, Upper = FALSE
  )
}


#' Copy of `smacof:::checkdiss` since not exported from `smacof`
#' @keywords internal
checkdiss <- function(diss) {
  if (any(sapply(diss, function(d0) any(d0 < 0, na.rm = TRUE)))) {
    stop("Dissimilarities should be non-negative!")
  }
}


#' Copy of `smacof:::initWeights` since not exported from `smacof`
#' @keywords internal
initWeights <- function(diss) {
  if (!is.list(diss)) {
    n <- attr(diss, "Size")
    ww <- matrix(1, n, n)
    ww[is.na(as.matrix(diss))] <- 0 ## blank out missings
    return(as.dist(ww))
  } else {
    n <- attr(diss[[1]], "Size")
    m <- length(diss)
    ww <- repList(matrix(1, n, n), m)
    for (i in 1:m) {
      wwi <- ww[[i]]
      wwi[is.na(as.matrix(diss[[i]]))] <- 0
      ww[[i]] <- as.dist(wwi)
    }
    return(ww)
  }
}


#' Copy of `smacof:::initConf` since not exported from `smacof`
#' @keywords internal
initConf <- function(init, diss, n, p, inddiff = FALSE) {
  if (inddiff) diss <- as.dist(apply(simplify2array(lapply(diss, as.matrix)), c(1, 2), sum, na.rm = TRUE))

  if (length(init) == 1) {
    if (init == "torgerson") {
      meandiss <- mean(diss, na.rm = TRUE) ## mean dissimilarity for imputation
      diss1 <- as.matrix(diss)
      diss1[is.na(diss1)] <- meandiss
      x <- smacof::torgerson(diss1, p = p)
      init <- "dummy"
    }
    if (init == "random") {
      x <- matrix(runif(n * p, min = -1), ncol = p)
    }
  }
  if (is.data.frame(init)) init <- as.matrix(init)
  if (is.matrix(init)) {
    x <- as.matrix(init)
    if (any(dim(x) != c(n, p))) stop(paste0("Dimension of the starting configuration matrix needs to be ", n, " times ", p, "!"))
  }
  return(x)
}


#' Copy of `smacof:::bmat` since not exported from `smacof`
#' @keywords internal
bmat <- function(diss, wgths, d, eps = 1E-12) {
  z <- ifelse(d < eps, 1, 0)
  b <- as.matrix((wgths * diss * (1 - z)) / (d + z))
  r <- rowSums(b)
  return(diag(r) - b)
}


#' Copy of `smacof:::spp` since not exported from `smacof`
#' @keywords internal
spp <- function(dhat, confdiss, wgths) {
  resmat <- as.matrix(wgths) * as.matrix(dhat - confdiss)^2 # point stress
  diag(resmat) <- NA
  spp <- colMeans(resmat, na.rm = TRUE)
  spp <- spp / sum(spp) * 100
  names(spp) <- colnames(resmat) <- rownames(resmat) <- attr(dhat, "Labels")
  return(list(spp = spp, resmat = resmat))
}


#' Copy of `smacof:::appendList` since not exported from `smacof`
#' @keywords internal
appendList <- function(x, a) {
  return(c(x, list(a)))
}


#' Copy of `smacof:::repList` since not exported from `smacof`
#' @keywords internal
repList <- function(x, n) {
  z <- list()
  for (i in 1:n) {
    z <- c(z, list(x))
  }
  return(z)
}
