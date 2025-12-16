#' SMACOF for Individual Differences
#'
#' @description Variant of `smacof::smacofIndDiff()` but with regularization
#'   (`weightlam`) in case of singularity in weight matrices (e.g., due to
#'   missingness) and optional weight scaling (`weights`) which multiples
#'   weights for each data source by some amount.
#'
#' @inheritParams smacof::smacofIndDiff
#' @param weights Vector of length equal to the number of distance matrices.
#'   If provided, weights for each distance matrix are multiplied by the
#'   corresponding value in `weights`. Default is `NULL` which will give
#'   each distance matrix the same weight (i.e., scaled by 1).
#' @param weightlam Regularization parameter for weight matrices to avoid
#'   singularity matrices. Default is 0 (i.e., no regularization). If
#'   regularization is needed, set to a small value (e.g., 1e-6).
#'
#' @export
smacofIndDiff <- function(delta, ndim = 2, type = c("ratio", "interval", "ordinal", "mspline"),
                          constraint = c("indscal", "idioscal", "identity"),
                          weights = NULL, weightmat = NULL, weightlam = 0,
                          init = "torgerson", ties = "primary",
                          verbose = FALSE, modulus = 1, itmax = 100, eps = 1e-6,
                          spline.degree = 2, spline.intKnots = 2) {
  type <- match.arg(type, c("ratio", "interval", "ordinal", "mspline"), several.ok = FALSE)
  constraint <- match.arg(constraint, c("indscal", "idioscal", "identity"), several.ok = FALSE)

  diss <- delta
  p <- ndim
  if (constraint == "indscal") constraint <- "diagonal"
  constr <- constraint
  if (!is.list(diss)) diss <- list(diss)
  if ((is.matrix(diss[[1]])) || (is.data.frame(diss[[1]]))) diss <- lapply(diss, strucprep)
  checkdiss(diss) ## sanity check

  ## --- weight matrix
  if (is.null(weightmat)) wgths <- initWeights(diss) else wgths <- weightmat
  if (!is.list(wgths)) {
    wgths <- list(wgths)
    if (length(wgths) != length(diss)) {
      wgths <- sapply(diss, function(wwr) {
        return(wgths)
      })
    }
  }
  if ((is.matrix(wgths[[1]])) || (is.data.frame(wgths[[1]]))) wgths <- lapply(wgths, strucprep)

  ## --- replace missings by 0
  diss <- lapply(diss, function(mm) {
    mm[is.na(mm)] <- 0
    return(mm)
  })
  wgths <- lapply(wgths, function(mm) {
    mm[is.na(mm)] <- 0
    return(mm)
  })

  ## --- scale weights by multiplier
  wgths_unscaled <- wgths
  if (!is.null(weights)) {
    wgths <- purrr::map2(wgths, weights, ~ .x * .y)
  }

  ## --- Prepare for optimal scaling
  trans <- type
  if (trans == "ratio") {
    trans <- "none"
  } else if (trans == "ordinal" & ties == "primary") {
    trans <- "ordinalp"
  } else if (trans == "ordinal" & ties == "secondary") {
    trans <- "ordinals"
  } else if (trans == "ordinal" & ties == "tertiary") {
    trans <- "ordinalt"
  } else if (trans == "spline") {
    trans <- "mspline"
  }
  disobj <- list()
  for (i in 1:length(diss)) {
    disobj[[i]] <- smacof::transPrep(diss[[i]], trans = trans, spline.intKnots = spline.intKnots, spline.degree = spline.degree)
  }
  ## --- end optimal scaling prep


  n <- attr(diss[[1]], "Size")
  if (p > (n - 1)) stop("Maximum number of dimensions is n-1!")

  nn <- n * (n - 1) / 2
  m <- length(diss) ## number of diss matrices
  itel <- 1

  if (is.null(attr(diss[[1]], "Labels"))) {
    for (i in 1:m) attr(diss[[i]], "Labels") <- paste(1:n)
  }

  dr <- list()
  wr <- list()
  vr <- list()
  dh <- list()

  for (j in 1:m) { # initialize weights, V, norm d as lists
    wr <- appendList(wr, vmat(wgths[[j]], weightlam))
    vr <- appendList(vr, myGenInv(wr[[j]]))
    dh <- appendList(dh, normDissN(diss[[j]], wgths_unscaled[[j]], 1))
  }
  xr <- list() # configurations as list
  sold <- sf1 <- sf2 <- 0 # stress init

  ## --- starting values
  aconf <- initConf(init, diss, n, p, inddiff = TRUE)


  bconf <- repList(diag(p), m) # 1-matrix
  for (j in 1:m) {
    xr[[j]] <- aconf %*% bconf[[j]] # same starting values for all ways
    dr[[j]] <- dist(xr[[j]]) # configuration distances
    sf1 <- sf1 + sum(wgths[[j]] * dr[[j]] * dh[[j]], na.rm = TRUE)
    sf2 <- sf2 + sum(wgths[[j]] * dr[[j]]^2, na.rm = TRUE)
  }

  lb <- sf1 / sf2 # normalization constant
  aconf <- lb * aconf
  for (j in 1:m) { # normalize X, D, compute stress
    # aconf <- lb*aconf
    xr[[j]] <- lb * xr[[j]]
    dr[[j]] <- lb * dr[[j]]
    sold <- sold + sum(wgths[[j]] * (dh[[j]] - dr[[j]])^2, na.rm = TRUE)
  }

  #--------------- begin majorization ------------------
  repeat  {
    br <- list()
    yr <- list()
    er <- list()
    sunc <- 0
    for (j in 1:m) { # compute B, Y,
      br <- appendList(br, bmat(dh[[j]], wgths[[j]], dr[[j]]))
      yr <- appendList(yr, vr[[j]] %*% (br[[j]] %*% xr[[j]]))
      er <- appendList(er, dist(yr[[j]]))
      sunc <- sunc + sum(wgths[[j]] * (dh[[j]] - er[[j]])^2, na.rm = TRUE)
    }
    scon <- sunc

    #--------- impose constraints ---------
    if (!is.null(constr)) {
      scon <- 0
      er <- list()

      #-- same configurations across ways, configuration weights I
      if (constr == "identity") {
        z <- matrix(0, n, p)
        u <- matrix(0, n, n)
        for (j in 1:m) {
          z <- z + wr[[j]] %*% yr[[j]]
          u <- u + wr[[j]]
        }
        aconf <- myGenInv(u) %*% z
        yr <- repList(aconf, m)
      }

      #-- configuration weights diagonal INDSCAL
      if (constr == "diagonal") {
        aux0 <- matrix(0, n, p)
        for (j in 1:m) {
          aux1 <- diag(crossprod(aconf, wr[[j]] %*% yr[[j]]))
          aux2 <- diag(crossprod(aconf, wr[[j]] %*% aconf))
          if ((length(aux1) == 1) && (length(aux2) == 1)) {
            bconf[[j]] <- diag(aux1 / aux2, nrow = 1, ncol = 1)
          } else {
            bconf[[j]] <- diag(aux1 / aux2)
          }
          aux0 <- aux0 + (wr[[j]] %*% yr[[j]] %*% bconf[[j]])
        }
        for (s in 1:p) {
          aux1 <- matrix(0, n, n)
          for (j in 1:m) aux1 <- aux1 + (bconf[[j]][s, s]^2) * wr[[j]]
          aconf[, s] <- myGenInv(aux1) %*% aux0[, s]
        }
        for (j in 1:m) yr[[j]] <- aconf %*% bconf[[j]]
      }

      #-- no constraints, idioscal ----
      if (constr == "idioscal") {
        aux0 <- matrix(0, n, p)
        auxk <- matrix(0, (n * p), (n * p))
        for (j in 1:m) {
          aux1 <- crossprod(aconf, wr[[j]] %*% yr[[j]])
          aux2 <- crossprod(aconf, wr[[j]] %*% aconf)
          auxb <- solve(aux2, aux1)
          bconf[[j]] <- auxb
          auxc <- crossprod(t(auxb))
          aux0 <- aux0 + (wr[[j]] %*% yr[[j]] %*% t(auxb))
          auxk <- auxk + kronecker(auxc, wr[[j]])
        }
        auxv <- kronecker(diag(p), matrix((1 / n), n, n))
        aconf <- matrix(solve(auxk + auxv, as.vector(aux0)), n, p)
        for (j in 1:m) yr[[j]] <- aconf %*% bconf[[j]]
      }

      for (j in 1:m) {
        er <- appendList(er, dist(yr[[j]]))
        scon <- scon + sum(wgths[[j]] * (dh[[j]] - er[[j]])^2, na.rm = TRUE) # constraint stress computation
      }
    }
    #-------- end constraints -----------

    snon <- scon

    snon <- 0
    dh <- list()
    for (j in 1:m) {
      do <- smacof::transform(er[[j]], disobj[[j]], w = wgths_unscaled[[j]], normq = nn)$res
      dh <- appendList(dh, do)
      snon <- snon + sum(wgths[[j]] * (dh[[j]] - er[[j]])^2)
    }


    if (verbose) {
      cat(
        "Iteration: ", formatC(itel, width = 3, format = "d"),
        " (Normalized) Stress: ", formatC(c(snon / (m * nn)), digits = 8, width = 12, format = "f"),
        "\n"
      )
    }

    if (((sold - snon) < (eps*m*nn)) || (itel == itmax)) break() # convergence

    xr <- yr
    dr <- er
    sold <- snon
    itel <- itel + 1
  }
  #------------ end majorization -----------------


  names(dh) <- names(er) <- names(yr) <- names(delta)
  cnames <- paste("D", 1:p, sep = "")
  for (i in 1:length(yr)) {
    colnames(yr[[i]]) <- cnames
    rownames(yr[[i]]) <- labels(diss[[i]])
    rownames(bconf[[i]]) <- colnames(bconf[[i]]) <- cnames
    dh[[i]] <- structure(dh[[i]], Size = n, call = quote(as.dist.default(m = b)), class = "dist", Diag = FALSE, Upper = FALSE)
    attr(dh[[i]], "Labels") <- attr(er[[i]], "Labels") <- labels(diss[[i]])
  }


  colnames(aconf) <- cnames
  rnames <- rownames(as.matrix(delta[[1]]))
  rownames(aconf) <- rnames
  names(bconf) <- names(dh)

  snon <- (snon / m) / nn # stress normalization  nn <- n*(n-1)/2, m number lists
  stress <- sqrt(snon)

  confdiss <- rep(list(NULL), m)
  for (j in 1:m) { # initialize weights, V, norm d as lists
    confdiss[[j]] <- normDissN(er[[j]], wgths_unscaled[[j]], 1)
  }


  ## stress-per-point
  # spoint <- list()
  spps <- matrix(0, m, n)
  rss <- NULL
  for (j in 1:m) {
    spointj <- spp(dh[[j]], confdiss[[j]], wgths[[j]])
    spps[j, ] <- spointj$spp ## SPP per subject
    rss[j] <- sum(spointj$resmat[lower.tri(spointj$resmat)]) ## RSS per subject
  }
  colnames(spps) <- rnames
  rownames(spps) <- names(delta)
  spp <- colMeans(spps) ## SPP
  rss <- sum(rss) ## total RSS (raw stress)

  ## stress per subject (unnormalized and normalized by weight)
  sps <- NULL
  sps_norm <- NULL
  for (j in 1:m) {
    sps[j] <- sum(wgths[[j]] * (dh[[j]] - er[[j]])^2)
    sps_norm[j] <- sps[j] / sum(wgths[[j]])
  }
  sps <- sps / sum(sps) * 100
  sps_norm <- sps_norm / sum(sps_norm) * 100
  names(sps) <- rownames(spps)
  names(sps_norm) <- rownames(spps)

  if (itel == itmax) warning("Iteration limit reached! Increase itmax argument!")

  # rename gspace and cweights dimnames
  colnames(aconf) <- paste0("CoMDS", 1:ncol(aconf))
  bconf <- purrr::map(
    bconf,
    function(.x) {
      colnames(.x) <- paste0("CoMDS", 1:ncol(.x))
      rownames(.x) <- paste0("CoMDS", 1:nrow(.x))
      return(.x)
    }
  )

  # return configurations, configuration distances, normalized observed distances
  result <- list(
    # delta = diss, dhat = dh, confdist = confdiss, conf = yr,
    gspace = aconf, cweights = bconf, stress = stress, rss = rss,
    # weightmat = wgths, resmat = spointj$resmat,
    spp = spp, spps = spps, sps = sps, sps_norm = sps_norm, ndim = p,
    model = "Consensus MDS", niter = itel, nobj = n, type = type,
    constraint = constraint, call = match.call()
  )
  class(result) <- "coMDS"
  result
}


#' SMACOF for Individual Differences with local neighborhoods
#'
#' @description Variant of `smacof::smacofIndDiff()` but for local MDS with
#'   local neighborhoods.
#'
#' @inheritParams smacof::smacofIndDiff
#' @param weights Vector of length equal to the number of distance matrices.
#'   If provided, weights for each distance matrix are multiplied by the
#'   corresponding value in `weights`. Default is `NULL` which will give
#'   each distance matrix the same weight (i.e., scaled by 1).
#' @param weightlam Regularization parameter for weight matrices to avoid
#'   singularity matrices. Default is 0 (i.e., no regularization). If
#'   regularization is needed, set to a small value (e.g., 1e-6).
#' @param tau Regularization parameter.
#' @param percentile Percentile used to define local neighborhoods. Larger
#'   percentiles result in larger neighborhoods.
#'
#' @export
localSmacofIndDiff <- function(delta, tau, percentile, ndim = 2,
                               type = c("ratio", "interval", "ordinal", "mspline"),
                               constraint = c("indscal", "idioscal", "identity"),
                               weights = NULL, weightlam = 0,
                               init = "torgerson", ties = "primary",
                               verbose = FALSE, modulus = 1, itmax = 100, eps = 1e-6,
                               spline.degree = 2, spline.intKnots = 2) {
  type <- match.arg(type, c("ratio", "interval", "ordinal", "mspline"), several.ok = FALSE)
  constraint <- match.arg(constraint, c("indscal", "idioscal", "identity"), several.ok = FALSE)

  diss <- delta
  p <- ndim
  if (constraint == "indscal") constraint <- "diagonal"
  constr <- constraint
  if (!is.list(diss)) diss <- list(diss)
  if ((is.matrix(diss[[1]])) || (is.data.frame(diss[[1]]))) diss <- lapply(diss, strucprep)
  checkdiss(diss) ## sanity check



  ## --- weight matrix
  ngbhd <- get_neighbors_list(delta, percentile)
  wgths <- init_local_weights(ngbhd)
  wgths <- lapply(wgths, strucprep)
  # tval <- get_t(delta, wgths, tau, percentile)
  tvals <- get_ts(delta, wgths, tau, percentile)
  # if (is.null(weightmat)) wgths <- initWeights(diss) else wgths <- weightmat
  # if (!is.list(wgths)) {
  # wgths <- list(wgths)
  # if (length(wgths) != length(diss)) wgths <- sapply(diss, function(wwr) return(wgths))
  # }
  # if ((is.matrix(wgths[[1]])) || (is.data.frame(wgths[[1]]))) wgths <- lapply(wgths, strucprep)

  ############################ add one matrix

  ones_matrix <- initWeights(diss[[1]])

  #############################

  ## --- replace missings by 0
  diss <- lapply(diss, function(mm) {
    mm[is.na(mm)] <- 0
    return(mm)
  })
  wgths <- lapply(wgths, function(mm) {
    mm[is.na(mm)] <- 0
    return(mm)
  })
  unif_wgths <- lapply(wgths, function(mm) {
    mm[!is.na(mm)] <- 1
    mm[is.na(mm)] <- 0
    return(mm)
  })

  ## --- scale weights by multiplier
  # wgths_unscaled <- wgths
  if (!is.null(weights)) {
    wgths <- purrr::map2(wgths, weights, ~ .x * .y)
  }

  ## --- Prepare for optimal scaling
  trans <- type
  if (trans == "ratio") {
    trans <- "none"
  } else if (trans == "ordinal" & ties == "primary") {
    trans <- "ordinalp"
  } else if (trans == "ordinal" & ties == "secondary") {
    trans <- "ordinals"
  } else if (trans == "ordinal" & ties == "tertiary") {
    trans <- "ordinalt"
  } else if (trans == "spline") {
    trans <- "mspline"
  }
  disobj <- list()
  for (i in 1:length(diss)) {
    disobj[[i]] <- smacof::transPrep(diss[[i]], trans = trans, spline.intKnots = spline.intKnots, spline.degree = spline.degree)
  }
  ## --- end optimal scaling prep


  n <- attr(diss[[1]], "Size")
  if (p > (n - 1)) stop("Maximum number of dimensions is n-1!")

  nn <- n * (n - 1) / 2
  m <- length(diss) ## number of diss matrices
  itel <- 1

  if (is.null(attr(diss[[1]], "Labels"))) {
    for (i in 1:m) attr(diss[[i]], "Labels") <- paste(1:n)
  }

  dr <- list()
  wr <- list()
  vr <- list()
  dh <- list()

  for (j in 1:m) { # initialize weights, V, norm d as lists
    wr <- appendList(wr, vmat(wgths[[j]], weightlam))
    vr <- appendList(vr, myGenInv(wr[[j]]))
    dh <- appendList(dh, normDissN(diss[[j]], wgths[[j]], 1))
    # dh <- appendList(dh, normDissN(diss[[j]], unif_wgths[[j]], 1))
  }
  xr <- list() # configurations as list
  sold <- sf1 <- sf2 <- 0 # stress init

  ## --- starting values
  aconf <- initConf(init, diss, n, p, inddiff = TRUE)


  bconf <- repList(diag(p), m) # 1-matrix
  for (j in 1:m) {
    xr[[j]] <- aconf %*% bconf[[j]] # same starting values for all ways
    dr[[j]] <- dist(xr[[j]]) # configuration distances
    sf1 <- sf1 + sum(wgths[[j]] * dr[[j]] * dh[[j]], na.rm = TRUE)
    sf2 <- sf2 + sum(wgths[[j]] * dr[[j]]^2, na.rm = TRUE)
  }

  lb <- sf1 / sf2 # normalization constant
  aconf <- lb * aconf
  for (j in 1:m) { # normalize X, D, compute stress
    # aconf <- lb*aconf
    xr[[j]] <- lb * xr[[j]]
    dr[[j]] <- lb * dr[[j]]
    sold <- sold + sum(wgths[[j]] * (dh[[j]] - dr[[j]])^2, na.rm = TRUE) - tvals[j] * sum((ones_matrix - wgths[[j]]) * (dr[[j]]), na.rm = TRUE)
  }


  #--------------- begin majorization ------------------
  repeat  {
    br <- list()
    yr <- list()
    er <- list()
    sunc <- 0
    for (j in 1:m) { # compute B, Y,
      br <- appendList(br, loco_bmat(dh[[j]], wgths[[j]], dr[[j]], tvals[j], ones_matrix))
      yr <- appendList(yr, vr[[j]] %*% (br[[j]] %*% xr[[j]]))
      er <- appendList(er, dist(yr[[j]]))
      sunc <- sunc + sum(wgths[[j]] * (dh[[j]] - er[[j]])^2, na.rm = TRUE) - tvals[j] * sum((ones_matrix - wgths[[j]]) * er[[j]], na.rm = TRUE)
    }
    scon <- sunc

    #--------- impose constraints ---------
    if (!is.null(constr)) {
      scon <- 0
      er <- list()

      #-- same configurations across ways, configuration weights I
      if (constr == "identity") {
        z <- matrix(0, n, p)
        u <- matrix(0, n, n)
        for (j in 1:m) {
          z <- z + wr[[j]] %*% yr[[j]]
          u <- u + wr[[j]]
        }
        aconf <- myGenInv(u) %*% z
        yr <- repList(aconf, m)
      }

      #-- configuration weights diagonal INDSCAL
      if (constr == "diagonal") {
        aux0 <- matrix(0, n, p)
        for (j in 1:m) {
          aux1 <- diag(crossprod(aconf, wr[[j]] %*% yr[[j]]))
          aux2 <- diag(crossprod(aconf, wr[[j]] %*% aconf))
          bconf[[j]] <- diag(aux1 / aux2)
          aux0 <- aux0 + (wr[[j]] %*% yr[[j]] %*% bconf[[j]])
        }
        for (s in 1:p) {
          aux1 <- matrix(0, n, n)
          for (j in 1:m) aux1 <- aux1 + (bconf[[j]][s, s]^2) * wr[[j]]
          aconf[, s] <- myGenInv(aux1) %*% aux0[, s]
        }
        for (j in 1:m) yr[[j]] <- aconf %*% bconf[[j]]
      }

      #-- no constraints, idioscal ----
      if (constr == "idioscal") {
        aux0 <- matrix(0, n, p)
        auxk <- matrix(0, (n * p), (n * p))
        for (j in 1:m) {
          aux1 <- crossprod(aconf, wr[[j]] %*% yr[[j]])
          aux2 <- crossprod(aconf, wr[[j]] %*% aconf)
          auxb <- solve(aux2, aux1)
          bconf[[j]] <- auxb
          auxc <- crossprod(t(auxb))
          aux0 <- aux0 + (wr[[j]] %*% yr[[j]] %*% t(auxb))
          auxk <- auxk + kronecker(auxc, wr[[j]])
        }
        auxv <- kronecker(diag(p), matrix((1 / n), n, n))
        aconf <- matrix(solve(auxk + auxv, as.vector(aux0)), n, p)
        for (j in 1:m) yr[[j]] <- aconf %*% bconf[[j]]
      }

      for (j in 1:m) {
        er <- appendList(er, dist(yr[[j]]))
        scon <- scon + sum(wgths[[j]] * (dh[[j]] - er[[j]])^2, na.rm = TRUE) - tvals[j] * sum((ones_matrix - wgths[[j]]) * er[[j]], na.rm = TRUE) # constraint stress computation
      }
    }
    #-------- end constraints -----------

    snon <- scon

    snon <- 0
    dh <- list()
    for (j in 1:m) {
      do <- smacof::transform(er[[j]], disobj[[j]], w = wgths[[j]], normq = nn)$res
      # do <- smacof::transform(er[[j]], disobj[[j]], w = unif_wgths[[j]], normq = nn)$res
      dh <- appendList(dh, do)
      snon <- snon + sum(wgths[[j]] * (dh[[j]] - er[[j]])^2) - tvals[j] * sum((ones_matrix - wgths[[j]]) * er[[j]], na.rm = TRUE)
    }


    if (verbose) {
      cat(
        "Iteration: ", formatC(itel, width = 3, format = "d"),
        " (Normalized) Stress: ", formatC(c(snon / (m * nn)), digits = 8, width = 12, format = "f"),
        "\n"
      )
    }

    if ((abs(sold - snon) < (eps*m*nn)) || (itel == itmax)) break() # convergence

    xr <- yr
    dr <- er
    sold <- snon
    itel <- itel + 1
  }
  #------------ end majorization -----------------


  names(dh) <- names(er) <- names(yr) <- names(delta)
  cnames <- paste("D", 1:p, sep = "")
  for (i in 1:length(yr)) {
    colnames(yr[[i]]) <- cnames
    rownames(yr[[i]]) <- labels(diss[[i]])
    rownames(bconf[[i]]) <- colnames(bconf[[i]]) <- cnames
    dh[[i]] <- structure(dh[[i]], Size = n, call = quote(as.dist.default(m = b)), class = "dist", Diag = FALSE, Upper = FALSE)
    attr(dh[[i]], "Labels") <- attr(er[[i]], "Labels") <- labels(diss[[i]])
  }


  colnames(aconf) <- cnames
  rnames <- rownames(as.matrix(delta[[1]]))
  rownames(aconf) <- rnames
  names(bconf) <- names(dh)

  snon <- (snon / m) / nn # stress normalization  nn <- n*(n-1)/2, m number lists
  stress <- sqrt(abs(snon))

  confdiss <- rep(list(NULL), m)
  for (j in 1:m) { # initialize weights, V, norm d as lists
    confdiss[[j]] <- normDissN(er[[j]], wgths[[j]], 1)
    # confdiss[[j]] <- normDissN(er[[j]], unif_wgths[[j]], 1)
  }

  ## stress-per-point
  # spoint <- list()
  spps <- matrix(0, m, n)
  rss <- NULL
  for (j in 1:m) {
    # if (is.null(weights)) {
    #   wgth_scale <- 1
    # } else {
    #   wgth_scale <- weights[j]
    # }
    spointj <- spp(dh[[j]], confdiss[[j]], wgths[[j]])
    # spointj <- spp(dh[[j]], confdiss[[j]], wgth_scale * unif_wgths[[j]])
    spps[j, ] <- spointj$spp ## SPP per subject
    rss[j] <- sum(spointj$resmat[lower.tri(spointj$resmat)]) ## RSS per subject
  }
  colnames(spps) <- rnames
  rownames(spps) <- names(delta)
  spp <- colMeans(spps) ## SPP
  rss <- sum(rss) ## total RSS (raw stress)

  ## stress per subject
  sps <- NULL
  sps_norm <- NULL
  for (j in 1:m) {
    # if (is.null(weights)) {
    #   wgth_scale <- 1
    # } else {
    #   wgth_scale <- weights[j]
    # }
    sps[j] <- sum(wgths[[j]] * (dh[[j]] - er[[j]])^2)
    # sps[j] <- wgth_scale * sum(unif_wgths[[j]] * (dh[[j]] - er[[j]])^2) / sum(unif_wgths[[j]])
    sps_norm[j] <- sps[j] / sum(wgths[[j]])
  }
  sps <- sps / sum(sps) * 100
  sps_norm <- sps_norm / sum(sps_norm) * 100
  names(sps) <- rownames(spps)
  names(sps_norm) <- rownames(spps)

  if (itel == itmax) warning("Iteration limit reached! Increase itmax argument!")

  # rename gspace and cweights dimnames
  colnames(aconf) <- paste0("LoCoMDS", 1:ncol(aconf))
  bconf <- purrr::map(
    bconf,
    function(.x) {
      colnames(.x) <- paste0("LoCoMDS", 1:ncol(.x))
      rownames(.x) <- paste0("LoCoMDS", 1:nrow(.x))
      return(.x)
    }
  )

  # return configurations, configuration distances, normalized observed distances
  result <- list(
    # delta = diss, dhat = dh, confdist = confdiss, conf = yr,
    gspace = aconf, cweights = bconf, stress = stress, rss = rss,
    # weightmat = wgths, resmat = spointj$resmat,
    spp = spp, spps = spps, sps = sps, sps_norm = sps_norm, ndim = p,
    model = "Local Consensus MDS", tau = tau, percentile = percentile,
    niter = itel, nobj = n, type = type,
    constraint = constraint, call = match.call()
  )
  class(result) <- "locoMDS"
  return(result)
}
