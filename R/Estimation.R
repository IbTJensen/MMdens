# Kernels
k <- function(r) ifelse(abs(r) < 1, 3/4*(1-r^2), 0)
# k_b <- function(r, b) ifelse(b == 0, Inf, ifelse(abs(r) < b, k(r / b) / b, 0))
k_b <- function(r, b) {
  if (b == 0) {
    return(rep(Inf, length(r)))
  } else {
    return(ifelse(r < b, k(r / b) / b, 0))
  }
}

# # Transformations
# expit_R <- function(x, R) R*(1/(1+exp(-x)))
# logit_R <- function(x, R) log(x/(R-x))

#' @importFrom Rcpp sourceCpp
#' @useDynLib MMdens, .registration = TRUE
NULL

#' @importFrom spatstat.geom crosspairs
#' @importFrom data.table data.table setcolorder
#' @importFrom spatstat.explore edge.Ripley
table_construct <- function(X, Z) {
  Neighbours <- crosspairs(X, Z, rmax = Inf)
  Z_dt <- data.table(v_x = Z$x, v_y = Z$y, Z_v = Z$marks)
  info_dt <- data.table(
    u_x = Neighbours$xi,
    u_y = Neighbours$yi,
    v_x = Neighbours$xj,
    v_y = Neighbours$yj,
    dist = Neighbours$d
  )
  info_dt <- merge(info_dt, Z_dt, by = c("v_x", "v_y"))
  setcolorder(info_dt, c("u_x", "u_y", "v_x", "v_y", "Z_v", "dist"))
  suppressWarnings(
    info_dt[,
      e := edge.Ripley(X = ppp(x = v_x, y = v_y, window = X$window), r = dist)
    ]
  )
  return(info_dt)
}

#' @importFrom spatstat.geom area
hatc0 <- function(info_dt, X, Z, r, b, divisor){
  N_tau <- Z$n
  lambda <- X$n/area(X$window)

  info_dt <- info_dt[dist > r-b & dist < r+b]
  info_dt[,k:=k_b(r-dist, b)]

  # Note: Markus' kode og besrkivelse bruger 2*pi*r, mens Rasmus'
  # besrkivelse bruger ||u-v||. Sidstenævnte svarer til at skifte
  # 2*pi*r ud med dist nedenfor. Førstenævnte giver dog et estimat i tråd med
  # forventningen, mens sidstenævnte ikke gør.
  # info_dt[,X_sum_terms:=(k*e)/(lambda*(2*pi*r))]
  if (divisor == "r") {
    info_dt[, X_sum_terms := Z_v * (k * e) / (lambda * (2 * pi * r))]
  }
  if (divisor == "dist") {
    info_dt[, X_sum_terms := Z_v * (k * e) / (lambda * (2 * pi * dist))]
  }
  # info_sum <- info_dt[,.(sum_terms = sum(Z_v*X_sum_terms)), list(v_x, v_y)]
  # c0 <- sum(info_sum$sum_terms)/N_tau
  c0 <- sum(info_dt$X_sum_terms)/N_tau
  return(c0)
}

#' @importFrom EstimationTools gauss_quad
Mise_est <- function(info_dt, X, Z, b, R, divisor, fast) {
  info_dt_R <- info_dt[dist < R]
  info_dt_R <- info_dt_R[order(dist)]
  N_tau <- Z$n
  lambda <- X$n / area(X$window)

  if (fast == FALSE) {
    # Extracting the indicies in info_dt_R relevant for the sum used to estimate c_0^-(u,v)(||u-v||)
    idx_dist_pairs <- Index_selection(as.matrix(info_dt_R[, c(1:4, 6)]), b)
    idx_dist_pairs <- lapply(1:length(idx_dist_pairs), function(i) {
      list(dist = info_dt_R$dist[i], idx = idx_dist_pairs[[i]])
    })

    # Caluclate the terms in the sum to estimate c0^-(u,v)(dist) for each dist
    c0_sum_terms <- lapply(
      idx_dist_pairs,
      function(x) {
        info_dt_R$Z_v[x$idx] *
          k_b(x$dist - info_dt_R$dist[x$idx], b) *
          info_dt_R$e[x$idx] /
          (lambda * 2 * pi * x$dist)
      }
    )
    rm(idx_dist_pairs)
    gc()

    # sum over all terms to esimate c0^-(u,v)(||u-v||) for every point-pair (u,v)
    c0_cv_list <- lapply(c0_sum_terms, sum)
    c0_cv <- unlist(c0_cv_list) / N_tau
    rm(c0_sum_terms)
    gc()

    info_dt_R[, c0_cv := c0_cv]
    Mise_term2 <- sum(info_dt_R$Z_v * info_dt_R$e * info_dt_R$c0_cv / lambda) /
      N_tau
  }

  # Alternative for term 2
  if (fast == TRUE) {
    c0 <- compute_c0_cpp(
      dist = info_dt_R$dist,
      Z_v = info_dt_R$Z_v,
      e = info_dt_R$e,
      lambda = lambda,
      b = b,
      N_tau = N_tau
    )

    info_dt_R[, c0 := c0]
    info_dt_R[,
      c0_cv := 1 /
        (N_tau - 1) *
        (N_tau * c0 - Z_v * k_b(0, b) * e / (2 * pi * dist * lambda))
    ]
    Mise_term2 <- sum(info_dt_R$Z_v * info_dt_R$e * info_dt_R$c0_cv / lambda) /
      (2 * pi * N_tau)
  }

  # Quadrature for term 1
  hatc <- function(info_dt, X, Z, r, b) {
    sapply(r, function(x) hatc0(info_dt, X, Z, x, b, divisor))
  }

  hatcsq <- function(r) {
    return(hatc(info_dt, X, Z, r, b)^2 * r)
  }

  Mise_term1 <- gauss_quad(
    fun = hatcsq,
    lower = 0,
    upper = R
  )

  return(Mise_term1 - 2 * Mise_term2)
}

bandwidth_selection_optim <- function(info_dt, X, Z, R, b_init = NULL, divisor, fast) {
  MISE_est_fct <- function(b) Mise_est(info_dt, X, Z, b = b, R, divisor, fast)
  # MISE_est_fct <- function(b) Mise_est(info_dt, X, Z, b = expit_R(b, R), R)

  if (!is.null(b_init)) {
    O <- optim(
      par = b_init,
      # par = ifelse(is.null(b_init), 0, b_init), # Note that logit_R(0) = R/2
      fn = MISE_est_fct,
      lower = 0,
      upper = R,
      method = "Brent"
    )
    return(O$par)
  } else {
    O <- optimise(
      f = MISE_est_fct,
      interval = c(0, R),
      maximum = FALSE
    )
    return(O$minimum)
  }
}

bandwidth_selection_grid <- function(info_dt, X, Z, R, grid, fast) {
  stopifnot(class(grid) %in% c("numeric", "double", "integer"))
  grid_search <- lapply(grid, function(b) Mise_est(info_dt, X, Z, b, R, divisor, fast))
  grid_search <- unlist(grid_search)
  best_idx <- which.min(grid_search)
  grid[best_idx]
}

#' Estimates covariance between point pattern and spatial covariance with
#' a kernel estimation procedure. Bandwidth selection is employed to selecet
#' the optimal bandwidth.
#' @param X A point pattern represented by a ppp object.
#' @param Z The spatial covariate represented as a marked ppp object where the
#' coordinates represents the spatial locations where the coviaraite is
#' measured, and the marks represents the measured values.
#' @param R The range over which the error of the esimator is assessed. This
#' is used to calculate the Mean Integrated Squared Error (MISE) used for
#' bandwidth selection. Must be a positive number.
#' @param r Vector of distance(s) for which to esimate the spatial covariance.
#' @param b_init Initial value for optimisation of the bandwidth. If left empty,
#' R/2 is used. If Grid is non-NULL, b_init should be left as NULL.
#' @param grid Grid to evaulate the MISE over. The bandwidth with the lowest MISE
#' will then be selected for the mixed moment estimator. If set to NULL the MISE
#' is found with numeric optimsation optimised. Using grid can sometimes save
#' considerable computational resources, but unless this is relevant consideration,
#' it is recommended to leave grid = NULL.
#' @param divisor Option to choose if r or ||u-v|| should be used in the divisor.
#' In the former case set divisor = "r", and in the latter case set divisor = "dist".
#' @param fast Should the leave-one-out cross-validation use a fast approximation? 
#' True by default.
#' @return Returns a list containing the estimated covariances (c0), the
#' distances at which the covariances are estimated (r), and the selected
#' bandwidth (b).
#' @export
SpatCovarEst <- function(X, Z, R, r, b_init = NULL, grid = NULL, divisor = "dist", fast = TRUE) {
  if (class(X) != "ppp" | class(Z) != "ppp") {
    stop("X and Z must be a spatstat ppp class")
  }

  if (!(class(b_init) %in% c("NULL", "numeric", "double", "integer"))) {
    stop("b_init must be either NULL or a number")
  }

  if (!is.null(b_init)) {
    if (b_init <= 0 | b_init > R) {
      stop("b_init must be between 0 and R")
    }
  }

  if (!(class(grid) %in% c("NULL", "numeric", "double", "integer"))) {
    stop("grid must be either NULL or a vector of numbers")
  }

  if (!is.null(grid)) {
    if (min(grid) <= 0 | max(grid) > R) {
      stop("All elements in grid must be between 0 and R")
    }
  }

  if(!(divisor %in% c("r", "dist"))){
    stop("Argument divisor must be 'r' or 'dist'.")
  }

  info_dt <- table_construct(X, Z)
  if (is.null(grid)) {
    b <- bandwidth_selection_optim(info_dt, X, Z, R, b_init, divisor, fast)
  } else {
    b <- bandwidth_selection_grid(info_dt, X, Z, R, grid, divisor, fast)
  }
  if (length(r) == 1) {
    c0 <- hatc0(info_dt, X, Z, r, b, divisor)
  }
  if (length(r) > 1) {
    c0 <- sapply(r, function(x) hatc0(info_dt, X, Z, x, b, divisor))
  }
  out <- list(c0 = c0, r = r, b = b)
}

#' Estimates covariance between point pattern and spatial covariance with
#' a kernel estimation procedure for a user-specified bandwidth.
#' @param X A point pattern represented by a ppp object.
#' @param Z The spatial covariate represented as a marked ppp object where the
#' coordinates represents the spatial locations where the coviaraite is
#' measured, and the marks represents the measured values.
#' @param r Vector of distance(s) for which to esimate the spatial covariance.
#' @param b The selected bandwidth.
#' @param divisor Option to choose if r or ||u-v|| should be used in the divisor.
#' In the former case set divisor = "r", and in the latter case set divisor = "dist".
#' @return Returns a list containing the estimated covariances (c0), the
#' distances at which the covariances are estimated (r), and the selected
#' bandwidth (b).
#' @export
SpatCovarEstFixed <- function(X, Z, r, b, divisor = "r") {
  if (class(X) != "ppp" | class(Z) != "ppp") {
    stop("X and Z must be a spatstat ppp class")
  }

  if(!(divisor %in% c("r", "dist"))){
    stop("Argument divisor must be 'r' or 'dist'.")
  }

  info_dt <- table_construct(X, Z)
  if (length(r) == 1) {
    c0 <- hatc0(info_dt, X, Z, r, b, divisor)
  }
  if (length(r) > 1) {
    c0 <- sapply(r, function(x) hatc0(info_dt, X, Z, x, b, divisor))
  }
  out <- list(c0 = c0, r = r, b = b)
}
