#' Monte Carlo k-NN Bias–Variance Toolkit
#'
#' Core utilities for a Monte Carlo study that decomposes the mean squared error (MSE)
#' of k-nearest neighbors (k-NN) regression into bias and variance under a known data
#' generating process (DGP).
#'
#' The DGP follows:
#' \deqn{Y = f(X) + \varepsilon, \quad X \sim \mathrm{Unif}(0,1)^d,\ \varepsilon \sim \mathcal{N}(0,\sigma^2).}
#'
#' @name knn_mc_toolkit
#' @docType package
NULL


#' Generate covariates X from the DGP
#'
#' Draws covariates according to \eqn{X \sim \mathrm{Unif}(0,1)^d}.
#'
#' @param n Integer. Training sample size.
#' @param d Integer. Dimension (1 or 2 in this project).
#' @param seed Optional integer. If provided, sets the random seed.
#'
#' @return A numeric matrix with \code{n} rows and \code{d} columns.
#'
#' @examples
#' x1 <- dgp_x(n = 5, d = 1, seed = 1)
#' x2 <- dgp_x(n = 5, d = 2, seed = 1)
#' dim(x1); dim(x2)
#'
#' @export
dgp_x <- function(n, d = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  matrix(stats::runif(n * d), nrow = n, ncol = d)
}


#' Deterministic target function f(x) for the project DGP
#'
#' Deterministic target function used in the project DGP.
#' 
#' Two options are provided:
#' \itemize{
#'   \item \code{"1d"}: \eqn{f(x) = \sin(2\pi x)}
#'   \item \code{"2d"}: \eqn{f(x_1,x_2) = \sin(\sqrt{x_1^2 + x_2^2})}
#' }
#'
#' @param x Numeric vector (for \code{"1d"})
#' or a numeric matrix with 2 columns (for \code{"2d"}).
#'
#' @param type Character. One of \code{"1d"} or \code{"2d"}.
#'
#' @return Numeric vector of \eqn{f(x)} values.
#'
#' @examples
#' f1 <- dgp_f(seq(0, 1, length.out = 5), type = "1d")
#' x2 <- cbind(c(0, 1), c(0, 1))
#' f2 <- dgp_f(x2, type = "2d")
#' f1; f2
#'
#' @export
dgp_f <- function(x, type = c("1d", "2d")) {
  type <- match.arg(type)

  if (type == "1d") {
    x <- as.numeric(x)
    return(sin(2 * pi * x))
  }

  # type == "2d"
  x <- as.matrix(x)
  if (ncol(x) != 2) stop("For type='2d', x must have exactly 2 columns.")
  r <- sqrt(x[, 1]^2 + x[, 2]^2)
  sin(r)
}


#' Generate responses Y from the DGP
#'
#' Generates \eqn{Y = f(X) + \varepsilon} with \eqn{\varepsilon \sim \mathcal{N}(0,\sigma^2)}.
#'
#' @param x Numeric matrix of covariates as returned by \code{dgp_x()}.
#' @param sigma Numeric. Noise standard deviation \eqn{\sigma}.
#' @param type Character. \code{"1d"} or \code{"2d"} for the chosen target function.
#' @param seed Optional integer. If provided, sets the random seed for noise generation.
#'
#' @return A numeric vector \code{y} of length \code{nrow(x)}.
#'
#' @examples
#' x <- dgp_x(10, d = 1, seed = 1)
#' y <- dgp_y(x, sigma = 0.2, type = "1d", seed = 2)
#' length(y)
#'
#' @export
dgp_y <- function(x, sigma, type = c("1d", "2d"), seed = NULL) {
  type <- match.arg(type)
  x <- as.matrix(x)
  if (!is.null(seed)) set.seed(seed)

  mu <- dgp_f(if (ncol(x) == 1) x[, 1] else x, type = type)
  mu + stats::rnorm(nrow(x), mean = 0, sd = sigma)
}


#' k-NN regression prediction at test points
#'
#' Computes k-nearest neighbors regression predictions
#' 
#' For each test point, the prediction is the average response
#' of the k closest training observations in Euclidean distance.
#' \deqn{\hat f(x_0) = \frac{1}{k}\sum_{i \in \mathcal{N}_k(x_0)} y_i,}
#' where \eqn{\mathcal{N}_k(x_0)} denotes the indices of the \eqn{k} closest
#' training points to \eqn{x_0} in Euclidean distance.
#'
#' @param x_train Numeric matrix of training covariates (n x d).
#' @param y_train Numeric vector of training responses (length n).
#' @param x_test Numeric matrix of test covariates (m x d) or numeric vector (length m) for 1d.
#' @param k Integer. Number of neighbors.
#'
#' @return Numeric vector of predictions of length m.
#'
#' @examples
#' x <- matrix(c(0, 0.5, 1), ncol = 1)
#' y <- c(0, 1, 2)
#' pred <- knn_predict(x, y, x_test = c(0.1, 0.9), k = 1)
#' pred
#'
#' @export
knn_predict <- function(x_train, y_train, x_test, k) {
  x_train <- as.matrix(x_train)
  y_train <- as.numeric(y_train)

  if (is.vector(x_test)) {
    x_test <- matrix(as.numeric(x_test), ncol = ncol(x_train))
  } else {
    x_test <- as.matrix(x_test)
  }

  n <- nrow(x_train)
  if (k > n) stop("k cannot exceed the number of training observations.")

  # Compute squared Euclidean distances efficiently for d = 1 or d = 2 (general d also works).
  # For each test point, find indices of k smallest distances.
  preds <- numeric(nrow(x_test))
  for (j in seq_len(nrow(x_test))) {
    diffs <- sweep(x_train, 2, x_test[j, ], FUN = "-")
    d2 <- rowSums(diffs^2)
    nn_idx <- order(d2)[seq_len(k)]
    preds[j] <- mean(y_train[nn_idx])
  }
  preds
}


#' Create an evaluation grid on [0, 1]^d
#'
#' Constructs a regular grid over \eqn{[0,1]^d} for pointwise evaluation.
#'
#' @param grid_size Integer. Points per axis. For \code{d = 1}, total points equal \code{grid_size}.
#'                  For \code{d = 2}, total points equal \code{grid_size^2}.
#' @param d Integer. Dimension (1 or 2).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{x_grid}: numeric matrix of grid points.
#'   \item \code{axis}: numeric vector used for each axis (useful for plotting in 2d).
#' }
#'
#' @examples
#' g1 <- make_grid(grid_size = 5, d = 1)
#' g2 <- make_grid(grid_size = 5, d = 2)
#' nrow(g1$x_grid); nrow(g2$x_grid)
#'
#' @export
make_grid <- function(grid_size = 101, d = 1) {
  axis <- seq(0, 1, length.out = grid_size)
  if (d == 1) {
    return(list(x_grid = matrix(axis, ncol = 1), axis = axis))
  }
  if (d == 2) {
    gg <- expand.grid(axis, axis)
    return(list(x_grid = as.matrix(gg), axis = axis))
  }
  stop("This project supports d=1 or d=2.")
}


#' Monte Carlo estimate of bias–variance decomposition for k-NN regression
#'
#' Runs B simulations, fits k-NN each time, and evaluates predictions on a grid.
#'
#' Definitions (at each grid point x):
#' \itemize{
#'   \item \eqn{\mathrm{Bias}(x) = \mathbb{E}[\hat f(x)] - f(x)}
#'   \item \eqn{\mathrm{Var}(x) = \mathrm{Var}(\hat f(x))}
#'   \item \eqn{\mathrm{MSE}_f(x) = \mathrm{Bias}(x)^2 + \mathrm{Var}(x)} (for estimating f)
#'   \item \eqn{\mathrm{MSE}_Y(x) = \mathrm{MSE}_f(x) + \sigma^2} (for predicting Y; adds irreducible error)
#' }
#'
#' @param n Integer. Training sample size.
#' @param k Integer. Number of neighbors.
#' @param B Integer. Number of Monte Carlo repetitions.
#' @param sigma Numeric. Noise standard deviation.
#' @param d Integer. Dimension (1 or 2).
#' @param type Character. \code{"1d"} or \code{"2d"}.
#' @param grid_size Integer. Grid resolution (see \code{make_grid()}).
#' @param seed Optional integer. If provided, sets the random seed for reproducibility.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{settings}: simulation settings.
#'   \item \code{grid_axis}: axis values used for the grid.
#'   \item \code{pointwise}: data frame of grid coordinates, truth, and Monte Carlo estimates of
#'         bias/variance/MSE components at each grid point.
#' }
#'
#' @examples
#' res <- mc_knn_decomp(n = 50, k = 5, B = 50, sigma = 0.2, d = 1, type = "1d", grid_size = 51, seed = 1)
#' names(res)
#' head(res$pointwise)
#'
#' @export
mc_knn_decomp <- function(n = 200,
                          k = 10,
                          B = 500,
                          sigma = 0.2,
                          d = 1,
                          type = c("1d", "2d"),
                          grid_size = 101,
                          seed = NULL) {
  type <- match.arg(type)
  if (!is.null(seed)) set.seed(seed)

  grid <- make_grid(grid_size = grid_size, d = d)
  x_grid <- grid$x_grid
  f_true <- dgp_f(if (d == 1) x_grid[, 1] else x_grid, type = type)

  # Store predictions: B x m
  m <- nrow(x_grid)
  preds <- matrix(NA_real_, nrow = B, ncol = m)

  for (b in seq_len(B)) {
    x_train <- dgp_x(n = n, d = d)
    y_train <- dgp_y(x_train, sigma = sigma, type = type)
    preds[b, ] <- knn_predict(x_train, y_train, x_test = x_grid, k = k)
  }

  f_hat_mean <- colMeans(preds)
  bias <- f_hat_mean - f_true
  variance <- apply(preds, 2, stats::var)

  mse_f <- bias^2 + variance
  mse_y <- mse_f + sigma^2

  pointwise <- data.frame(
    x1 = x_grid[, 1],
    x2 = if (d == 2) x_grid[, 2] else NA_real_,
    f_true = f_true,
    f_hat_mean = f_hat_mean,
    bias = bias,
    bias2 = bias^2,
    variance = variance,
    mse_f = mse_f,
    mse_y = mse_y
  )

  list(
    settings = list(n = n, k = k, B = B, sigma = sigma, d = d, type = type, grid_size = grid_size, seed = seed),
    grid_axis = grid$axis,
    pointwise = pointwise
  )
}


#' Compute an "optimal k" curve by scanning k values
#'
#' Runs \code{mc_knn_decomp()} for each k in \code{k_values} and summarizes performance
#' by averaging MSE across grid points (an empirical approximation to integrated MSE).
#'
#' @param k_values Integer vector. Candidate k values.
#' @param n,B,sigma,d,type,grid_size,seed Passed to \code{mc_knn_decomp()}.
#' @param target Character. \code{"f"} to minimize \eqn{MSE_f}, \code{"y"} to minimize \eqn{MSE_Y}.
#'
#' @return A data frame with k and the average MSE over grid points.
#'
#' @examples
#' curve <- mc_knn_curve(k_values = c(1, 3, 5), n = 50, B = 50, sigma = 0.2, d = 1, type = "1d", grid_size = 51, seed = 1)
#' curve
#'
#' @export
mc_knn_curve <- function(k_values,
                         n = 200,
                         B = 500,
                         sigma = 0.2,
                         d = 1,
                         type = c("1d", "2d"),
                         grid_size = 101,
                         seed = NULL,
                         target = c("y", "f")) {
  type <- match.arg(type)
  target <- match.arg(target)

  out <- lapply(seq_along(k_values), function(i) {
    k <- k_values[i]
    res <- mc_knn_decomp(n = n, k = k, B = B, sigma = sigma, d = d, type = type, grid_size = grid_size,
                         seed = if (is.null(seed)) NULL else seed + i)
    pw <- res$pointwise
    avg <- if (target == "y") mean(pw$mse_y) else mean(pw$mse_f)
    data.frame(k = k, avg_mse = avg)
  })

  do.call(rbind, out)
}


#' Plot MSE decomposition components (1d)
#'
#' Convenience plotter for \code{d = 1}.
#'
#' @param res Output from \code{mc_knn_decomp()} with \code{d=1}.
#' @param show Character vector indicating which curves to show. Any of:
#'   \code{"bias2"}, \code{"variance"}, \code{"mse_f"}, \code{"mse_y"}.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect (a base R plot).
#'
#' @examples
#' res <- mc_knn_decomp(n = 50, k = 5, B = 50, sigma = 0.2, d = 1, type = "1d", grid_size = 51, seed = 1)
#' plot_mse_1d(res)
#'
#' @export
plot_mse_1d <- function(res, show = c("bias2", "variance", "mse_y")) {
  pw <- res$pointwise
  x <- pw$x1

  # base plot scaffold
  y_all <- unlist(pw[show])
  plot(x, y_all[seq_along(x)], type = "n", xlab = "x", ylab = "Error component",
       main = sprintf("k-NN MSE decomposition (k=%s, n=%s, B=%s)", res$settings$k, res$settings$n, res$settings$B))

  for (nm in show) {
    lines(x, pw[[nm]])
  }
  legend("topright", legend = show, lty = 1, bty = "n")
  invisible(NULL)
}
