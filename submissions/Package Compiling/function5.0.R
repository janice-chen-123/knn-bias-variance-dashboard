#' Generate covariates X from the DGP
#'
#' Draws covariates according to \eqn{X \sim \mathrm{Unif}(0,1)^d}.
#'
#' @param n Integer. Training sample size.
#' @param d Integer. Dimension (1 or 2).
#' @param seed Optional integer. If provided, sets the random seed.
#'
#' @return A numeric matrix with \code{n} rows and \code{d} columns.
#' @export
dgp_x <- function(n, d = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  matrix(runif(n * d), nrow = n, ncol = d)
}


#' Deterministic target function f(x)
#'
#' Two options:
#' \itemize{
#'   \item \code{"1d"}: \eqn{f(x) = \sin(2\pi x)}
#'   \item \code{"2d"}: \eqn{f(x_1,x_2) = \sin(\sqrt{x_1^2 + x_2^2})}
#' }
#'
#' @param x Numeric vector (1d) or matrix with 2 columns (2d).
#' @param type Character. One of \code{"1d"} or \code{"2d"}.
#'
#' @return Numeric vector.
#' @export
dgp_f <- function(x, type = c("1d", "2d")) {
  type <- match.arg(type)
  
  if (type == "1d") {
    x <- as.numeric(x)
    return(sin(2 * pi * x))
  }
  
  x <- as.matrix(x)
  if (ncol(x) != 2) stop("For type='2d', x must have exactly 2 columns.")
  r <- sqrt(x[, 1]^2 + x[, 2]^2)
  sin(r)
}


#' Generate responses Y from the DGP
#'
#' Generates \eqn{Y = f(X) + \varepsilon},
#' where \eqn{\varepsilon \sim \mathcal{N}(0,\sigma^2)}.
#'
#' @param x Numeric matrix of covariates.
#' @param sigma Noise standard deviation.
#' @param type \code{"1d"} or \code{"2d"}.
#' @param seed Optional integer.
#'
#' @return Numeric vector.
#' @export
dgp_y <- function(x, sigma, type = c("1d", "2d"), seed = NULL) {
  type <- match.arg(type)
  x <- as.matrix(x)
  if (!is.null(seed)) set.seed(seed)
  
  mu <- dgp_f(if (ncol(x) == 1) x[, 1] else x, type = type)
  mu + rnorm(nrow(x), mean = 0, sd = sigma)
}


#' k-NN regression prediction
#'
#' Computes k-nearest neighbors regression predictions.
#'
#' @param x_train Numeric matrix (n x d).
#' @param y_train Numeric vector (length n).
#' @param x_test Numeric matrix (m x d) or vector (1d).
#' @param k Integer.
#'
#' @return Numeric vector of length m.
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
  if (k > n) stop("k cannot exceed number of training observations.")
  
  preds <- numeric(nrow(x_test))
  
  for (j in seq_len(nrow(x_test))) {
    diffs <- sweep(x_train, 2, x_test[j, ], "-")
    d2 <- rowSums(diffs^2)
    nn_idx <- order(d2)[seq_len(k)]
    preds[j] <- mean(y_train[nn_idx])
  }
  
  preds
}


#' Create evaluation grid on \eqn{[0,1]^d}
#'
#' @param grid_size Integer.
#' @param d Dimension (1 or 2).
#'
#' @return List with x_grid and axis.
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
  
  stop("Only d = 1 or 2 supported.")
}


#' Monte Carlo Bias–Variance decomposition for k-NN
#'
#' @param n Training sample size.
#' @param k Number of neighbors.
#' @param B Monte Carlo repetitions.
#' @param sigma Noise sd.
#' @param d Dimension.
#' @param type Target function.
#' @param grid_size Grid resolution.
#' @param seed Optional integer.
#'
#' @return List with settings, grid_axis, and pointwise data frame.
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
  
  m <- nrow(x_grid)
  preds <- matrix(NA_real_, nrow = B, ncol = m)
  
  for (b in seq_len(B)) {
    x_train <- dgp_x(n = n, d = d)
    y_train <- dgp_y(x_train, sigma = sigma, type = type)
    preds[b, ] <- knn_predict(x_train, y_train, x_grid, k)
  }
  
  f_hat_mean <- colMeans(preds)
  bias <- f_hat_mean - f_true
  variance <- apply(preds, 2, var)
  
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
    settings = list(n = n, k = k, B = B, sigma = sigma,
                    d = d, type = type, grid_size = grid_size, seed = seed),
    grid_axis = grid$axis,
    pointwise = pointwise
  )
}


#' Compute average MSE curve across k values
#'
#' @param k_values Integer vector of candidate k values.
#' @param n Training sample size.
#' @param B Monte Carlo repetitions.
#' @param sigma Noise standard deviation.
#' @param d Dimension (1 or 2).
#' @param type Target function ("1d" or "2d").
#' @param grid_size Grid resolution.
#' @param seed Optional integer for reproducibility.
#' @param target Character: "y" to minimize MSE_Y,
#'               "f" to minimize MSE_f.
#'
#' @return Data frame with columns k and avg_mse.
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
    
    res <- mc_knn_decomp(
      n = n, k = k, B = B, sigma = sigma,
      d = d, type = type, grid_size = grid_size,
      seed = if (is.null(seed)) NULL else seed + i
    )
    
    pw <- res$pointwise
    avg <- if (target == "y") mean(pw$mse_y) else mean(pw$mse_f)
    
    data.frame(k = k, avg_mse = avg)
  })
  
  do.call(rbind, out)
}


#' Plot MSE decomposition (1d) using ggplot2
#'
#' @param res Output of mc_knn_decomp() with d = 1.
#' @param show Character vector of components.
#'
#' @return ggplot object.
#' @export
plot_mse_1d <- function(res,
                        show = c("bias2", "variance", "mse_y")) {
  
  if (res$settings$d != 1) {
    stop("plot_mse_1d() requires d = 1.")
  }
  
  pw <- res$pointwise
  df <- pw[, c("x1", show)]
  names(df)[1] <- "x"
  
  df_long <- tidyr::pivot_longer(
    df,
    cols = -x,
    names_to = "component",
    values_to = "value"
  )
  
  ggplot2::ggplot(df_long,
                  ggplot2::aes(x = x,
                               y = value,
                               color = component)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(
      title = sprintf(
        "k-NN MSE decomposition (k=%s, n=%s, B=%s)",
        res$settings$k,
        res$settings$n,
        res$settings$B
      ),
      x = "x",
      y = "Error component",
      color = "Component"
    ) +
    ggplot2::theme(legend.position = "top")
}

#' Compute average Bias^2 / Variance / MSE components across k values
#'
#' This is a convenience wrapper around \code{mc_knn_decomp()} that returns
#' the **grid-averaged** components for each candidate \code{k}. It's useful
#' for making the classic bias–variance tradeoff plot versus \code{k}.
#'
#' @param k_values Integer vector of candidate k values.
#' @param n Training sample size.
#' @param B Monte Carlo repetitions.
#' @param sigma Noise standard deviation.
#' @param d Dimension (1 or 2).
#' @param type Target function ("1d" or "2d").
#' @param grid_size Grid resolution.
#' @param seed Optional integer for reproducibility.
#'
#' @return Data frame with one row per k and columns:
#'   \code{k}, \code{avg_bias2}, \code{avg_variance}, \code{avg_mse_f}, \code{avg_mse_y}.
#' @export
mc_knn_curve_components <- function(k_values,
                                    n = 200,
                                    B = 500,
                                    sigma = 0.2,
                                    d = 1,
                                    type = c("1d", "2d"),
                                    grid_size = 101,
                                    seed = NULL) {
  
  type <- match.arg(type)
  
  out <- lapply(seq_along(k_values), function(i) {
    
    k <- k_values[i]
    
    res <- mc_knn_decomp(
      n = n, k = k, B = B, sigma = sigma,
      d = d, type = type, grid_size = grid_size,
      seed = if (is.null(seed)) NULL else seed + i
    )
    
    pw <- res$pointwise
    
    data.frame(
      k = k,
      avg_bias2 = mean(pw$bias2),
      avg_variance = mean(pw$variance),
      avg_mse_f = mean(pw$mse_f),
      avg_mse_y = mean(pw$mse_y)
    )
  })
  
  do.call(rbind, out)
}


#' Plot bias–variance tradeoff vs k using ggplot2
#'
#' @param curve_df Output of \code{mc_knn_curve_components()}.
#' @param show Character vector of components to plot.
#'
#' @return ggplot object.
#' @export
plot_curve_components <- function(curve_df,
                                  show = c("avg_bias2", "avg_variance", "avg_mse_y")) {
  
  needed <- c("k", show)
  if (!all(needed %in% names(curve_df))) {
    stop("curve_df must contain columns: ", paste(needed, collapse = ", "))
  }
  
  df <- curve_df[, needed, drop = FALSE]
  
  df_long <- tidyr::pivot_longer(
    df,
    cols = -k,
    names_to = "component",
    values_to = "value"
  )
  
  ggplot2::ggplot(df_long, ggplot2::aes(x = k, y = value,
                                        color = component, shape = component)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_shape_manual(values = c(
      "avg_bias2" = 16,
      "avg_variance" = 15,
      "avg_mse_y" = 17
    )) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(
      title = "Bias–variance tradeoff across k (averaged over x)",
      x = "k",
      y = "Average over x",
      color = "Component",
      shape = "Component"
    ) +
    ggplot2::theme(legend.position = "top")
}
