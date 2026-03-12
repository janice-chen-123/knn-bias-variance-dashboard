#' @importFrom stats rnorm runif var
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_line geom_point labs theme_minimal scale_shape_manual theme
#' @importFrom rlang .data
NULL

# Sort 2D grid values before plotting
order_surface_df <- function(df) {
  df[order(df$x2, df$x1), , drop = FALSE]
}

#' Generate covariates X from the DGP
#'
#' Draws covariates according to X ~ Unif(0,1)^d.
#'
#' @param n Integer. Training sample size.
#' @param d Integer. Dimension (1 or 2).
#' @param seed Optional integer. If provided, sets the random seed.
#'
#' @return A numeric matrix with n rows and d columns.
#' @export
dgp_x <- function(n, d = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  matrix(runif(n * d), nrow = n, ncol = d)
}

#' Deterministic target function f(x)
#'
#' Two options:
#'   "1d": f(x) = sin(2*pi*x)
#'   "2d": f(x1,x2) = sin(pi*x1) * cos(pi*x2)
#'
#' @param x Numeric vector (1d) or matrix with 2 columns (2d).
#' @param type Character. One of "1d" or "2d".
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
  sin(pi * x[, 1]) * cos(pi * x[, 2])
}

#' Generate responses Y from the DGP
#'
#' Generates Y = f(X) + eps, where eps ~ N(0, sigma^2).
#'
#' @param x Numeric matrix of covariates.
#' @param sigma Noise standard deviation.
#' @param type "1d" or "2d".
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

  if (is.null(dim(x_test))) {
    x_test <- matrix(as.numeric(x_test), ncol = ncol(x_train))
  } else {
    x_test <- as.matrix(x_test)
  }

  if (ncol(x_test) != ncol(x_train)) {
    stop("x_test and x_train must have the same number of columns.")
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
    gg <- expand.grid(x1 = axis, x2 = axis)
    return(list(x_grid = as.matrix(gg), axis = axis))
  }

  stop("Only d = 1 or 2 supported.")
}

#' Monte Carlo Bias-Variance decomposition for k-NN
#'
#' Estimate bias, variance, and MSE of k-NN predictions on a grid using Monte Carlo simulation.
#'
#' @param n Training sample size.
#' @param k Number of neighbors.
#' @param B Number of Monte Carlo repetitions.
#' @param sigma Noise standard deviation.
#' @param d Dimension (1 or 2).
#' @param type "1d" or "2d".
#' @param grid_size Number of grid points per axis.
#' @param seed Optional random seed.
#'
#' @return A list containing simulation settings, grid axis, and pointwise results.
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

# Compute average MSE over a sequence of k values
#' @noRd
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

# Plot 1D pointwise error components
#' @noRd
plot_mse_1d <- function(res,
                        show = c("bias2", "variance", "mse_y"),
                        colors = NULL) {
  if (res$settings$d != 1) stop("plot_mse_1d() requires d = 1.")

  pw <- res$pointwise
  df <- pw[, c("x1", show), drop = FALSE]
  names(df)[1] <- "x"

  df_long <- tidyr::pivot_longer(
    df,
    cols = -x,
    names_to = "component",
    values_to = "value"
  )

  p <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(x = .data$x, y = .data$value, color = .data$component)
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(
      title = sprintf(
        "k-NN pointwise error decomposition (k = %s, n = %s, B = %s)",
        res$settings$k, res$settings$n, res$settings$B
      ),
      x = "x",
      y = "Error component",
      color = "Component"
    ) +
    ggplot2::theme(legend.position = "top")

  if (!is.null(colors)) p <- p + ggplot2::scale_color_manual(values = colors)
  p
}

# Compute average bias, variance, and MSE across k
#' @noRd
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

# Plot bias-variance tradeoff curves
#' @noRd
plot_curve_components <- function(curve_df,
                                  show = c("avg_bias2", "avg_variance", "avg_mse_y"),
                                  colors = NULL) {
  needed <- c("k", show)
  if (!all(needed %in% names(curve_df))) {
    stop("curve_df must contain columns: ", paste(needed, collapse = ", "))
  }

  df <- curve_df[, needed, drop = FALSE]
  df_long <- tidyr::pivot_longer(df, cols = -k, names_to = "component", values_to = "value")

  p <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(x = .data$k, y = .data$value, color = .data$component, shape = .data$component)
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(
      title = "Bias-variance tradeoff across k",
      x = "k",
      y = "Average over grid",
      color = "Component",
      shape = "Component"
    ) +
    ggplot2::theme(legend.position = "top")

  if (!is.null(colors)) p <- p + ggplot2::scale_color_manual(values = colors)
  p
}

# Simulate one fitted k-NN model
simulate_one_fit <- function(n = 200,
                             k = 10,
                             sigma = 0.2,
                             d = 1,
                             type = c("1d", "2d"),
                             grid_size = 101,
                             seed = NULL) {
  type <- match.arg(type)
  if (!is.null(seed)) set.seed(seed)

  x_train <- dgp_x(n = n, d = d)
  y_train <- dgp_y(x_train, sigma = sigma, type = type)

  grid <- make_grid(grid_size = grid_size, d = d)
  x_grid <- grid$x_grid

  f_true <- if (d == 1) dgp_f(x_grid[, 1], type = type) else dgp_f(x_grid, type = type)
  f_hat <- knn_predict(x_train, y_train, x_grid, k = k)

  fit <- data.frame(
    x1 = x_grid[, 1],
    x2 = if (d == 2) x_grid[, 2] else NA_real_,
    f_true = f_true,
    f_hat = f_hat
  )

  train <- data.frame(
    x1 = x_train[, 1],
    x2 = if (d == 2) x_train[, 2] else NA_real_,
    y = y_train
  )

  list(
    settings = list(n = n, k = k, sigma = sigma, d = d,
                    type = type, grid_size = grid_size, seed = seed),
    grid_axis = grid$axis,
    fit = fit,
    train = train
  )
}

# Compare Monte Carlo MSE with theoretical MSE
mc_knn_compare_curve <- function(k_values,
                                 n = 200,
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
  f_true <- if (d == 1) dgp_f(x_grid[, 1], type = type) else dgp_f(x_grid, type = type)

  out <- lapply(seq_along(k_values), function(i) {
    k <- k_values[i]
    preds <- matrix(NA_real_, nrow = B, ncol = nrow(x_grid))
    mc_mse_each <- numeric(B)

    if (!is.null(seed)) set.seed(seed + i)

    for (b in seq_len(B)) {
      x_train <- dgp_x(n = n, d = d)
      y_train <- dgp_y(x_train, sigma = sigma, type = type)

      pred_b <- knn_predict(x_train, y_train, x_grid, k = k)
      preds[b, ] <- pred_b

      y_test_b <- f_true + stats::rnorm(length(f_true), mean = 0, sd = sigma)
      mc_mse_each[b] <- mean((y_test_b - pred_b)^2)
    }

    bias2_avg <- mean((colMeans(preds) - f_true)^2)
    variance_avg <- mean(apply(preds, 2, stats::var))
    theoretical_mse <- bias2_avg + variance_avg + sigma^2

    data.frame(
      k = k,
      monte_carlo_mse = mean(mc_mse_each),
      theoretical_mse = theoretical_mse
    )
  })

  do.call(rbind, out)
}

# Plot 1D fit against the true function
plot_fit_vs_truth_1d <- function(res,
                                 truth_col = "black",
                                 est_col = "red",
                                 point_col = "grey40") {
  if (res$settings$d != 1) stop("plot_fit_vs_truth_1d() requires d = 1.")

  fit_df <- res$fit
  train_df <- res$train

  ggplot2::ggplot() +
    ggplot2::geom_point(
      data = train_df,
      ggplot2::aes(x = .data$x1, y = .data$y),
      color = point_col, alpha = 0.6, size = 2
    ) +
    ggplot2::geom_line(
      data = fit_df,
      ggplot2::aes(x = .data$x1, y = .data$f_true, color = "True f(x)"),
      linewidth = 1.2
    ) +
    ggplot2::geom_line(
      data = fit_df,
      ggplot2::aes(x = .data$x1, y = .data$f_hat, color = "Estimated fhat(x)"),
      linewidth = 1.2
    ) +
    ggplot2::scale_color_manual(values = c("True f(x)" = truth_col, "Estimated fhat(x)" = est_col)) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(title = "k-NN fit versus true function", x = "x", y = "Function value", color = NULL) +
    ggplot2::theme(legend.position = "top")
}

# Plot true and estimated surfaces in 2D
plot_fit_vs_truth_2d <- function(res) {
  if (res$settings$d != 2) stop("plot_fit_vs_truth_2d() requires d = 2.")

  fit_df <- order_surface_df(res$fit)

  long_df <- tidyr::pivot_longer(
    fit_df[, c("x1", "x2", "f_true", "f_hat")],
    cols = c("f_true", "f_hat"),
    names_to = "surface",
    values_to = "value"
  )

  long_df$surface <- factor(
    long_df$surface,
    levels = c("f_true", "f_hat"),
    labels = c("True surface", "Estimated surface")
  )

  ggplot2::ggplot(long_df, ggplot2::aes(x = x1, y = x2, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::coord_fixed() +
    ggplot2::facet_wrap(~ surface) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "#F7F7F7",
      high = "#B2182B",
      midpoint = 0
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::labs(
      title = "True surface vs estimated surface",
      x = "x1",
      y = "x2",
      fill = "Value"
    )
}

# Plot one 2D error surface
plot_mse_surface_2d <- function(res, component = c("bias2", "variance", "mse_y")) {
  component <- match.arg(component)
  if (res$settings$d != 2) stop("plot_mse_surface_2d() requires d = 2.")

  pw <- order_surface_df(res$pointwise)
  pw$z <- pw[[component]]

  ggplot2::ggplot(pw, ggplot2::aes(x = x1, y = x2, fill = z)) +
    ggplot2::geom_tile() +
    ggplot2::coord_fixed() +
    ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#08306B") +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::labs(
      title = paste("2D surface of", component),
      x = "x1",
      y = "x2",
      fill = component
    )
}

# Plot Monte Carlo and theoretical MSE curves
plot_mse_compare_curve <- function(compare_df,
                                   mc_col = "steelblue",
                                   true_col = "firebrick") {
  needed <- c("k", "monte_carlo_mse", "theoretical_mse")
  if (!all(needed %in% names(compare_df))) {
    stop("compare_df must contain columns: k, monte_carlo_mse, theoretical_mse")
  }

  df_long <- tidyr::pivot_longer(
    compare_df,
    cols = c("monte_carlo_mse", "theoretical_mse"),
    names_to = "curve",
    values_to = "value"
  )

  ggplot2::ggplot(
    df_long,
    ggplot2::aes(x = .data$k, y = .data$value, color = .data$curve, shape = .data$curve)
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_color_manual(
      values = c("monte_carlo_mse" = mc_col, "theoretical_mse" = true_col),
      labels = c("monte_carlo_mse" = "Monte Carlo MSE", "theoretical_mse" = "Theoretical MSE")
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(
      title = "Monte Carlo MSE versus theoretical MSE",
      x = "k",
      y = "Average MSE",
      color = NULL,
      shape = NULL
    ) +
    ggplot2::theme(legend.position = "top")
}
