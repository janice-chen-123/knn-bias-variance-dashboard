#' @importFrom stats rnorm runif var
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_line geom_point labs theme_minimal scale_shape_manual theme
#' @importFrom rlang .data
NULL

# Sort 2D grid values before plotting
order_surface_df <- function(df) {
  df[order(df$x2, df$x1), , drop = FALSE]
}

# ---- helper for pretty labels (no underscores) + LaTeX ----
safe_TeX <- function(x) {
  if (requireNamespace("latex2exp", quietly = TRUE)) {
    latex2exp::TeX(x)
  } else {
    gsub("\\\\", "", gsub("\\$","", x))
  }
}

.component_label_map <- function() {
  c(
    bias2 = safe_TeX("$\\mathrm{Bias}^2(x)$"),
    variance = safe_TeX("$\\mathrm{Var}(\\hat{f}(x))$"),
    mse_f = safe_TeX("$\\mathrm{MSE}_f(x)$"),
    mse_y = safe_TeX("$\\mathrm{MSE}_Y(x)$"),
    avg_bias2 = safe_TeX("Average $\\mathrm{Bias}^2$"),
    avg_variance = safe_TeX("Average Variance"),
    avg_mse_f = safe_TeX("Average $\\mathrm{MSE}_f$"),
    avg_mse_y = safe_TeX("Average $\\mathrm{MSE}_Y$")
  )
}

.component_title_map <- function() {
  c(
    bias2 = "Bias^2 shown as a 3D surface",
    variance = "Variance shown as a 3D surface",
    mse_f = "Estimation MSE shown as a 3D surface",
    mse_y = "Response MSE shown as a 3D surface"
  )
}

#' Generate covariates X from the DGP
#' @export
dgp_x <- function(n, d = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  matrix(runif(n * d), nrow = n, ncol = d)
}

#' Deterministic target function f(x)
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
#' @export
dgp_y <- function(x, sigma, type = c("1d", "2d"), seed = NULL) {
  type <- match.arg(type)
  x <- as.matrix(x)
  if (!is.null(seed)) set.seed(seed)
  
  mu <- dgp_f(if (ncol(x) == 1) x[, 1] else x, type = type)
  mu + rnorm(nrow(x), mean = 0, sd = sigma)
}

#' k-NN regression prediction
#' @export
knn_predict <- function(x_train, y_train, x_test, k) {
  
  x_train <- as.matrix(x_train)
  y_train <- as.numeric(y_train)
  
  if (is.null(dim(x_test))) {
    x_test <- matrix(as.numeric(x_test), ncol = ncol(x_train))
  } else {
    x_test <- as.matrix(x_test)
  }
  
  if (ncol(x_test) != ncol(x_train)) stop("x_test and x_train must have the same number of columns.")
  
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

#' Create evaluation grid on [0,1]^d
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

#' Plot MSE decomposition (1d) using ggplot2
#' @export
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
  
  lab <- .component_label_map()
  breaks <- show
  
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
      x = safe_TeX("$x$"),
      y = safe_TeX("Error component"),
      color = NULL
    ) +
    ggplot2::theme(legend.position = "top")
  
  if (!is.null(colors)) {
    p <- p + ggplot2::scale_color_manual(
      values = colors,
      breaks = breaks,
      labels = unname(lab[breaks])
    )
  } else {
    p <- p + ggplot2::scale_color_discrete(
      breaks = breaks,
      labels = unname(lab[breaks])
    )
  }
  
  p
}

#' Compute average Bias^2 / Variance / MSE components across k values
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

#' Plot bias-variance tradeoff vs k using ggplot2
#' @export
plot_curve_components <- function(curve_df,
                                  show = c("avg_bias2", "avg_variance", "avg_mse_y"),
                                  colors = NULL) {
  
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
  
  lab <- .component_label_map()
  breaks <- show
  
  p <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(
      x = .data$k,
      y = .data$value,
      color = .data$component
    )
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(
      title = "Bias-variance tradeoff across k",
      x = safe_TeX("$k$"),
      y = "Average over grid",
      color = NULL
    ) +
    ggplot2::theme(legend.position = "top")
  
  if (!is.null(colors)) {
    p <- p + ggplot2::scale_color_manual(
      values = colors,
      breaks = breaks,
      labels = unname(lab[breaks])
    )
  } else {
    p <- p + ggplot2::scale_color_discrete(
      breaks = breaks,
      labels = unname(lab[breaks])
    )
  }
  
  p
}


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
  
  f_true <- if (d == 1) {
    dgp_f(x_grid[, 1], type = type)
  } else {
    dgp_f(x_grid, type = type)
  }
  
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
    settings = list(
      n = n, k = k, sigma = sigma, d = d,
      type = type, grid_size = grid_size, seed = seed
    ),
    grid_axis = grid$axis,
    fit = fit,
    train = train
  )
}

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
  
  f_true <- if (d == 1) {
    dgp_f(x_grid[, 1], type = type)
  } else {
    dgp_f(x_grid, type = type)
  }
  
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
      ggplot2::aes(x = .data$x1, y = .data$f_true, color = "true"),
      linewidth = 1.2
    ) +
    ggplot2::geom_line(
      data = fit_df,
      ggplot2::aes(x = .data$x1, y = .data$f_hat, color = "est"),
      linewidth = 1.2
    ) +
    ggplot2::scale_color_manual(
      values = c("true" = truth_col, "est" = est_col),
      breaks = c("true", "est"),
      labels = c(
        safe_TeX("True $f(x)$"),
        safe_TeX("Estimated $\\hat{f}(x)$")
      )
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(
      title = "k-NN fit versus true function",
      x = safe_TeX("$x$"),
      y = safe_TeX("$y$"),
      color = NULL
    ) +
    ggplot2::theme(legend.position = "top")
}

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
  
  ggplot2::ggplot(long_df, ggplot2::aes(x = .data$x1, y = .data$x2, fill = .data$value)) +
    ggplot2::geom_tile() +
    ggplot2::coord_fixed() +
    ggplot2::facet_wrap(~ surface) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::labs(
      title = "True surface vs estimated surface",
      x = safe_TeX("$x_1$"),
      y = safe_TeX("$x_2$"),
      fill = safe_TeX("Value")
    )
}

plot_mse_surface_2d <- function(res, component = c("bias2", "variance", "mse_y")) {
  component <- match.arg(component)
  
  if (res$settings$d != 2) stop("plot_mse_surface_2d() requires d = 2.")
  
  pw <- order_surface_df(res$pointwise)
  pw$z <- pw[[component]]
  
  lab <- .component_label_map()
  ttl <- .component_title_map()
  
  ggplot2::ggplot(pw, ggplot2::aes(x = .data$x1, y = .data$x2, fill = .data$z)) +
    ggplot2::geom_tile() +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::labs(
      title = unname(ttl[component]),
      x = safe_TeX("$x_1$"),
      y = safe_TeX("$x_2$"),
      fill = unname(lab[component])
    )
}

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
      values = c(
        "monte_carlo_mse" = mc_col,
        "theoretical_mse" = true_col
      ),
      labels = c(
        "monte_carlo_mse" = safe_TeX("Monte Carlo MSE"),
        "theoretical_mse" = safe_TeX("Theoretical MSE")
      )
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(
      title = "Monte Carlo MSE versus theoretical MSE",
      x = safe_TeX("$k$"),
      y = safe_TeX("Average MSE"),
      color = NULL,
      shape = NULL
    ) +
    ggplot2::theme(legend.position = "top")
}

surface_matrix <- function(df, value_col) {
  df <- order_surface_df(df)
  x_vals <- sort(unique(df$x1))
  y_vals <- sort(unique(df$x2))
  z_mat <- matrix(df[[value_col]], nrow = length(y_vals), ncol = length(x_vals), byrow = TRUE)
  list(x = x_vals, y = y_vals, z = z_mat)
}

plot_fit_vs_truth_3d <- function(res) {
  if (res$settings$d != 2) stop("plot_fit_vs_truth_3d() requires d = 2.")
  if (!requireNamespace("plotly", quietly = TRUE)) stop("Package 'plotly' is required for 3D plots.")

  surf_true <- surface_matrix(res$fit, "f_true")
  surf_hat  <- surface_matrix(res$fit, "f_hat")

  p_true <- plotly::plot_ly(
    x = surf_true$x,
    y = surf_true$y,
    z = surf_true$z,
    type = "surface",
    colorscale = "Blues",
    showscale = FALSE
  ) |> plotly::layout(
    title = list(text = "True function shown as a 3D surface"),
    scene = list(
      xaxis = list(title = "Input 1"),
      yaxis = list(title = "Input 2"),
      zaxis = list(title = "True Function value")
    )
  )

  p_hat <- plotly::plot_ly(
    x = surf_hat$x,
    y = surf_hat$y,
    z = surf_hat$z,
    type = "surface",
    colorscale = "Blues",
    showscale = TRUE
  ) |> plotly::layout(
    title = list(text = "Estimated k-NN function shown as a 3D surface"),
    scene = list(
      xaxis = list(title = "Input 1"),
      yaxis = list(title = "Input 2"),
      zaxis = list(title = "Estimated function value")
    )
  )

  plotly::subplot(p_true, p_hat, nrows = 1, titleX = TRUE, titleY = TRUE, margin = 0.05)
}

plot_mse_surface_3d <- function(res, component = c("bias2", "variance", "mse_y")) {
  component <- match.arg(component)
  if (res$settings$d != 2) stop("plot_mse_surface_3d() requires d = 2.")
  if (!requireNamespace("plotly", quietly = TRUE)) stop("Package 'plotly' is required for 3D plots.")
  
  pw <- order_surface_df(res$pointwise)
  surf <- surface_matrix(pw, component)
  ttl <- .component_title_map()
  zlab <- c(
    bias2 = "Bias^2",
    variance = "Variance",
    mse_f = "Estimation MSE",
    mse_y = "Response MSE"
  )
  
  plotly::plot_ly(
    x = surf$x,
    y = surf$y,
    z = surf$z,
    type = "surface",
    colorscale = "Viridis"
  ) |> plotly::layout(
    title = list(text = unname(ttl[component])),
    scene = list(
      xaxis = list(title = "Input 1"),
      yaxis = list(title = "Input 2"),
      zaxis = list(title = unname(zlab[component]))
    )
  )
}