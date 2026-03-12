## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all()
} else {
  stop("Please install 'devtools' first: install.packages('devtools')")
}

## -----------------------------------------------------------------------------
x <- dgp_x(n = 10, d = 1, seed = 1)
y <- dgp_y(x, sigma = 0.2, type = "1d", seed = 2)
head(cbind(x, y))

## -----------------------------------------------------------------------------
res_1d <- mc_knn_decomp(
  n = 200, k = 15, B = 200,
  sigma = 0.2, d = 1, type = "1d",
  grid_size = 101, seed = 123
)

str(res_1d$settings)
head(res_1d$pointwise)

## ----fig.width=8, fig.height=5------------------------------------------------
plot_mse_1d(res_1d, show = c("bias2", "variance", "mse_y"))

## ----fig.width=8, fig.height=5, message=FALSE, warning=FALSE------------------
library(ggplot2)

pw <- res_1d$pointwise
df <- data.frame(x = pw$x1, f_true = pw$f_true, f_hat_mean = pw$f_hat_mean)

ggplot(df, aes(x = x)) +
  geom_line(aes(y = f_true, color = "True f(x)"), linewidth = 1.2) +
  geom_line(aes(y = f_hat_mean, color = "Monte Carlo mean E[f_hat(x)]"), linewidth = 1.2) +
  scale_color_manual(values = c(
    "True f(x)" = "#0000FF",
    "Monte Carlo mean E[f_hat(x)]" = "#FF00FF"
  )) +
  labs(
    title = "True regression function vs Monte Carlo mean of k-NN predictions",
    x = "x", y = "Value", color = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")

## -----------------------------------------------------------------------------
k_values <- c(1, 3, 5, 10, 15, 25, 40)
curve <- mc_knn_curve(
  k_values = k_values,
  n = 200, B = 200, sigma = 0.2,
  d = 1, type = "1d",
  grid_size = 101, seed = 10,
  target = "y"
)
curve

## ----fig.width=8, fig.height=5------------------------------------------------
curve_comp <- mc_knn_curve_components(
  k_values = k_values,
  n = 200, B = 200, sigma = 0.2,
  d = 1, type = "1d",
  grid_size = 101, seed = 10
)

plot_curve_components(curve_comp, show = c("avg_bias2", "avg_variance", "avg_mse_y"))

## ----include=FALSE------------------------------------------------------------
res_2d <- mc_knn_decomp(
  n = 200, k = 20, B = 100,
  sigma = 0.2, d = 2, type = "2d",
  grid_size = 31, seed = 99
)
head(res_2d$pointwise)

