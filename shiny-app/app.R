library(shiny)
library(ggplot2)
library(tidyr)
library(rlang)

source("function_final.R")

ui <- fluidPage(
  tags$head(
    tags$style(HTML("\n      .recalculating, .shiny-busy, .shiny-output-error { opacity: 1 !important; }\n      .shiny-output-error { color: #c0392b; font-weight: 600; }\n      .well { background-color: #f7f7f7; }\n    "))
  ),

  titlePanel("Monte Carlo Study of the Bias-Variance Tradeoff in k-NN"),

  sidebarLayout(
    sidebarPanel(
      width = 4,
      tabsetPanel(
        id = "side_tabs",
        tabPanel(
          "Simulation",
          br(),
          numericInput("seed", "Random seed", value = 380, min = 1, step = 1),
          radioButtons(
            "dimension",
            "Dimension",
            choices = c("1D" = "1d", "2D" = "2d"),
            selected = "1d",
            inline = TRUE
          ),
          numericInput("n", "Training sample size (n)", value = 120, min = 30, max = 500, step = 10),
          sliderInput("k", "Neighbors (k)", min = 1, max = 120, value = 10, step = 1),
          sliderInput("sigma", "Noise level (sigma)", min = 0.10, max = 0.50, value = 0.20, step = 0.05),
          numericInput("B", "Monte Carlo repetitions (B)", value = 80, min = 20, max = 500, step = 10),
          sliderInput("grid_size", "Grid size", min = 21, max = 81, value = 41, step = 10),
          sliderInput("kmax", "Maximum k for curve plots", min = 5, max = 30, value = 20, step = 1),
          helpText("This version keeps a few sliders for interaction, but leaves the heavier settings as numeric inputs so the app stays fast."),
          actionButton("run", "Run Simulation")
        ),
        tabPanel(
          "Colors",
          br(),
          selectInput(
            "truth_col", "Color: true function / true MSE",
            choices = c("black", "blue", "darkgreen", "purple", "firebrick"),
            selected = "black"
          ),
          selectInput(
            "est_col", "Color: estimated fit / MC MSE",
            choices = c("red", "steelblue", "darkorange", "purple", "darkgreen"),
            selected = "red"
          ),
          selectInput(
            "bias_col", "Color: Bias^2",
            choices = c("steelblue", "blue", "purple", "darkgreen"),
            selected = "steelblue"
          ),
          selectInput(
            "var_col", "Color: Variance",
            choices = c("darkorange", "pink", "firebrick", "brown"),
            selected = "darkorange"
          ),
          selectInput(
            "mse_col", "Color: MSE",
            choices = c("firebrick", "red", "black", "purple"),
            selected = "firebrick"
          )
        )
      )
    ),

    mainPanel(
      width = 8,
      tabsetPanel(
        tabPanel("Fit vs Truth", plotOutput("fit_plot", height = "520px")),
        tabPanel("Bias-Variance Decomposition", plotOutput("decomp_plot", height = "520px")),
        tabPanel("Tradeoff by k", plotOutput("tradeoff_plot", height = "520px")),
        tabPanel("MC MSE vs True MSE", plotOutput("compare_plot", height = "520px"))
      ),
      br(),
      tags$h4("Summary"),
      verbatimTextOutput("summary_text")
    )
  )
)

server <- function(input, output, session) {

  observe({
    n_now <- max(1, as.integer(input$n))

    updateSliderInput(
      session, "k",
      min = 1,
      max = n_now,
      value = min(as.integer(input$k), n_now)
    )

    updateSliderInput(
      session, "kmax",
      min = 5,
      max = min(30, n_now),
      value = min(as.integer(input$kmax), min(30, n_now))
    )
  })

  app_params <- eventReactive(input$run, {
    list(
      seed = as.integer(input$seed),
      d = ifelse(input$dimension == "1d", 1, 2),
      type = input$dimension,
      n = as.integer(input$n),
      k = as.integer(input$k),
      sigma = as.numeric(input$sigma),
      B = as.integer(input$B),
      grid_size = as.integer(input$grid_size),
      kmax = min(as.integer(input$kmax), as.integer(input$n))
    )
  })

  fit_res <- eventReactive(input$run, {
    p <- app_params()
    simulate_one_fit(
      n = p$n,
      k = p$k,
      sigma = p$sigma,
      d = p$d,
      type = p$type,
      grid_size = p$grid_size,
      seed = p$seed
    )
  })

  decomp_res <- eventReactive(input$run, {
    p <- app_params()
    mc_knn_decomp(
      n = p$n,
      k = p$k,
      B = p$B,
      sigma = p$sigma,
      d = p$d,
      type = p$type,
      grid_size = p$grid_size,
      seed = p$seed
    )
  })

  curve_res <- eventReactive(input$run, {
    p <- app_params()
    mc_knn_curve_components(
      k_values = 1:p$kmax,
      n = p$n,
      B = p$B,
      sigma = p$sigma,
      d = p$d,
      type = p$type,
      grid_size = p$grid_size,
      seed = p$seed
    )
  })

  compare_res <- eventReactive(input$run, {
    p <- app_params()
    mc_knn_compare_curve(
      k_values = 1:p$kmax,
      n = p$n,
      B = p$B,
      sigma = p$sigma,
      d = p$d,
      type = p$type,
      grid_size = p$grid_size,
      seed = p$seed
    )
  })

  output$fit_plot <- renderPlot({
    res <- fit_res()
    req(res)

    if (res$settings$d == 1) {
      plot_fit_vs_truth_1d(
        res,
        truth_col = input$truth_col,
        est_col = input$est_col
      )
    } else {
      plot_fit_vs_truth_2d(res)
    }
  })

  output$decomp_plot <- renderPlot({
    res <- decomp_res()
    req(res)

    if (res$settings$d == 1) {
      plot_mse_1d(
        res,
        show = c("bias2", "variance", "mse_y"),
        colors = c(
          "bias2" = input$bias_col,
          "variance" = input$var_col,
          "mse_y" = input$mse_col
        )
      )
    } else {
      plot_mse_surface_2d(res, component = "mse_y")
    }
  })

  output$tradeoff_plot <- renderPlot({
    req(curve_res())
    plot_curve_components(
      curve_res(),
      show = c("avg_bias2", "avg_variance", "avg_mse_y"),
      colors = c(
        "avg_bias2" = input$bias_col,
        "avg_variance" = input$var_col,
        "avg_mse_y" = input$mse_col
      )
    )
  })

  output$compare_plot <- renderPlot({
    req(compare_res())
    plot_mse_compare_curve(
      compare_res(),
      mc_col = input$est_col,
      true_col = input$truth_col
    )
  })

  output$summary_text <- renderPrint({
    req(compare_res(), curve_res(), app_params())

    p <- app_params()
    cmp <- compare_res()
    crv <- curve_res()

    cat(
      paste0(
        "Dimension: ", ifelse(p$d == 1, "1D", "2D"), "\n",
        "Seed: ", p$seed, "\n",
        "n = ", p$n,
        " , k = ", p$k,
        " , sigma = ", p$sigma,
        " , B = ", p$B, "\n",
        "Grid size = ", p$grid_size,
        " , max k = ", p$kmax, "\n\n",
        "Best k by Monte Carlo MSE: ", cmp$best_k_mc, "\n",
        "Best k by Theoretical MSE: ", cmp$best_k_true, "\n",
        "Best k by averaged tradeoff curve: ", crv$best_k_mse_y, "\n"
      )
    )
  })
}

shinyApp(ui = ui, server = server)
