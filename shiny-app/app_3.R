library(shiny)
library(ggplot2)
library(tidyr)
library(rlang)
library(shinycssloaders)

source("function_3.R")

ui <- fluidPage(
  withMathJax(),
  tags$head(
    tags$style(HTML("\n      body { background: #f7f7f5; color: #222; font-family: Georgia, 'Times New Roman', serif; }\n      .page-title { font-size: 34px; font-weight: 700; margin: 20px 0 4px 0; color: #1f2a44; }\n      .page-subtitle { font-size: 15px; color: #555; margin-bottom: 18px; max-width: 980px; }\n      .section-box { background: #fff; border: 1px solid #ddd8cf; border-radius: 8px; padding: 16px 18px; margin-bottom: 18px; }\n      .section-title { font-size: 20px; font-weight: 700; color: #1f2a44; margin: 0 0 12px 0; }\n      .subsection-title { font-size: 15px; font-weight: 700; text-transform: uppercase; letter-spacing: 0.04em; color: #5a4a2f; margin: 10px 0 10px 0; }\n      .note-text { font-size: 13px; color: #666; margin-top: 8px; }\n      .plot-box { background: #fff; border: 1px solid #ddd8cf; border-radius: 8px; padding: 8px 10px 0 10px; margin-bottom: 18px; }\n      .results-table { width: 100%; border-collapse: collapse; font-size: 15px; }\n      .results-table th, .results-table td { border-bottom: 1px solid #e7e1d8; padding: 10px 8px; vertical-align: top; }\n      .results-table th { width: 42%; color: #1f2a44; font-weight: 700; background: #fcfbf8; }\n      .summary-pre { margin: 0; white-space: pre-wrap; word-break: break-word; color: #333; font-size: 14px; }\n      .nav-tabs > li > a { color: #1f2a44; font-weight: 700; }\n      .nav-tabs > li.active > a, .nav-tabs > li.active > a:hover, .nav-tabs > li.active > a:focus { color: #1f2a44; background: #faf8f3; border-color: #ddd8cf #ddd8cf transparent; }\n      .btn-primary, .btn-default { background: #1f2a44; color: white !important; border: none; border-radius: 4px; font-weight: 700; }\n      .btn-primary:hover, .btn-default:hover { background: #162033; }\n      .control-label { font-weight: 700; color: #1f2a44; }\n      .shiny-output-error { color: #a12622; font-weight: 700; }\n    "))
  ),

  div(class = "page-title", "Monte Carlo Study of the Bias-Variance Tradeoff in k-NN"),
  div(
    class = "page-subtitle",
    "This interface studies how ", HTML("\\(k\\)"), ", training size, noise level, and dimension affect fitted k-NN regression curves and error decomposition under repeated simulation."
  ),

  fluidRow(
    column(
      width = 4,
      div(
        class = "section-box",
        div(class = "section-title", "Simulation Controls"),
        
        tabsetPanel(
          id = "control_tabs",
          type = "tabs",
          
          tabPanel(
            "Simulation",
            
            br(),
            div(class = "subsection-title", "Core settings"),
            numericInput("seed", "Random seed", value = 380, min = 1, step = 1),
            radioButtons("dimension", "Dimension", choices = c("1D" = "1d", "2D" = "2d"), selected = "1d"),
            numericInput("n", HTML("Training sample size \\(n\\)"), value = 120, min = 30, max = 400, step = 10),
            numericInput("B", HTML("Monte Carlo repetitions \\(B\\)"), value = 80, min = 20, max = 500, step = 10),
            
            div(class = "subsection-title", "Model settings"),
            sliderInput("k", HTML("Neighbors \\(k\\)"), min = 1, max = 120, value = 10, step = 1),
            sliderInput("sigma", HTML("Noise level \\(\\sigma\\)"), min = 0.10, max = 0.50, value = 0.20, step = 0.01),
            
            div(class = "subsection-title", "Display settings"),
            sliderInput("grid_size", "Grid size", min = 21, max = 81, value = 41, step = 10),
            sliderInput("kmax", HTML("Maximum \\(k\\) for curve plots"), min = 5, max = 30, value = 20, step = 1),
            
            p(class = "note-text",
              "For 2D simulations, smaller defaults are recommended because surface plots and repeated Monte Carlo calculations may take longer to compute."
            ),
            actionButton("run", "Run Simulation", class = "btn-primary")
          ),
          
          tabPanel(
            "Colors",
            
            br(),
            selectInput("truth_col", "True function / theoretical MSE", choices = c("black", "blue", "darkgreen"), selected = "black"),
            selectInput("est_col", "Estimated fit / Monte Carlo MSE", choices = c("red", "steelblue", "purple"), selected = "red"),
            selectInput("bias_col", "Bias component", choices = c("black", "darkgreen", "purple"), selected = "black"),
            selectInput("var_col", "Variance component", choices = c("darkorange", "pink", "firebrick"), selected = "darkorange"),
            selectInput("mse_col", "Overall error", choices = c("firebrick", "steelblue", "purple"), selected = "firebrick")
          )
        )
      )

    ),
    column(
      width = 8,
      div(
        class = "plot-box",
        tabsetPanel(
          tabPanel("Model Fit", withSpinner(plotOutput("fit_plot", height = "520px"), type = 6)),
          tabPanel("Bias-Variance Breakdown", withSpinner(plotOutput("decomp_plot", height = "520px"), type = 6)),
          tabPanel("k Selection", withSpinner(plotOutput("tradeoff_plot", height = "520px"), type = 6)),
          tabPanel("MSE Comparison", withSpinner(plotOutput("compare_plot", height = "520px"), type = 6)),
          tabPanel(
            "Model Notes",
            div(
              style = "padding: 14px 10px 8px 10px; max-width: 900px;",
              tags$h4("Model description"),
              tags$p("The app studies repeated-sample behavior of k-nearest neighbors regression."),
              tags$p("In one dimension, the target regression function is ", HTML("\\(f(x) = \\sin(2\\pi x)\\).")),
              tags$p("In two dimensions, the target surface is ", HTML("\\(f(x_1, x_2) = \\sin(\\pi x_1) \\cos(\\pi x_2)\\).")),
              tags$ul(
                tags$li(HTML("\\(n\\): number of training observations")),
                tags$li(HTML("\\(k\\): number of neighbors used in k-NN prediction")),
                tags$li(HTML("\\(\\sigma\\): standard deviation of the additive noise term")),
                tags$li(HTML("\\(B\\): number of Monte Carlo repetitions used for bias, variance, and MSE estimation")),
                tags$li("Grid size: number of evaluation locations per axis used for plotting and averaging")
              )
            )
          )
        )
      ),
      div(
        class = "section-box",
        div(class = "section-title", "Results Overview"),
        tableOutput("results_table")
      ),
      div(
        class = "section-box",
        div(class = "section-title", "Run Summary"),
        withSpinner(verbatimTextOutput("summary_text"), type = 4)
      )
    )
  )
)

server <- function(input, output, session) {
  observe({
    n_now <- max(1, as.integer(input$n))
    updateSliderInput(session, "k", min = 1, max = n_now, value = min(as.integer(input$k), n_now))
    updateSliderInput(session, "kmax", min = 5, max = min(30, n_now), value = min(as.integer(input$kmax), min(30, n_now)))
  })

  observeEvent(input$n, {
    n_now <- max(1, as.integer(input$n))
    
    new_k_max <- n_now
    new_k_val <- min(as.integer(input$k), new_k_max)
    
    new_kmax_max <- min(30, n_now)
    new_kmax_val <- min(as.integer(input$kmax), new_kmax_max)
    
    freezeReactiveValue(input, "k")
    freezeReactiveValue(input, "kmax")
    
    updateSliderInput(
      session, "k",
      min = 1,
      max = new_k_max,
      value = new_k_val
    )
    
    updateSliderInput(
      session, "kmax",
      min = 5,
      max = new_kmax_max,
      value = new_kmax_val
    )
  }, ignoreInit = FALSE)
  
  observeEvent(input$dimension, {
    if (input$dimension == "2d") {
      updateNumericInput(session, "n", value = 80)
      updateNumericInput(session, "B", value = 20)
      updateSliderInput(session, "grid_size", value = 21)
      updateSliderInput(session, "kmax", value = 10)
    }
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
    withProgress(message = "Computing fitted surface...", value = 0, {
      incProgress(0.1)
      res <- simulate_one_fit(n = p$n, k = p$k, sigma = p$sigma, d = p$d, type = p$type, grid_size = p$grid_size, seed = p$seed)
      incProgress(0.9)
      res
    })
  })

  decomp_res <- eventReactive(input$run, {
    p <- app_params()
    withProgress(message = "Computing decomposition...", value = 0, {
      incProgress(0.1)
      res <- mc_knn_decomp(n = p$n, k = p$k, B = p$B, sigma = p$sigma, d = p$d, type = p$type, grid_size = p$grid_size, seed = p$seed)
      incProgress(0.9)
      res
    })
  })

  curve_res <- eventReactive(input$run, {
    p <- app_params()
    withProgress(message = "Computing tradeoff curve across k...", value = 0, {
      incProgress(0.1)
      res <- mc_knn_curve_components(k_values = 1:p$kmax, n = p$n, B = p$B, sigma = p$sigma, d = p$d, type = p$type, grid_size = p$grid_size, seed = p$seed)
      incProgress(0.9)
      res
    })
  })

  compare_res <- eventReactive(input$run, {
    p <- app_params()
    withProgress(message = "Comparing Monte Carlo MSE with theoretical MSE...", value = 0, {
      incProgress(0.1)
      res <- mc_knn_compare_curve(k_values = 1:p$kmax, n = p$n, B = p$B, sigma = p$sigma, d = p$d, type = p$type, grid_size = p$grid_size, seed = p$seed)
      incProgress(0.9)
      res
    })
  })

  output$fit_plot <- renderPlot({
    res <- fit_res()
    req(res)
    if (res$settings$d == 1) {
      plot_fit_vs_truth_1d(res, truth_col = input$truth_col, est_col = input$est_col)
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
        colors = c("bias2" = input$bias_col, "variance" = input$var_col, "mse_y" = input$mse_col)
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
      colors = c("avg_bias2" = input$bias_col, "avg_variance" = input$var_col, "avg_mse_y" = input$mse_col)
    )
  })

  output$compare_plot <- renderPlot({
    req(compare_res())
    plot_mse_compare_curve(compare_res(), mc_col = input$est_col, true_col = input$truth_col)
  })

  safe_best_k <- function(df, metric_col) {
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NA_integer_)
    if (!metric_col %in% names(df)) return(NA_integer_)
    v <- df[[metric_col]]
    if (length(v) == 0 || all(is.na(v))) return(NA_integer_)
    idx <- which.min(replace(v, is.na(v), Inf))
    if (length(idx) == 0 || is.na(idx)) return(NA_integer_)
    df$k[idx]
  }

  fmt_num <- function(x, digits = 4) {
    if (length(x) == 0 || is.na(x) || !is.finite(x)) return("Not available")
    format(round(x, digits), nsmall = digits)
  }

  best_values <- reactive({
    req(compare_res(), curve_res(), app_params())
    cmp <- compare_res()
    crv <- curve_res()
    list(
      best_mc = safe_best_k(cmp, "monte_carlo_mse"),
      best_true = safe_best_k(cmp, "theoretical_mse"),
      best_tradeoff = safe_best_k(crv, "avg_mse_y"),
      min_mc = if (!is.data.frame(cmp) || !"monte_carlo_mse" %in% names(cmp)) NA_real_ else suppressWarnings(min(cmp$monte_carlo_mse, na.rm = TRUE)),
      min_true = if (!is.data.frame(cmp) || !"theoretical_mse" %in% names(cmp)) NA_real_ else suppressWarnings(min(cmp$theoretical_mse, na.rm = TRUE))
    )
  })

  output$results_table <- renderTable({
    req(app_params(), best_values())
    p <- app_params()
    vals <- best_values()
    data.frame(
      Metric = c(
        "Current mode",
        "Best k by Monte Carlo MSE",
        "Best k by theoretical MSE",
        "Best k by averaged tradeoff curve",
        "Minimum Monte Carlo MSE",
        "Minimum theoretical MSE"
      ),
      Value = c(
        ifelse(p$d == 1, "1D", "2D"),
        ifelse(is.na(vals$best_mc), "Not available", as.character(vals$best_mc)),
        ifelse(is.na(vals$best_true), "Not available", as.character(vals$best_true)),
        ifelse(is.na(vals$best_tradeoff), "Not available", as.character(vals$best_tradeoff)),
        fmt_num(vals$min_mc),
        fmt_num(vals$min_true)
      ),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }, striped = FALSE, bordered = FALSE, hover = FALSE, spacing = "m", width = "100%", sanitize.text.function = identity)

  output$summary_text <- renderText({
    req(compare_res(), curve_res(), app_params())
    p <- app_params()
    vals <- best_values()
    paste0(
      "Dimension: ", ifelse(p$d == 1, "1D", "2D"), "\n",
      "Seed: ", p$seed, "\n",
      "n = ", p$n,
      " , k = ", p$k,
      " , sigma = ", p$sigma,
      " , B = ", p$B, "\n",
      "Grid size = ", p$grid_size,
      " , max k = ", p$kmax, "\n\n",
      "Best k by Monte Carlo MSE: ", ifelse(is.na(vals$best_mc), "Not available", vals$best_mc), "\n",
      "Best k by Theoretical MSE: ", ifelse(is.na(vals$best_true), "Not available", vals$best_true), "\n",
      "Best k by averaged tradeoff curve: ", ifelse(is.na(vals$best_tradeoff), "Not available", vals$best_tradeoff), "\n",
      "Minimum Monte Carlo MSE: ", fmt_num(vals$min_mc), "\n",
      "Minimum Theoretical MSE: ", fmt_num(vals$min_true)
    )
  })
}

shinyApp(ui = ui, server = server)
