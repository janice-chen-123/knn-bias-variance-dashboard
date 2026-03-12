library(shiny)
library(ggplot2)
library(tidyr)
library(rlang)

source("function.R")

ui <- fluidPage(
  titlePanel("Monte Carlo Study of the Bias-Variance Tradeoff in k-NN"),
  tags$p("This app explores how k-NN performance changes across different simulation settings."),
  
  sidebarLayout(
    position = "right",
    sidebarPanel(
      width = 3,
      tags$h4("Simulation Settings"),
      helpText("Set the simulation parameters below."),
      
      numericInput("seed", "Random seed", value = 380, min = 1, step = 1),
      
      selectInput(
        "dimension",
        "Dimension",
        choices = c("1D" = "1d", "2D" = "2d"),
        selected = "1d"
      ),
      
      numericInput("n", "Training sample size (n)", min = 30, max = 200, value = 100, step = 10),
      sliderInput("k", "Neighbors (k)", min = 1, max = 30, value = 10, step = 1),
      sliderInput("sigma", "Noise level (sigma)", min = 0.05, max = 0.5, value = 0.2, step = 0.05),
      numericInput("B", "Monte Carlo repetitions (B)", min = 20, max = 100, value = 50, step = 10),
      sliderInput("grid_size", "Grid size", min = 21, max = 51, value = 31, step = 10),
      numericInput("kmax", "Maximum k for curve plots", min = 5, max = 30, value = 20, step = 1),
      
      tags$hr(),
      
      tags$h4("Plot Appearance"),
      helpText("Choose colors for the plots."),
      
      selectInput(
        "truth_col", "Color: true function / true MSE",
        choices = c("black", "blue", "darkgreen", "purple", "firebrick"),
        selected = "black"
      ),
      
      selectInput(
        "est_col", "Color: estimated fit / MC MSE",
        choices = c("red", "steelblue", "pink", "darkgreen", "purple"),
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
      ),
      
      tags$hr(),
      helpText("Click the button to update all plots and summary."),
      actionButton("run", "Run Simulation")
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("Fit vs Truth", plotOutput("fit_plot", height = "520px")),
        tabPanel("Bias-Variance Decomposition", plotOutput("decomp_plot", height = "520px")),
        tabPanel("Tradeoff by k", plotOutput("tradeoff_plot", height = "520px")),
        tabPanel("MC MSE vs True MSE", plotOutput("compare_plot", height = "520px"))
      ),
      
      br(),
      tags$h4("Summary"),
      helpText("A short summary of the current simulation results."),
      verbatimTextOutput("summary_text")
    )
  )
)


server <- function(input, output, session) {
  
  observe({
    
    # keep k <= n
    updateSliderInput(
      session, "k",
      max = input$n,
      value = min(input$k, input$n)
    )
    
    updateNumericInput(
      session, "kmax",
      max = input$n,
      value = min(input$kmax, input$n)
    )
    
  })
  
  # collect parameters when the simulation is run
  app_params <- eventReactive(input$run, {
    
    list(
      seed = input$seed,
      d = ifelse(input$dimension == "1d", 1, 2),
      type = input$dimension,
      n = input$n,
      k = input$k,
      sigma = input$sigma,
      B = input$B,
      grid_size = input$grid_size,
      kmax = min(input$kmax, input$n)
    )
    
  })
  
  # single simulation fit
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
  
  # bias-variance decomposition
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
  
  # bias, variance and mse curves across k
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
  
  
  