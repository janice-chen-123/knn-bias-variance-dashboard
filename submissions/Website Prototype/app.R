library(shiny)
library(ggplot2)
library(tidyr)
library(rlang)

source("function.R")

ui <- fluidPage(
  titlePanel("Monte Carlo Study of the Bias-Variance Tradeoff in k-NN"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("seed", "Random seed", value = 380, min = 1, step = 1),
      
      selectInput(
        "dimension",
        "Dimension",
        choices = c("1D" = "1d", "2D" = "2d"),
        selected = "1d"
      ),
      
      sliderInput("n", "Training sample size (n)", min = 30, max = 500, value = 200, step = 10),
      sliderInput("k", "Neighbors (k)", min = 1, max = 50, value = 10, step = 1),
      sliderInput("sigma", "Noise level (sigma)", min = 0.05, max = 1, value = 0.2, step = 0.05),
      sliderInput("B", "Monte Carlo repetitions (B)", min = 50, max = 500, value = 200, step = 50),
      sliderInput("grid_size", "Grid size", min = 21, max = 101, value = 51, step = 10),
      sliderInput("kmax", "Maximum k for curve plots", min = 5, max = 100, value = 30, step = 1),
      
      tags$hr(),
      
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
      
      actionButton("run", "Run Simulation")
    ),
    
    mainPanel(
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

