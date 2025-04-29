# app.R
library(shiny)
library(shinydashboard)
library(dplyr)
library(tidyr)
library(DT)
library(ggplot2)
library(plotly)
library(matrixcalc)
library(MASS)
library(Matrix)
library(caret)
library(lme4)
library(nlme)
library(rsconnect)

# Source helper functions
source("helper_functions.R")

# Pre-train models function
train_models <- function(train_data) {
  frequencies <- c("T250", "T500", "T1K", "T2K", "T3K", "T4K", "T6K", "T8K")
  
  # Initial estimates
  initial_estimates <- calculate_initial_estimates(train_data, frequencies)
  theta <- initial_estimates$initial_theta
  empirical_sds <- initial_estimates$empirical_sds
  
  # Calculate naive QDA covariance matrices
  naive_train_cov <- lapply(1:nrow(theta), function(k) {
    phenotype_data <- train_data[train_data$Label == k, frequencies]
    cov(as.matrix(phenotype_data))
  })
  
  # Convert train data to long format for mixed model
  train_long <- train_data %>%
    pivot_longer(cols = all_of(frequencies),
                 names_to = "Frequency",
                 values_to = "Threshold") %>%
    mutate(Frequency = factor(gsub("T", "", Frequency)),
           Label = factor(Label),
           EAR = factor(EAR),
           SID = factor(SID),
           subtype1 = as.integer(Label == 1),
           subtype2 = as.integer(Label == 2),
           subtype3 = as.integer(Label == 3),
           subtype4 = as.integer(Label == 4))
  
  # Fit mixed model
  mix_model <- lmer(Threshold ~ 1 +
                      (1 | SID) +
                      (1 | SID:EAR) +
                      (0 + subtype1:Frequency | SID) +
                      (0 + subtype2:Frequency | SID) +
                      (0 + subtype3:Frequency | SID) +
                      (0 + subtype4:Frequency | SID),
                    data = train_long
  )
  
  # Extract variance components from mixed model
  cov_components <- extract_variance_components(mix_model)
  
  # Return pre-trained models
  return(list(
    frequencies = frequencies,
    initial_estimates = initial_estimates,
    theta = theta,
    empirical_sds = empirical_sds,
    naive_train_cov = naive_train_cov,
    cov_components = cov_components,
    train_data = train_data
  ))
}

# Check if pre-trained models exist, if not, create and save them
if (!file.exists("pretrained_models.rds")) {
  message("Pre-trained models not found. Training models...")
  
  # Load MUSC data
  musc_data <- read.csv("MUSCdata.csv")
  
  # Train the models
  pretrained_models <- train_models(musc_data)
  
  # Save the trained models
  saveRDS(pretrained_models, "pretrained_models.rds")
  
  message("Models trained and saved to 'pretrained_models.rds'")
} else {
  message("Pre-trained models found. Loading from file...")
}

# UI
ui <- dashboardPage(
  dashboardHeader(title = "Ear Type Classification"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    )
  ),
  
  dashboardBody(
    tabItems(
      # Dashboard Tab
      tabItem(tabName = "dashboard",
        fluidRow(
          box(
            title = "Data Input", width = 6, status = "primary",
            tabsetPanel(
              tabPanel("Upload CSV", 
                       fileInput("file", "Upload CSV File",
                                accept = c("text/csv", "text/comma-separated-values", ".csv")),
                       checkboxInput("header", "Header", TRUE),
                       actionButton("classify_file", "Run Classification", class = "btn-primary")
              ),
              tabPanel("Manual Entry",
                       fluidRow(
                         column(6, textInput("sid", "Subject ID", "")),
                         column(6, selectInput("ear", "Ear", choices = c("Left" = "L", "Right" = "R")))
                       ),
                       fluidRow(
                         column(3, numericInput("t250", "250 Hz", 0, min = -10, max = 120)),
                         column(3, numericInput("t500", "500 Hz", 0, min = -10, max = 120)),
                         column(3, numericInput("t1k", "1000 Hz", 0, min = -10, max = 120)),
                         column(3, numericInput("t2k", "2000 Hz", 0, min = -10, max = 120))
                       ),
                       fluidRow(
                         column(3, numericInput("t3k", "3000 Hz", 0, min = -10, max = 120)),
                         column(3, numericInput("t4k", "4000 Hz", 0, min = -10, max = 120)),
                         column(3, numericInput("t6k", "6000 Hz", 0, min = -10, max = 120)),
                         column(3, numericInput("t8k", "8000 Hz", 0, min = -10, max = 120))
                       ),
                       actionButton("add_entry", "Add Entry", class = "btn-success"),
                       hr(),
                       dataTableOutput("manual_data"),
                       actionButton("classify_manual", "Run Classification", class = "btn-primary")
              ),
              tabPanel("Paste Data",
                       tags$p("Paste data in CSV format below (with or without header):"),
                       tags$p("Format: SID,EAR,T250,T500,T1K,T2K,T3K,T4K,T6K,T8K"),
                       tags$textarea(id = "csv_input", rows = 10, cols = 50, 
                                    class = "form-control",
                                    placeholder = "SID,EAR,T250,T500,T1K,T2K,T3K,T4K,T6K,T8K\n101,L,10,15,15,20,30,40,50,55\n101,R,5,10,15,25,30,35,45,50"),
                       checkboxInput("paste_header", "First row is header", TRUE),
                       actionButton("classify_paste", "Run Classification", class = "btn-primary")
              )
            )
          ),
          
          box(
            title = "Model Settings", width = 6, status = "warning",
            fluidRow(
              column(6, 
                     checkboxInput("use_naive", "Use Naive QDA", TRUE),
                     checkboxInput("use_improved", "Use Improved QDA", TRUE)
              ),
              column(6,
                     selectInput("train_data", "Training Data", 
                                choices = c("Pre-trained MUSC Model" = "musc"),
                                selected = "musc")
              )
            ),
            hr(),
            h4("Description of Ear Types:"),
            tags$ul(
              tags$li(tags$strong("Type 1:"), "Normal hearing - flat pattern with low thresholds"),
              tags$li(tags$strong("Type 2:"), "Mild to moderate hearing loss - sloping pattern"),
              tags$li(tags$strong("Type 3:"), "Moderately-severe hearing loss - sloping pattern"),
              tags$li(tags$strong("Type 4:"), "Severe to profound hearing loss - steep sloping pattern")
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Classification Results", width = 12, status = "success",
            dataTableOutput("results_table")
          )
        ),
        
        fluidRow(
          box(
            title = "Audiogram", width = 6, status = "info",
            plotlyOutput("audiogram_plot", height = "400px")
          ),
          box(
            title = "Classification Probabilities", width = 6, status = "info",
            plotlyOutput("probability_plot", height = "400px")
          )
        )
      ),
      
      # About Tab
      tabItem(tabName = "about",
        box(
          title = "About This Application", width = 12,
          tags$p("This dashboard provides real-time classification of ear types using two methods:"),
          tags$ul(
            tags$li(tags$strong("Naive QDA:"), "Quadratic Discriminant Analysis treating each ear independently"),
            tags$li(tags$strong("Improved QDA:"), "Enhanced QDA that accounts for between-ear correlation structure")
          ),
          tags$p("The classification is based on the MUSC dataset, which contains audiometric data from patients with various types of hearing loss."),
          tags$p("The models are pre-trained for faster classification and better user experience."),
          tags$p("For more information about the methods, please refer to the original research paper.")
        )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  # Load pre-trained models
  pretrained_models <- reactive({
    readRDS("pretrained_models.rds")
  })
  
  # Store user-input data
  manual_data <- reactiveVal(data.frame(
    SID = numeric(0),
    EAR = character(0),
    T250 = numeric(0),
    T500 = numeric(0),
    T1K = numeric(0),
    T2K = numeric(0),
    T3K = numeric(0),
    T4K = numeric(0),
    T6K = numeric(0),
    T8K = numeric(0)
  ))
  
  # Current data for classification
  current_data <- reactiveVal(NULL)
  
  # Classification results
  classification_results <- reactiveVal(NULL)
  
  # Add manual entry to data table
  observeEvent(input$add_entry, {
    new_row <- data.frame(
      SID = as.numeric(input$sid),
      EAR = input$ear,
      T250 = input$t250,
      T500 = input$t500,
      T1K = input$t1k,
      T2K = input$t2k,
      T3K = input$t3k,
      T4K = input$t4k,
      T6K = input$t6k,
      T8K = input$t8k
    )
    
    manual_data(rbind(manual_data(), new_row))
  })
  
  # Display manual data table
  output$manual_data <- renderDT({
    datatable(manual_data(), options = list(pageLength = 5))
  })
  
  # Parse pasted CSV data
  parse_csv_data <- function() {
    csv_text <- input$csv_input
    if (is.null(csv_text) || csv_text == "") {
      return(NULL)
    }
    
    con <- textConnection(csv_text)
    df <- read.csv(con, header = input$paste_header)
    close(con)
    
    # Ensure column names are correct if no header
    if (!input$paste_header) {
      colnames(df) <- c("SID", "EAR", "T250", "T500", "T1K", "T2K", "T3K", "T4K", "T6K", "T8K")
    }
    
    return(df)
  }
  
  # Classify function to handle both methods
  run_classification <- function(test_data) {
    if (is.null(test_data) || nrow(test_data) == 0) {
      return(NULL)
    }
    
    # Get pre-trained models
    models <- pretrained_models()
    
    # Initialize results dataframe
    results <- data.frame(
      SID = test_data$SID,
      EAR = test_data$EAR
    )
    
    # Naive QDA classification
    if (input$use_naive) {
      # Add dummy Label column to test data (required by the function)
      test_data$Label <- 1  # Temporary value, not used in prediction
      
      naive_results <- naive_QDA(
        test_data = test_data,
        theta = models$theta,
        frequencies = models$frequencies,
        train_cov_matrices = models$naive_train_cov,
        train_data = models$train_data
      )
      
      # Extract and add naive predictions
      naive_preds_flat <- unlist(naive_results$classifications)
      
      results$naive_pred <- naive_preds_flat
      
      # Extract probabilities for each phenotype
      for (i in 1:4) {
        results[[paste0("naive_prob_", i)]] <- unlist(lapply(naive_results$probabilities, function(probs) {
          if (is.matrix(probs)) {
            return(probs[, i])
          } else {
            return(probs[i])
          }
        }))
      }
    }
    
    # Improved QDA classification
    if (input$use_improved) {
      improved_results <- improved_QDA(
        test_data = test_data,
        theta = models$theta,
        empirical_sds = models$empirical_sds,
        frequencies = models$frequencies,
        cov_components = models$cov_components,
        train_data = models$train_data
      )
      
      # Extract and add improved predictions
      improved_preds_flat <- unlist(improved_results$classifications)
      
      results$improved_pred <- improved_preds_flat
      
      # Extract probabilities for each phenotype
      for (i in 1:4) {
        results[[paste0("improved_prob_", i)]] <- unlist(lapply(improved_results$probabilities, function(probs) {
          if (is.matrix(probs)) {
            return(probs[, i])
          } else {
            return(probs[i])
          }
        }))
      }
    }
    
    return(results)
  }
  
  # Handle file upload classification
  observeEvent(input$classify_file, {
    req(input$file)
    
    df <- read.csv(input$file$datapath, header = input$header)
    
    # Validate data format
    required_cols <- c("SID", "EAR", "T250", "T500", "T1K", "T2K", "T3K", "T4K", "T6K", "T8K")
    if (!all(required_cols %in% colnames(df))) {
      showNotification("Uploaded file must have columns: SID, EAR, T250, T500, T1K, T2K, T3K, T4K, T6K, T8K", 
                       type = "error")
      return()
    }
    
    current_data(df)
    classification_results(run_classification(df))
  })
  
  # Handle manual entry classification
  observeEvent(input$classify_manual, {
    req(manual_data())
    
    if (nrow(manual_data()) == 0) {
      showNotification("Please add at least one entry before classifying", type = "warning")
      return()
    }
    
    current_data(manual_data())
    classification_results(run_classification(manual_data()))
  })
  
  # Handle pasted data classification
  observeEvent(input$classify_paste, {
    df <- parse_csv_data()
    
    if (is.null(df)) {
      showNotification("Please paste valid CSV data", type = "warning")
      return()
    }
    
    # Validate data format
    required_cols <- c("SID", "EAR", "T250", "T500", "T1K", "T2K", "T3K", "T4K", "T6K", "T8K")
    if (!all(required_cols %in% colnames(df))) {
      showNotification("Pasted data must have columns: SID, EAR, T250, T500, T1K, T2K, T3K, T4K, T6K, T8K", 
                       type = "error")
      return()
    }
    
    current_data(df)
    classification_results(run_classification(df))
  })
  
  # Display classification results table
  output$results_table <- renderDT({
    req(classification_results())
    
    results <- classification_results()
    
    # Format table for display
    display_df <- data.frame(
      SID = results$SID,
      EAR = results$EAR
    )
    
    # Add naive results if available
    if (input$use_naive) {
      display_df$Naive_QDA <- paste0("Type ", results$naive_pred)
      naive_probs <- apply(results[, grep("naive_prob_", colnames(results))], 1, function(row) {
        max_prob <- max(row)
        return(paste0(round(max_prob * 100), "%"))
      })
      display_df$Naive_Confidence <- naive_probs
    }
    
    # Add improved results if available
    if (input$use_improved) {
      display_df$Improved_QDA <- paste0("Type ", results$improved_pred)
      improved_probs <- apply(results[, grep("improved_prob_", colnames(results))], 1, function(row) {
        max_prob <- max(row)
        return(paste0(round(max_prob * 100), "%"))
      })
      display_df$Improved_Confidence <- improved_probs
    }
    
    datatable(display_df, options = list(pageLength = 10))
  })
  
  # Display audiogram
  output$audiogram_plot <- renderPlotly({
    req(current_data())
    
    test_data <- current_data()
    
    # Convert to long format for plotting
    audiogram_data <- test_data %>%
      pivot_longer(
        cols = c("T250", "T500", "T1K", "T2K", "T3K", "T4K", "T6K", "T8K"),
        names_to = "Frequency",
        values_to = "Threshold"
      ) %>%
      mutate(
        Frequency_val = case_when(
          Frequency == "T250" ~ 250,
          Frequency == "T500" ~ 500,
          Frequency == "T1K" ~ 1000,
          Frequency == "T2K" ~ 2000,
          Frequency == "T3K" ~ 3000,
          Frequency == "T4K" ~ 4000,
          Frequency == "T6K" ~ 6000,
          Frequency == "T8K" ~ 8000
        ),
        Subject_Ear = paste0("Subject ", SID, " (", EAR, ")")
      )
    
    # Create audiogram with all lines on same plot but optimized legend
    p <- ggplot(audiogram_data, aes(x = Frequency_val, y = Threshold, color = Subject_Ear, group = Subject_Ear)) +
      geom_line(size = 1) +
      geom_point(size = 2) +
      scale_x_log10(
        breaks = c(250, 500, 1000, 2000, 4000, 8000),
        labels = c("250", "500", "1K", "2K", "4K", "8K")
      ) +
      scale_y_reverse(breaks = seq(-10, 120, by = 20), limits = c(120, -10)) +  # Fixed y-axis limits with wider range
      labs(
        x = "Frequency (Hz)",
        y = "Hearing Level (dB HL)",
        title = "Audiogram",
        color = "Subject-Ear"  # Better legend title
      ) +
      theme_minimal() +
      theme(
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_line(color = "gray95"),
        legend.key.size = unit(0.8, "lines"),  # Smaller legend keys
        legend.text = element_text(size = 8)   # Smaller legend text
      )
    
    # Create plotly with optimized layout
    ggplotly(p, height = 500, width = 700) %>% 
      layout(
        legend = list(
          orientation = "v",         # Vertical legend
          x = 1.02,                  # Position just outside the plot
          y = 0.5,                   # Centered vertically
          xanchor = "left",
          font = list(size = 9),     # Smaller text
          tracegroupgap = 2          # Less space between items
        ),
        margin = list(r = 100)       # Extra right margin for legend
      ) %>%
      config(displayModeBar = TRUE, scrollZoom = TRUE)  # Enable zoom with scroll wheel
  })
  
  # Display probability plot
  output$probability_plot <- renderPlotly({
    req(classification_results())
    
    results <- classification_results()
    
    # Create a long-format dataframe for plotting
    plot_data <- data.frame()
    
    # Add naive QDA probabilities if available
    if (input$use_naive) {
      naive_probs <- results %>%
        select(SID, EAR, starts_with("naive_prob_")) %>%
        pivot_longer(
          cols = starts_with("naive_prob_"),
          names_to = "Type",
          values_to = "Probability"
        ) %>%
        mutate(
          Method = "Naive QDA",
          Type = as.integer(gsub("naive_prob_", "", Type))
        )
      
      plot_data <- rbind(plot_data, naive_probs)
    }
    
    # Add improved QDA probabilities if available
    if (input$use_improved) {
      improved_probs <- results %>%
        select(SID, EAR, starts_with("improved_prob_")) %>%
        pivot_longer(
          cols = starts_with("improved_prob_"),
          names_to = "Type",
          values_to = "Probability"
        ) %>%
        mutate(
          Method = "Improved QDA",
          Type = as.integer(gsub("improved_prob_", "", Type))
        )
      
      plot_data <- rbind(plot_data, improved_probs)
    }
    
    # Create subject/ear identifier
    plot_data$Subject_Ear <- paste0("Subject ", plot_data$SID, " (", plot_data$EAR, ")")
    
    # Create probability plot
    p <- ggplot(plot_data, aes(x = factor(Type), y = Probability, fill = Method)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      facet_wrap(~ Subject_Ear) +
      scale_x_discrete(name = "Ear Type", labels = c("Type 1", "Type 2", "Type 3", "Type 4")) +
      scale_y_continuous(name = "Probability", limits = c(0, 1)) +
      labs(title = "Classification Probabilities") +
      theme_minimal()
    
    ggplotly(p)
  })
}

# Run the app
shinyApp(ui, server)