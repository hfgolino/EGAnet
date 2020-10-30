library(shiny)
library(shinyalert)
library(shinyjs)
library(shinyBS)
library(EGAnet)

# Interface for EGAnet----
ui <- (
  
  navbarPage(title = "Exploratory Graph Analysis", id = "tabs",
             
             # Load Data Panel
             tabPanel(
               "Load Data",
               
               # Input
               sidebarPanel(
                 
                 # Previous objects
                 uiOutput("prev_objs_ui"),
                 
                 # Data upload
                 tags$div(fileInput("data", label = "Upload Variables",
                                    accept = c(".rds", ".csv", ".xls", ".xlsx", ".txt")), id = "data"),
                 
                 # Data Example
                 actionButton("data_example", label = "See Example", inline = TRUE),
                 
                 br(), br(),
                 
                 actionButton("load_data", label = "Load Data", inline = TRUE)
               ),
               
               # Output
               tagList(
                 
                 mainPanel(
                   
                   tableOutput("example_data_response")
                   
                 ),
                 
                 tags$footer(htmlOutput("EGAnet_cite"), align = "left", style = "position:absolute; top:0; width:30%; height:50px; color: black; margin-left: 20px; margin-top: 700px; z-index: 1000;")
                 
               )
               
             ),
             
             # Exploratory Graph Analysis Panel
             tabPanel(
               "Exploratory Graph Analysis",
               
               # Input
               sidebarPanel(
                 
                 # Number of participants (if necessary)
                 uiOutput("n"),
                 
                 # Network estimation method
                 selectInput("model", "Network Estimation Method", c("Graphical LASSO (glasso)", 
                                                                     "Triangulated Maximally Filtered Graph (TMFG)"
                 )),
                 
                 # Community detection algorithm
                 selectInput("algorithm", "Community Detection Algorithm", c("Walktrap",
                                                                             "Louvain"
                 )),
                 
                 # Check for unidimensionality
                 selectInput("unidimensional", "Check for Unidimensionality?", c("TRUE",
                                                                                "FALSE"
                 ), selected = "TRUE"),
                 
                 # Plot type
                 selectInput("plot_type", "Plotting Package", c("GGally",
                                                                "qgraph"
                 )),
                 
                 actionButton("run_ega", label = "Estimate EGA"),
                 
                 actionButton("save_ega", label = "Save EGA Results", inline = TRUE),
                 
                 actionButton("print_ega", label = "Print EGA Methods Section", inline = TRUE)
                 
               ),
               
               # Output
               tagList(
                 
                 mainPanel(
                   
                   plotOutput("ega_plot")
                   
                 ),
                 
                 tags$footer(htmlOutput("ega_cite"), align = "left", style = "position:absolute; top:0; width:30%; height:50px; color: black; margin-left: 20px; margin-top: 500px; z-index: 1000;")
                 
               )
               
             ),
             
            # Save and Reset Results Panel
            tabPanel(
              "Save and Reset All Results",
              
              # Input
              sidebarPanel(
                
                actionButton("save_master", label = "Save All Results"),
                
                br(), br(),
                
                actionButton("reset", label = "Clear Results")
              )
              
            ),
             
             # Use shinyalert
             shinyalert::useShinyalert(),
             
             # Use shinyjs
             shinyjs::useShinyjs()
             
  )
)