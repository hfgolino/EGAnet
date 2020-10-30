# Code for EGAnet----
server <- function(input, output, session)
{
  # Keep previous environment
  prev.env <<- ls(envir = globalenv())
  
  # Check if anything exists in previous environment
  if(length(prev.env) != 0){
    
    ## Get environment objects
    env.objs <<- names(prev.env)
    
    ## Render drop down for objects
    output$prev_objs_ui <- renderUI({
      selectInput("prev_obj", label = "Objects Detected in R's environment",
                  choices = c("", env.objs), selected = 1)
    })
  
  }
  
  # semna citation
  output$EGAnet_cite <- renderUI({
    
    year <<- unlist(strsplit(as.character(Sys.Date()), split = "-"))[1]
    
    HTML(
      
      paste('<b>Please cite:</b><br>
            Golino, H., & Christensen, A. P. (',
            year, 
            '). EGAnet: Exploratory Graph Analysis -- A framework for estimating the number of dimensions in multivariate data using network psychometrics. Retrived from <a href="https://cran.r-project.org/package=EGAnet">https://cran.r-project.org/package=EGAnet </a>',
            sep = "")
    )
    
  })
  
  ###################
  #### HIDE TABS ####
  ###################
  
  hideTab(inputId = "tabs", target = "Exploratory Graph Analysis")
  hideTab(inputId = "tabs", target = "Save and Reset All Results")
  
  ###########################
  #### HIDE SAVE BUTTONS ####
  ###########################
  
  shinyjs::hide("save_ega")
  shinyjs::hide("reset")
  shinyjs::hide("print_ega")
  
  #######################
  #### DATA EXAMPLES ####
  #######################
  
  # Data
  observeEvent(input$data_example,
               {
                 output$example_data_response <- renderTable({head(EGAnet::wmt2[,7:24])},
                                                             rownames = TRUE,
                                                             caption = "Participants by Variables",
                                                             caption.placement = getOption("xtable.caption.placement", "top")
                 )
               }
  )
  
  #############
  #### EGA ####
  #############
  
  # Load Data panel
  observeEvent(input$load_data,
               {
                 # Let user know
                 showNotification("Loading data...")
                 
                 if(!is.null(input$prev_obj)){
                   
                   # Load data from R environment
                   if(input$prev_obj != ""){
                     dat <<- get(input$prev_obj, envir = globalenv())
                   }
                   
                 }
                 
                 # Load preprocessed data
                 if(!is.null(input$data)){
                   dat <<- EGAnet:::read.data(input$data$datapath)
                 }
                 
                 # Load data from EGAnet package
                 if(!exists("dat")){
                   dat <<- EGAnet::wmt2[,7:24]
                 }
                 
                 # Check for number of participants
                 if(ncol(dat) == nrow(dat)){
                   
                   output$n <- renderUI({
                     numericInput("n", "Number of Participants", value = 0)
                   })
                   
                 }else{cases <<- nrow(dat)}
                 
                 # Show network estimation tab
                 showTab(inputId = "tabs", target = "Exploratory Graph Analysis")
                 
                 # Show save and reset tab
                 showTab(inputId = "tabs", target = "Save and Reset All Results")
                 
                 # Show save data button
                 shinyjs::show("save_data")
                 
                 # Print waiting message
                 # FOR R PACKAGE AND WEB
                 shinyalert::shinyalert(title = "Data Loaded Successfully",
                                        type = "info",
                                        showConfirmButton = TRUE)
                 
                 # Move on to network estimation tab
                 updateTabsetPanel(session, "tabs",
                                   selected = "Exploratory Graph Analysis")
                 
               }
  )
  
  # Network Estimation panel
  observeEvent(input$run_ega,
               {
                 # Let user know
                 showNotification("Estimating EGA...")
 
                 ## Grab network estimation method
                 network <<- switch(input$model,
                                    "Graphical LASSO (glasso)" = "glasso",
                                    "Triangulated Maximally Filtered Graph (TMFG)" = "TMFG"
                 )
                 
                 ## Grab community detection algorihtm
                 algo <<- switch(input$algorithm,
                                 "Walktrap" = "walktrap",
                                 "Louvain" = "louvain"
                 )
                 
                 ## Check for number of cases
                 if(!exists("cases")){
                   cases <<- input$n
                 }
                 
                 ## Run EGA
                 ega.res <<- EGA(
                   data = dat, n = cases,
                   uni = as.logical(input$unidimensional),
                   model = network, algorithm = algo,
                   plot.EGA = FALSE
                 )
                 
                 ## Prepare citations
                 
                 ## Citation list
                 cites_list <<- list()
                 
                 ### Header
                 header <<- "<b>Please cite:</b><br>"
                 
                 ### General citations
                 cites_list$golino2017ega <- 'Golino, H., & Epskamp, S. (2017). Exploratory graph analysis: A new approach for estimating the number of dimensions in psychological research. <em>PLoS ONE</em>, <em>12</em>, e0174035. <a href="https://doi.org/10.1371/journal.pone.0174035">https://doi.org/10.1371/journal.pone.0174035</a>'
                 cites_list$golino2020ega <- 'Golino, H., Shi, D., Christensen, A. P., Garrido, L. E., Nieto, M. D., Sadana, R., & Thiyagarajan, J. A. (2020). Investigating the performance of Exploratory Graph Analysis and traditional techniques to identify the number of latent factors: A simulation and tutorial. <em>Psychological Methods</em>, <em>25</em>, 292-320. <a href="https://doi.org/10.1037/met0000255">https://doi.org/10.1037/met0000255</a>'
                 
                 ### Model and algorithm specific citations
                 #### Model
                 if(network == "glasso"){
                   
                   cites_list$epskamp2018glasso <- 'Epskamp, S., & Fried, E. I. (2018). A tutorial on regularized partial correlation networks. <em>Psychological Methods</em>, <em>23</em>, 617-634. <a href="https://doi.org/10.1037/met0000167">https://doi.org/10.1037/met0000167</a>'
                   cites_list$friedman2008glasso <- 'Friedman, J., Hastie, T., & Tibshirani, R. (2008). Sparse inverse covariance estimation with the graphical lasso. <em>Biostatistics</em>, <em>9</em>, 432-441. <a href="https://doi.org/10.1093/biostatistics/kxm045">https://doi.org/10.1093/biostatistics/kxm045</a>'
                   
                 }else if(network == "TMFG"){
                   
                   cites_list$massara2016tmfg <- 'Massara, G. P., Di Matteo, T., & Aste, T. (2016). Network filtering for big data: Triangulated maximally filtered graph. <em>Journal of Complex Networks</em>, <em>5</em>, 161-178. <a href="https://doi.org/10.1093/comnet/cnw015">https://doi.org/10.1093/comnet/cnw015</a>'
                   
                 }
                 
                 #### Algorithm
                 if(algo == "walktrap"){
                   
                   cites_list$pons2006walktrap <- 'Pons, P., & Latapy, M. (2006). Computing communities in large networks using random walks. <em>Journal of Graph Algorithms and Applications</em>, <em>10</em>, 191-218. <a href="https://doi.org/10.7155/jgaa.00185">https://doi.org/10.7155/jgaa.00185</a>'
                   
                 }else if(algo == "louvain"){
                   
                   cites_list$blondel2008louvain <- 'Blondel, V. D., Guillaume, J.-L., Lambiotte, R., & Lefebvre, E. (2008). Fast unfolding of communities in large networks. <em>Journal of Statistical Mechanics: Theory and Experiment</em>, <em>2008</em>, P10008. <a href="https://doi.org/10.1088/1742-5468/2008/10/P10008">https://doi.org/10.1088/1742-5468/2008/10/P10008</a>'
                   cites_list$christensen2020communities <- 'Christensen, A. P., & Golino, H. (2020). Estimating factors with psychometric networks: A Monte Carlo simulation comparing community detection algorithms. <em>PsyArXiv</em>. <a href="https://doi.org/10.31234/osf.io/hz89e">https://doi.org/10.31234/osf.io/hz89e</a>'
                   
                 }
                 
                 ## Formatted EGA citations
                 format_ega_cites <<- paste(cites_list[order(names(cites_list))], collapse = "<br><br>")
                 
                 # Output citation list
                 output$ega_cite <- renderUI({
                   
                   HTML(
                     
                     paste(header, format_ega_cites, collapse = "")
                     
                   )
                   
                 })
                 
                 ## Render semantic networks plot
                 output$ega_plot <- renderPlot({
                   plot(ega.res, plot.type = input$plot_type)
                 })
                 
                 # Show buttons
                 shinyjs::show("save_ega")
                 shinyjs::show("reset")
                 shinyjs::show("print_ega")
                 
               }
  )
  
  # Print EGA Methods section
  observeEvent(input$print_ega,
               {
                 EGA.methods.section(ega.res)
               })
  
  # Reset
  observeEvent(input$reset,
               {
                 
                 shinyalert::shinyalert(title = "Are you sure?",
                                        text = "You are about to erase your results\n(Data and saved results will not be erased)",
                                        type = "error",
                                        showConfirmButton = TRUE,
                                        showCancelButton = TRUE,
                                        callbackR = function(x)
                                        {
                                          if(x)
                                          {
                                            showNotification("Results cleared")
                                            
                                            # Refresh tables and plots
                                            output$ega_plot <- renderPlot({})
                                            
                                            # Network Estimation tab
                                            updateSelectInput(session = session,
                                                              inputId = "model",
                                                              label = "Network Estimation Method",
                                                              choices = c("Graphical LASSO (glasso)", 
                                                                          "Triangulated Maximally Filtered Graph (TMFG)")
                                            )
                                            
                                            updateSelectInput(session = session,
                                                              inputId = "algorithm",
                                                              label = "Community Detection Algorithm",
                                                              choices = c("Walktrap",
                                                                          "Louvain")
                                            )
                                            
                                            updateSelectInput(session = session,
                                                              inputId = "unidimensional",
                                                              label = "Check for Unidimensionality?",
                                                              choices = c("TRUE",
                                                                          "FALSE")
                                            )
                                            
                                            updateSelectInput(session = session,
                                                              inputId = "plot_type",
                                                              label = "Plotting Package",
                                                              choices = c("GGally",
                                                                          "qgraph")
                                            )
                                          
                                            # Hide tabs
                                            hideTab(inputId = "tabs", target = "Save and Reset All Results")
                                            
                                            # Hide save buttons
                                            shinyjs::hide("save_ega")
                                            
                                          }
                                        })
               })
  
  
  # Save events
  ## Data
  #observeEvent(input$save_data,
  #             {
  #               
  #               # Allow user to type name for object
  #               shinyalert::shinyalert(
  #                 title = "Save Data",
  #                 text = "Enter name for object:",
  #                 type = "input",
  #                 callbackR = function(value){
  #                   
  #                   # Get name for object
  #                   res.name <<- value
  #                   
  #                   # Add name to previous environment so it's not removed
  #                   prev.env <<- c(prev.env, res.name)
  #                   
  #                   # Create list
  #                   saveList <<- list()
  #                   
  #                   if(exists("dat", envir = globalenv()))
  #                   {saveList$data <<- dat}
  #                   
  #                   # Assign save list to result name
  #                   assign(
  #                     x = res.name,
  #                     value = saveList,
  #                     envir = globalenv()
  #                   )
  #                   
  #                   # Let user know save was successful
  #                   shinyalert::shinyalert(
  #                     title = "Save Successful",
  #                     text = paste("Data was saved as '", res.name, "'", sep = ""),
  #                     type = "info"
  #                   )
  #                 }
  #               )
  #               
  #             }
  #)
  
  ## Networks
  observeEvent(input$save_ega,
               {
                 
                 # Allow user to type name for object
                 shinyalert::shinyalert(
                   title = "Save Networks",
                   text = "Enter name for object:",
                   type = "input",
                   callbackR = function(value){
                     
                     # Get name for object
                     res.name <<- value
                     
                     # Add name to previous environment so it's not removed
                     prev.env <<- c(prev.env, res.name)
                     
                     # Create list
                     saveList <<- list()
                     
                     if(exists("ega.res", envir = globalenv()))
                     {saveList$EGA <<- ega.res}
                     
                     # Assign save list to result name
                     assign(
                       x = res.name,
                       value = saveList,
                       envir = globalenv()
                     )
                     
                     # Let user know save was successful
                     shinyalert::shinyalert(
                       title = "Save Successful",
                       text = paste("EGA results were saved as '", res.name, "'", sep = ""),
                       type = "info"
                     )
                   }
                 )
                 
               }
  )
  
  ## Master save
  observeEvent(input$save_master,
               {
                 
                 # Allow user to type name for object
                 shinyalert::shinyalert(
                   title = "Save All Results",
                   text = "Enter name for object:",
                   type = "input",
                   callbackR = function(value){
                     
                     # Get name for object
                     res.name <<- value
                     
                     # Add name to previous environment so it's not removed
                     prev.env <<- c(prev.env, res.name)
                     
                     # Create list
                     saveList <<- list()
                     
                     if(exists("dat", envir = globalenv()))
                     {saveList$data <<- dat}
                     
                     if(exists("ega.res", envir = globalenv()))
                     {saveList$EGA <<- ega.res}
                     
                     # Assign save list to result name
                     assign(
                       x = res.name,
                       value = saveList,
                       envir = globalenv()
                     )
                     
                     # Let user know save was successful
                     shinyalert::shinyalert(
                       title = "Save Successful",
                       text = paste("Spreading Activation Analyses were saved as '", res.name, "'", sep = ""),
                       type = "info"
                     )
                   }
                 )
                 
               }
  )
  
  
  
  onStop(function(x)
  {
    # Save results into condensed list
    resultShiny <<- list()
    
    if(exists("dat", envir = globalenv()))
    {resultShiny$data <<- dat}
  
    if(exists("ega.res", envir = globalenv()))
    {resultShiny$EGA <<- ega.res}
    
    # Remove all other variables from global environment
    rm(list = ls(envir = globalenv())[-match(c("resultShiny", prev.env), ls(globalenv()))], envir = globalenv())
    
    # Remove plots from user view
    if(!is.null(dev.list()))
    {dev.off()}
  }
  )
  
}