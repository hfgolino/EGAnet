compare.EGA.plots <- function(..., rows, columns)
{
  # Obtain object list
  object.list <- list(...)
  
  # Identify EGA objects
  EGA.idx <- grep("EGA", unlist(lapply(object.list, class)))
  
  # Obtain EGA objects only
  object.list <- object.list[EGA.idx]
  
  # Check for at least two objects
  if(length(object.list) < 2){
    stop("EGA plot comparisons require two or more EGA objects")
  }
  
  # Initialize plot list
  plots.ega <- list()
  plots.ega[[1]] <- plot(object.list[[1]])
  
  # Loop through matching 
  for(i in 2:length(object.list)){
    plots.ega[[i]] <- compare.EGA(object.list[[1]], object.list[[i]])[[2]]
  }
  
  # Set up grid return
  ggpubr::ggarrange(plotlist = plots.ega, ncol = columns, nrow = rows)
  
}