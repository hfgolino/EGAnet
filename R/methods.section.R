#' Automated Methods Section for \code{\link{EGAnet}} Objects
#' 
#' @description This function accepts \code{\link[EGAnet]{EGA}} objects
#' and generates a Methods section for your analysis. The output is
#' an HTML page containing the descriptions of the methods and parameters
#' as well as a Reference section for appropriate citation.
#' 
#' @param ... \code{\link{EGAnet}} objects.
#' Available methods (more methods will be added soon!):
#' 
#' \itemize{
#' 
#' \item{\code{\link[EGAnet]{EGA}}}
#' {Exploratory graph analysis}
#' 
#' \item{\code{\link[EGAnet]{bootEGA}}}
#' {Bootstrap exploratory graph analysis}
#' 
#' }
#' 
#' @param stats Methods section for statistics in \code{\link{EGAnet}}.
#' Multiple statistics can be input.
#' Available statistics:
#' 
#' \itemize{
#' 
#' \item{\code{\link[EGAnet]{net.loads}}}
#' {Network loadings. Requires \code{\link[EGAnet]{EGA}} object to be input}
#' 
#' \item{\code{\link[EGAnet]{net.scores}}}
#' {Network scores. Requires \code{\link[EGAnet]{EGA}} object to be input}
#' 
#' \item{\code{\link[EGAnet]{dimStability}}}
#' {Structural consistency. Requires \code{\link[EGAnet]{bootEGA}} object to be input}
#' 
#' \item{\code{\link[EGAnet]{itemStability}}}
#' {Item stability. Requires \code{\link[EGAnet]{bootEGA}} object to be input}
#' 
#' }
#'
#' @return Automated HTML Methods section in your default browser
#' 
#' @examples
#' \donttest{# Estimate EGA
#' ## plot.type = "qqraph" used for CRAN checks
#' ## plot.type = "GGally" is the default
#' ega.wmt <- EGA(data = wmt2[,7:24], plot.type = "qgraph")
#' }
#' # EGA Methods section
#' if(interactive()){
#' methods.section(ega.wmt)
#' }
#' 
#' \donttest{# Estimate standardized network loadings
#' wmt.loads <- net.loads(ega.wmt)$std
#' }
#' # EGA Methods section with network loadings
#' if(interactive()){
#' methods.section(ega.wmt, stats = "net.loads")
#' }
#' 
#' \donttest{# bootEGA example
#' ## plot.type = "qqraph" used for CRAN checks
#' ## plot.type = "GGally" is the default
#' boot.wmt <- bootEGA(data = wmt2[,7:24], iter = 500, plot.type = "qgraph",
#' type = "parametric", ncores = 2)
#' }
#' # EGA and bootEGA Methods section
#' if(interactive()){
#' methods.section(ega.wmt, boot.wmt)
#' }
#' 
#' \donttest{# Estimate structural consistency
#' sc.wmt <- dimStability(boot.wmt, orig.wc = ega.wmt$wc, item.stability = TRUE)
#' }
#' # EGA and bootEGA Methods section with structural consistency and item stability
#' if(interactive()){
#' methods.section(boot.wmt, stats = c("dimStability", "itemStability"))
#' }
#' 
#' # EGA with network loadings and
#' # bootEGA Methods section with structural consistency and item stability
#' if(interactive()){
#' methods.section(ega.wmt, boot.wmt, stats = c("net.loads", "dimStability", "itemStability"))
#' }
#' 
#' @importFrom utils packageVersion browseURL
#'
#' @export
#'
# Methods Section----
# Updated 15.02.2021
methods.section <- function(..., stats = c("net.loads", "net.scores",
                                           "dimStability", "itemStability"))
{
  # All statistics
  all.stats <- c("net.loads", "net.scores",
                 "dimStability", "itemStability",
                 "entropyFit", "tefi", "vn.entropy",
                 "ergoInfo", "mctest.ergoInfo")
  
  # Get stats
  if(missing(stats)){
    stats <- ""
  }else{
    stats <- match.arg(stats, several.ok = TRUE)
  }
  
  # All objects
  all.objects <- list(...)
  
  # All methods
  all.methods <- c("EGA", "EGA.fit", "bootEGA", "dynEGA", "LCT", "UVA")
  
  # Initialize Methods section
  section <- list()
  
  # Initialize ordering
  ordering <- list()
  
  # Loop through objects
  for(i in 1:length(all.objects)){
    
    # Check for EGA object
    if(!class(all.objects[[i]]) %in% all.methods){
      
      # Give warning that object is not found in EGA methods
      warning(paste("Object input", names(all.objects)[i], "is not an EGA methods class. This object has been skipped."))
      
    }else{
      
      # Get Methods section function
      section[[i]] <- switch(class(all.objects[[i]]),
                             "EGA" = EGA.methods.section(all.objects[[i]],
                                                         ifelse("net.loads" %in% stats, TRUE, FALSE),
                                                         ifelse("net.scores" %in% stats, TRUE, FALSE)),
                             "EGA.fit" = EGA.fit.methods.section(all.objects[[i]]),
                             "bootEGA" = bootEGA.methods.section(all.objects[[i]],
                                                                 ifelse("dimStability" %in% stats, TRUE, FALSE),
                                                                 ifelse("itemStability" %in% stats, TRUE, FALSE)),
                             "dynEGA" = dynEGA.methods.section(all.objects[[i]]),
                             "UVA" = UVA.methods.section(all.objects[[i]])
      )
      
      # Ensure ordering of sections
      ordering[[i]] <- switch(class(all.objects[[i]]),
                              "UVA" = 1,
                              "EGA" = 2,
                              "EGA.fit" = 2,
                              "dynEGA" = 2,
                              "bootEGA" = 3
      )
      
      # Change name of section
      names(section)[i] <- class(all.objects[[i]])
      
    }
  }
  
  # Unlist ordering
  ordering <- unlist(ordering)
  ordering <- order(ordering)
  
  # Order objects and sections
  all.objects <- all.objects[ordering]
  section <- section[ordering]
  
  # Markdown YAML
  YAML <- c("---",
            "title: EGA Methods Section",
            "output: html_document",
            "---"
  )
  
  # Year
  year <- unlist(strsplit(as.character(Sys.Date()), split = "\\-"))[1]
  
  # Version
  version <- packageVersion("EGAnet")
  
  # Version Year
  version.year <- packageDescription("EGAnet")
  version.year <- unlist(strsplit(as.character(version.year$Date), split = "\\-"))[1]
  
  # References
  refs <- list()
  
  refs$golino <- paste("Golino, H., & Christensen, A. P. (", version.year, "). ",
                        "EGAnet: Exploratory Graph Analysis -- A framework for estimating the number of dimensions in multivariate data using network psychometrics. ",
                        "Retrieved from https://cran.r-project.org/package=EGAnet",
                        sep = "")
  
  refs$r <- paste("R Core Team (", version.year, "). ",
                  "<em>R: A language and environment for statistical computing.</em> R Foundatin for Statistical Computing, Vienna, Austria. ",
                  "Retrieved from https://www.R-project.org",
                  sep = "")
  
  ## Data analysis section
  data.analysis.header <- "## Data Analysis"
  
  ## Data analysis text
  ### Methods and statistics
  if(length(section) == 1){
    
    data.analysis.text <- paste("&emsp;",
      names(section),
      " was applied using the *EGAnet* package ",
      "(version ", version, "; Golino & Christensen, ", version.year, ") in R (version ",
      paste(R.version$major, R.version$minor, sep = "."), "; R Core Team, ", year, ").",
      sep = "")
      
  }else if(length(section) == 2){
    
    data.analysis.text <- paste("&emsp;",
      names(section)[1], " and ", names(section)[2],
      " were applied using the *EGAnet* package ",
      "(version ", version, "; Golino & Christensen, ", version.year, ") in R (version ",
      paste(R.version$major, R.version$minor, sep = "."), "; R Core Team, ", year, ").",
      sep = "")
    
  }else{
    
    data.analysis.text <- paste("&emsp;",
      paste(paste(names(section)[1:(length(section) - 1)], sep = "", collapse = ", "), ", and ", names(section)[length(section)], sep = ""),
      " were applied using the *EGAnet* package ",
      "(version ", version, "; Golino & Christensen, ", version.year, ") in R (version ",
      paste(R.version$major, R.version$minor, sep = "."), "; R Core Team, ", year, ").",
      sep = "")
    
  }
  
  ### Plots
  plots <- numeric(length(all.objects))
  
  for(i in 1:length(all.objects)){
    
    # Names of objects (lowercase)
    name <- tolower(names(all.objects[[i]]))
    
    # Check for plot
    idx <- grep("plot", name)
    
    if(length(idx) != 0){
      plots[i] <- as.numeric(ggplot2::is.ggplot(all.objects[[i]][[idx]])) + 1
      # Add one for:
      # 0 = no plot
      # 1 = qgraph
      # 2 = ggplot2/GGally
    }
    
  }
  
  ## Check for plot names
  plot.names <- names(section[plots != 0])
  
  ## Plot references
  plot.refs <- list()
  
  ### Years
  GGally.version <- packageVersion("GGally")
  GGally.year <- unlist(strsplit(as.character(packageDescription("GGally")$Date), split = "\\-"))[1]
  ggplot2.version <- packageVersion("ggplot2")
  ggplot2.year <- unlist(strsplit(as.character(packageDescription("ggplot2")$Date), split = "\\-"))[1]
  qgraph.version <- packageVersion("qgraph")
  qgraph.year <- unlist(strsplit(as.character(packageDescription("qgraph")$Date), split = "\\-"))[1]
  
  ## Check for all ggplot2
  if(all(plots[plots != 0] != 2)){
    
    ### Plots
    if(length(plot.names) == 1){
      
      plot.text <- paste(plot.names,
                         " and associated results were visualized using the *GGally* ",
                         "(version ", GGally.version, "; Schloerke et al., ", GGally.year, "), ",
                         "*ggplot2* ", "(version ", ggplot2.version, "; Wickham, ", ggplot2.year, "), ",
                         "and *qgraph* ", "(version ", qgraph.version, "; Epskamp et al., 2012) ",
                         "packages in R.",
                         sep = "")
      
      
    }else{
      
      plot.text <- paste(paste(plot.names, collapse = ", "),
                         ", and associated results were visualized using the *GGally* ",
                         "(version ", GGally.version, "; Schloerke et al., ", GGally.year, "), ",
                         "*ggplot2* ", "(version ", ggplot2.version, "; Wickham, ", ggplot2.year, "), ",
                         "and *qgraph* ", "(version ", qgraph.version, "; Epskamp et al., 2012) ",
                         "packages in R.",
                         sep = "")
      
    }
    
    plot.refs$schloerke <- paste("Schloerke, Cook, Larmarange, Briatte, & Marbach (", GGally.year, "). ",
                                 "GGally: Extention to 'ggplot2'. ",
                                 "Retrieved from https://cran.r-project.org/package=GGally",
                                 sep = "")
    
    plot.refs$wickham <- paste("Wickham, H. (", ggplot2.year, "). ",
                               "*ggplot2: Elegant graphics for data analysis.* ",
                               "New York, NY: Springer. ",
                               "Retrieved from https://ggplot2-book.org",
                               sep = "")
    
    data.analysis.text <- paste(data.analysis.text, plot.text, sep = " ")
    
  }else{
    
    ### Plots
    if(length(plot.names) == 1){
      
      plot.text <- paste(plot.names,
                         " and associated results were visualized using the *GGally* ",
                         "(version ", GGally.version, "; Schloerke et al., ", GGally.year, ") and ",
                         "*ggplot2* ", "(version ", ggplot2.version, "; Wickham, ", ggplot2.year, ") ",
                         "packages in R.",
                         sep = "")
      
      
    }else{
      
      plot.text <- paste(paste(plot.names, collapse = ", "),
                         ", and associated results were visualized using the *GGally* ",
                         "(version ", GGally.version, "; Schloerke et al., ", GGally.year, ") and ",
                         "*ggplot2* ", "(version ", ggplot2.version, "; Wickham, ", ggplot2.year, ") ",
                         "packages in R.",
                         sep = "")
      
    }
    
    plot.refs$schloerke <- paste("Schloerke, B., Cook, D., Larmarange, J., Briatte, F., & Marbach, M. (", GGally.year, "). ",
                                 "GGally: Extention to 'ggplot2'. ",
                                 "Retrieved from https://cran.r-project.org/package=GGally",
                                 sep = "")
    
    plot.refs$wickham <- paste("Wickham, H. (", ggplot2.year, "). ",
                               "*ggplot2: Elegant graphics for data analysis.* ",
                               "New York, NY: Springer. ",
                               "Retrieved from https://ggplot2-book.org",
                               sep = "")
    
    data.analysis.text <- paste(data.analysis.text, plot.text, sep = " ")
    
  }
  
  ## References
  references.header <- "## References"
  
  ## Get all other references
  for(i in 1:length(section)){
    refs <- c(refs, section[[i]]$references)
  }
  
  ## Add in plotting references
  refs <- c(refs, plot.refs)
  
  ## Order alphabetically
  references.text <- refs[order(names(refs))]
  
  references.text <- paste(references.text, collapse = "\n\n")
  
  ### Loop through sections
  sections <- paste(lapply(section, function(x){x$text}), collapse = "\n")
  
  markobj <- paste(YAML,
                   sections,
                   data.analysis.header, data.analysis.text,
                   references.header, references.text,
                   sep = "\n")
  
  tempHTML <- paste(tempdir(), "EGAnet_method.html", sep = "\\")
  tempHTML <- gsub("\\\\", "/", tempHTML)
  
  markdown::markdownToHTML(text = knitr::knit(text = markobj), output = tempHTML)
  
  browseURL(tempHTML)
  
}
