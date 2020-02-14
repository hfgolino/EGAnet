#' Changes Variable Names to Descriptions for \code{\link[EGAnet]{node.redundant}} Objects
#' @description Using a key, this function changes the variable names in
#' the \code{\link[EGAnet]{node.redundant}} output to descriptions
#'
#' @param node.redundant.obj A \code{\link[EGAnet]{node.redundant}} object
#'
#' @param key Character vector.
#' A vector with variable descriptions that correspond
#' to the order of variables from the data used as input into the
#' \code{\link[EGAnet]{node.redundant}} function
#'
#' @return Returns a list:
#'
#' \item{redundant}{Vectors nested within the list corresponding
#' to redundant nodes with the name of object in the list}
#'
#' \item{data}{Returns original data}
#'
#' \item{weights}{Returns weights determine by weighted topological overlap
#' or partial correlations}
#'
#' \item{key}{Returns original key}
#'
#' @examples
#' # obtain SAPA items
#' items <- psychTools::spi[,-c(1:10)]
#'
#' \donttest{
#' # weighted topological overlap
#' redund <- node.redundant(items, method = "wTO", type = "adapt")
#'
#' # partial correlation
#' redund <- node.redundant(items, method = "pcor", type = "adapt")
#'
#' # check redundancies
#' key.ind <- match(colnames(items), as.character(psychTools::spi.dictionary$item_id))
#' key <- as.character(psychTools::spi.dictionary$item[key.ind])
#'
#' # change names in redundancy output to questionnaire item description
#' named.nr <- node.redundant.names(redund, key)
#' }
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
# Node Redundant Naming Function
# Updated 13.02.2020
node.redundant.names <- function(node.redundant.obj, key)
{
  # Check for node.redundant object class
  if(class(node.redundant.obj) != "node.redundant")
  {stop("A 'node.redundant' object must be used as input")}

  # Obtain and remove data from node redundant object
  data <- node.redundant.obj$data

  # Check that columns match key
  if(ncol(data) != length(as.vector(key)))
  {stop("Number of columns in data does not match the length of 'key'")}

  # Names of node.redundant object
  nr.names <- names(node.redundant.obj$redundant)

  # Key names
  key.names <- colnames(data)

  # Key change
  key.chn <- key

  for(i in 1:length(nr.names))
  {
    # Target redundant node
    target.r <- match(names(node.redundant.obj$redundant)[i],key.names)

    # Replace item name with description
    names(node.redundant.obj$redundant)[i] <- as.character(key.chn[target.r])

    # Target other nodes
    target.o <- match(node.redundant.obj$redundant[[i]],key.names)

    # Replace item names with description
    node.redundant.obj$redundant[[i]] <- as.character(key.chn[target.o])
  }

  names(key) <- colnames(data)
  node.redundant.obj$key <- key

  return(node.redundant.obj)
}
#----
