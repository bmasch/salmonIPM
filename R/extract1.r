#' Extracts fitted Stan object
#'
#' Improved function for extracting fitted Stan objects.
#' 
#' @param object fitted Stan object
#' @param par parameters/states to return
#' 
#' @export
extract1 <- function(object, par)
{
  extract(object, par)[[1]]
}
