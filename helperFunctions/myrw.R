
#' @export
#' Read tab delimited file
myread <- function(file, ...){
  return(read.table(file, header = TRUE, as.is = TRUE, sep = "\t", quote = "", comment = "", ...))
}

#' @export
#' Write tab delimited file
mywrite <- function(df, file, ...){
  write.table(df, file, sep = "\t", quote = FALSE, row.names = FALSE, ...)
}
