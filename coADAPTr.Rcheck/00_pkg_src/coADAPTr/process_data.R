#' process_data
#'
#' @return A data frame with the desired data set
#' @export
#'
#' @examples import<- process_data()
#' @aliases file_path
process_data <- function() {
  file_path <- file.choose()
  df <- read.xlsx(file_path)
  return(df)
}
