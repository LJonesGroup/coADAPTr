#' import_data
#'
#' @return A data frame with the desired data set
#' @export
#'
#' @examples import<- import_data()
#' @aliases import_data
import_data <- function() {
  file_path <- file.choose()
  df <- read.xlsx(file_path)
  return(df)
}
