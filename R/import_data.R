#' Import data from an Excel file (step 3 )
#'
#' @return A data frame with the desired data set
#' @export
#'
#' @examples import <- import_data()
#' @aliases import_data
import_data <- function() {
  library(openxlsx)

  file_path <- file.choose()
  df <- read.xlsx(file_path, check.names = FALSE)

  colnames(df) <- gsub("\\.", " ", colnames(df))

  return(df)
}
