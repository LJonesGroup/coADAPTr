#' Import data from an Excel file (step 3 )
#'
#' @return A data frame with the desired data set
#' @export
#'
#' @examples import <- import_data()
#' @aliases import_data
import_data <- function() {
  library(openxlsx)

  cat("Please select the excel file that corresponds to your sequence searched data: \n")
  file_path <- file.choose()
  df <- read.xlsx(file_path, check.names = FALSE)

  # Replace dots in column names with spaces
  colnames(df) <- gsub("\\.", " ", colnames(df))

  return(df)
}
