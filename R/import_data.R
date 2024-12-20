#' Import data from an Excel file (step 3 )
#'
#' @return A data frame with the desired data set
#' @export
#'
#' @examples import <- import_data()
#' @aliases import_data
import_data <- function() {
  library(openxlsx)


  cat("Please select the Excel or .tsv file that corresponds to your sequence searched data: \n")
  file_path <- file.choose()

  # Detecting file extension for proper data frame creation
  filename <- basename(file_path)
  split_result <- strsplit(filename, "\\.")[[1]]
  filename <- split_result[1]
  extension <- split_result[2]

  if(extension == "tsv"){
    print("Input file is tab-delim type")
    if(filename == "psm"){
      print("Analysing psm.tsv file (FragPipe-DDA)")
      df <- read.csv(file_path, sep = "\t",check.names = FALSE)
    }else if(grepl("report", filename, fixed = TRUE)){
      print("Analysing psm.tsv file (FragPipe-DIA)")
      df <- read.csv(file_path, sep = "\t",check.names = FALSE)
    }

  }else{
    print("Input file is Excel type")
    df <- read.xlsx(file_path, check.names = FALSE)
  }

  # Replace dots in column names with spaces
  colnames(df) <- gsub("\\.", " ", colnames(df))

  return(df)
}

