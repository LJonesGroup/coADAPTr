process_data <- function() {
  # Ask user to select the file
  file_path <- file.choose()

  # Your data processing code here
  # For example, reading the Excel file
  df <- read.xlsx(file_path)  # Reading Excel file

  # Additional data processing steps can be added

  return(df)  # Returning processed data
}
