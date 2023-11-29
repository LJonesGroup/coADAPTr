FASTA_file <- function() {
  # Ask user to select the file
  FASTA_path <- file.choose()

  # Your data processing code here
  # For example, reading the Excel file
  FASTA <- read.fasta(FASTA_path)  # Reading Excel file

  # Additional data processing steps can be added

  return(FASTA)  # Returning processed data
}

