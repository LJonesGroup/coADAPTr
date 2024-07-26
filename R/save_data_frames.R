#' Save the Data Frames from your Analysis (Step 16B)
#'
#' @param output_directory the directory to save the data frames
#' @param ... data frames to save
#'
#' @return data frames saved in the desired file path
#' @export
#'
#' @examples save_data_frames(file_path, file_output,
#' TotalsTable = TotalsTable, quant_graph_df_pep = quant_graph_df_pep,
#' quant_graph_df_res = quant_graph_df_res)
#' @aliases save_data_frames
#'
save_data_frames <- function(output_directory, ...) {
  # Prompt the user to input the file name
  cat("Enter the file name for the resulting tables of data that were generated during analysis (without extension): ")
  file_name <- readline()

  # Remove any leading or trailing whitespace
  file_name <- trimws(file_name)

  # Check if the file name is empty, if so, use a default name
  if (nchar(file_name) == 0) {
    file_name <- "output"
  }

  # Create a list of data frames
  data_frames <- list(...)

  # Iterate over each data frame and save as separate Excel files
  for (i in seq_along(data_frames)) {
    # Get the current data frame
    df <- data_frames[[i]]

    # Create the output file name
    output_file_name <- paste0(file_name, "_", names(data_frames)[i])

    # Create the output file path in the output directory
    output_file_path <- file.path(output_directory, paste0(output_file_name, ".xlsx"))

    # Save the data frame as an Excel file
    writexl::write_xlsx(df, output_file_path)

    # Print a message to indicate successful saving
    cat(names(data_frames)[i], "saved as", output_file_path, "\n")
  }
}


