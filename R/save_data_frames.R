#' save_data_frames
#'
#' @param file_path
#' @param output_directory
#' @param ...
#'
#' @return data frames saved in the desired file path
#' @export
#'
#' @examples save_data_frames(file_path, file_output,
#' TotalsTable = TotalsTable, quant_graph_df_pep = quant_graph_df_pep,
#' quant_graph_df_res = quant_graph_df_res)
save_data_frames <- function(file_path, output_directory, ...) {
  file_name <- basename(file_path)

  file_name <- tools::file_path_sans_ext(file_name)

  data_frames <- list(...)

  for (i in seq_along(data_frames)) {

    df <- data_frames[[i]]

    output_file_name <- paste0(file_name, "_", names(data_frames)[i])

    output_file_path <- file.path(output_directory, paste0(output_file_name, ".xlsx"))

    writexl::write_xlsx(df, output_file_path)

    cat(names(data_frames)[i], "saved as", output_file_path, "\n")
  }
}


