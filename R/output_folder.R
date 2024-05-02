#' output_folder
#' @param file_output The file path to save the output
#' @return The path to the folder where the output will be saved.
#' @export
#'
#' @examples output_folder()
#' @aliases output_folder
output_folder <- function() {
  file_output <- choose.dir()
  return(file_output)
}
