#' output_folder
#'
#' @return The path to the folder where the output will be saved.
#' @export
#'
#' @examples output_folder()
#' @aliases output
output_folder <- function() {
  file_output <- choose.dir()
  return(selected_folder)
}
