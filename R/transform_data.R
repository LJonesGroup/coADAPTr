

#' Transform the Horizontal TMT Data to a Linear Format
#'
#' @param data
#'
#' @return A data frame that has transformed the input data to a linear format
#' @export
#'
#' @examples transformed_data <- transform_data(renamed_data)
#' @aliases transform_data
transform_data <- function(data) {
  # Pivot the data into long format, splitting columns containing ":"
  transformed_data <- pivot_longer(data, cols = contains(":"), names_to = "Spectrum File", values_to = "Precursor Abundance")

  return(transformed_data)
}

