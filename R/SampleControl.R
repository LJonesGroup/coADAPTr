#' Identify Sample (Laser) or Control (No Laser) Files (Step4)
#'
#' @param pd_data Sequence searched data to which a column can be added to indicate
#' Sample vs Control files by referencing the name indicated in the spectrum file.
#' Ensure that laser irradiated samples are indicated by "L" in the spectrum file.
#' Meanwhile control samples are indicated by "NL" in the spectrum file.
#'
#' @return A column indicating Sample vs Control
#' @export
#'
#' @examples pd_data <- SampleControl(pd_data)
#' @aliases SampleControl
SampleControl <- function(pd_data) {
  pd_data$SampleControl <- ifelse(grepl("NL", pd_data$`Spectrum File`), "Control", "Sample")
  return(pd_data)
}
