#' locate_startend_res
#'
#' @param raw_data data frame with original data to map peptide locations on
#' protein sequences
#'
#' @return A data frame with with the start and ending peptide location
#' @export
#'
#' @examples new_df <- locate_startend_res(raw_data)
#' @aliases locate_startend_res
locate_startend_res <- function(raw_data){

  uniqueMPA <- unique(raw_data[, c("Master Protein Accessions",'Sequence')])
  uniqueMPA <- as.data.frame(uniqueMPA)

  raw_data <- merge(raw_data, FASTA, by = "UniprotID")

  index <- str_locate(raw_data$protein_sequence, raw_data$Sequence)
  raw_data <- cbind(raw_data, index)  # TODO: memory efficiency?

  raw_data$peptide<- paste(raw_data$start,"-", raw_data$end)

  raw_data$mod_count <- str_count(raw_data$Modifications, "\\(.*?\\)")
  raw_data$mod_count <- ifelse(raw_data$MOD == "Unoxidized", 0, raw_data$mod_count)

  raw_data$mod_res <- ifelse(raw_data$ModPositionN>0, raw_data$start + raw_data$ModPositionN - 1, NA)

  raw_data$Res<- paste(raw_data$ModPositionL, raw_data$mod_res)

  return(raw_data)
}
