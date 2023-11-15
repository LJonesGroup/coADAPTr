###########################################################################
#The function aims to ocate the start and end position of the Modification on the protein sequence
locate_startend_res <- function(raw_data){

  uniqueMPA <- unique(raw_data[, c("Master Protein Accessions",'Sequence')])
  uniqueMPA <- as.data.frame(uniqueMPA)

  # merge MPA and fasta dataframes together
  raw_data <- merge(raw_data, FASTA, by = "UniprotID")

  # get index of peptide sequence from protein sequence
  index <- str_locate(raw_data$protein_sequence, raw_data$Sequence)
  raw_data <- cbind(raw_data, index)  # TODO: memory efficiency?

  #annotate peptide
  raw_data$peptide<- paste(raw_data$start,"-", raw_data$end)

  raw_data$mod_count <- str_count(raw_data$Modifications, "\\(.*?\\)")
  raw_data$mod_count <- ifelse(raw_data$MOD == "Unoxidized", 0, raw_data$mod_count)

  # get modified residue and number
  raw_data$mod_res <- ifelse(raw_data$ModPositionN>0, raw_data$start + raw_data$ModPositionN - 1, NA)

  #combine the mod AA and res
  raw_data$Res<- paste(raw_data$ModPositionL, raw_data$mod_res)


  return(raw_data)
}
