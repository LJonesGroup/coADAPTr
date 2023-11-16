area_calculations_pep <- function(df_in) {
  df_out <- df_in %>%
    group_by(MasterProteinAccessions, Sequence, SampleControl, MOD) %>%
    reframe(TotalArea = sum(`Precursor Abundance`) ) %>%
    ungroup() %>%
    pivot_wider(
      id_cols = c("MasterProteinAccessions", "Sequence"),
      names_from = c("SampleControl", "MOD"),
      values_from = "TotalArea",
      values_fill = NA
    )

  df_out <- df_out %>%
    rename(
      Control_OxidizedArea = Control_Oxidized,    # Rename OldColumn1 to NewColumn1
      Control_UnoxidizedArea = Control_Unoxidized,    # Rename OldColumn2 to NewColumn2
      Sample_OxidizedArea = Sample_Oxidized,     # Rename OldColumn3 to NewColumn3
      Sample_UnoxidizedArea = Sample_Unoxidized     # Rename OldColumn3 to NewColumn3
    )

  #create conditions so that the bare minimum criteria for FPOP can be calculated
  #specifically cases where sample oxidized area is detected and control oxidized area re detected.
  #Cases where that peptide was not detected in the sample unoxidized area and control oxidized are will be turned to 0
  ###*****MAKE SURE THIS IS OK WITH LISA WANT TO AVOID IMPLICIT BIAS
  #The peptide at least has to be detected by the mass spec in order to calculate control total area?
  df_out$Control_OxidizedArea <- ifelse(df_out$Control_UnoxidizedArea > 0 & is.na(df_out$Control_OxidizedArea), 0, df_out$Control_OxidizedArea)

  df_out$Sample_UnoxidizedArea <- ifelse(df_out$Sample_OxidizedArea > 0 & is.na(df_out$Sample_UnoxidizedArea), 0, df_out$Sample_UnoxidizedArea)

  #Calculating the extent of modification (EOM) for the sample and control data for each peptide

  df_out$TotalSampleArea <- rowSums(df_out[, c("Sample_OxidizedArea", "Sample_UnoxidizedArea")])
  df_out$TotalControlArea <- rowSums(df_out[, c("Control_OxidizedArea", "Control_UnoxidizedArea")])


  #Calculating the extent of modification (EOM) by subtracting background oxidation (control)
  df_out$EOMSample <- df_out$Sample_OxidizedArea / df_out$TotalSampleArea
  df_out$EOMControl <- df_out$Control_OxidizedArea / df_out$TotalControlArea
  df_out$EOM <- df_out$EOMSample - df_out$EOMControl

  #Calculating N or the number of times the petpdide was detected by the MS
  df_out$N <- df_in %>%
    group_by(MasterProteinAccessions, Sequence) %>%
    count()

  df_out$N$Sequence <- NULL
  df_out$N$MasterProteinAccessions <- NULL
  colnames(df_out)[12] <- "N"

  #Calcculating standard deviation
  test <- pd_data_fasta_merged %>%
    group_by(MasterProteinAccessions, Sequence) %>%
    reframe(sdprep = sd(`Precursor Abundance`))
  df_out$sdprep<- test$sdprep
  df_out$SD<- df_out$SD <- df_out$sdprep/(df_out$TotalSampleArea+df_out$TotalControlArea)

  # Renaming the 12th column to "N"


  return(df_out)
}
