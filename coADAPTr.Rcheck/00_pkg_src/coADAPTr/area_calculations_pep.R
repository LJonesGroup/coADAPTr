#'area_calculations_pep
#'
#' @param df_in A dataframe with the following columns: MasterProteinAccessions, Sequence, SampleControl, MOD, Precursor Abundance
#'
#' @return a dataframe containing the values for the extent of modification calculation
#' @export
#'
#' @examples
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
      Control_OxidizedArea = Control_Oxidized,
      Control_UnoxidizedArea = Control_Unoxidized,
      Sample_OxidizedArea = Sample_Oxidized,
      Sample_UnoxidizedArea = Sample_Unoxidized
    )

  df_out$Control_OxidizedArea <- ifelse(df_out$Control_UnoxidizedArea > 0 & is.na(df_out$Control_OxidizedArea), 0, df_out$Control_OxidizedArea)

  df_out$Sample_UnoxidizedArea <- ifelse(df_out$Sample_OxidizedArea > 0 & is.na(df_out$Sample_UnoxidizedArea), 0, df_out$Sample_UnoxidizedArea)


  df_out$TotalSampleArea <- rowSums(df_out[, c("Sample_OxidizedArea", "Sample_UnoxidizedArea")])
  df_out$TotalControlArea <- rowSums(df_out[, c("Control_OxidizedArea", "Control_UnoxidizedArea")])


  df_out$EOMSample <- df_out$Sample_OxidizedArea / df_out$TotalSampleArea
  df_out$EOMControl <- df_out$Control_OxidizedArea / df_out$TotalControlArea
  df_out$EOM <- df_out$EOMSample - df_out$EOMControl

  df_out$N <- df_in %>%
    group_by(MasterProteinAccessions, Sequence) %>%
    count()

  df_out$N$Sequence <- NULL
  df_out$N$MasterProteinAccessions <- NULL
  colnames(df_out)[12] <- "N"


  test <- pd_data_fasta_merged %>%
    group_by(MasterProteinAccessions, Sequence) %>%
    reframe(sdprep = sd(`Precursor Abundance`))
  df_out$sdprep<- test$sdprep
  df_out$SD<- df_out$SD <- df_out$sdprep/(df_out$TotalSampleArea+df_out$TotalControlArea)


  return(df_out)
}
