#'Calculate the Extent of Modification at the Peptide Level (Step 9)
#'
#' @param df_in The annotated data frame merged with the parsed FASTA file.
#'
#' @return a data frame containing the values for the extent of modification
#' calculation
#' @export
#'
#' @examples EOM<- area_calculations_pep(pd_data_annotated)
#' @aliases area_calculations_pep
area_calculations_pep <- function(df_in) {
  # Group by the necessary columns including 'Condition'
  df_out <- df_in %>%
    group_by(MasterProteinAccessions, Sequence, SampleControl, MOD, Condition) %>%
    reframe(TotalArea = sum(`Precursor Abundance`, na.rm = TRUE)) %>%
    ungroup() %>%
    pivot_wider(
      id_cols = c("MasterProteinAccessions", "Sequence", "Condition"),
      names_from = c("SampleControl", "MOD"),
      values_from = "TotalArea",
      values_fill = NA
    )

  # Rename columns
  df_out <- df_out %>%
    rename(
      Control_OxidizedArea = Control_Oxidized,
      Control_UnoxidizedArea = Control_Unoxidized,
      Sample_OxidizedArea = Sample_Oxidized,
      Sample_UnoxidizedArea = Sample_Unoxidized
    )


  # Calculate Total Areas
  df_out$TotalSampleArea <- rowSums(df_out[, c("Sample_OxidizedArea", "Sample_UnoxidizedArea")], na.rm = TRUE)
  df_out$TotalControlArea <- rowSums(df_out[, c("Control_OxidizedArea", "Control_UnoxidizedArea")], na.rm = TRUE)

  # Calculate EOM values
  df_out$EOMSample <- df_out$Sample_OxidizedArea / df_out$TotalSampleArea
  df_out$EOMControl <- df_out$Control_OxidizedArea / df_out$TotalControlArea
  df_out$EOM <- df_out$EOMSample - df_out$EOMControl

  # Calculate count (N)
  df_out <- df_out %>%
    left_join(df_in %>%
                group_by(MasterProteinAccessions, Sequence, Condition) %>%
                summarize(N = n()),
              by = c("MasterProteinAccessions", "Sequence", "Condition"))

  #Calculate Variance and Standard Variation
  # Group by the necessary columns including 'Condition'
  test <- df_in %>%
    group_by(MasterProteinAccessions, Sequence, SampleControl, MOD, Condition) %>%
    reframe(TotalVar = var(`Precursor Abundance`, na.rm = TRUE)) %>%
    ungroup() %>%
    pivot_wider(
      id_cols = c("MasterProteinAccessions", "Sequence", "Condition"),
      names_from = c("SampleControl", "MOD"),
      values_from = "TotalVar",
      values_fill = NA
    )

  # Rename columns
  test <- test %>%
    rename(
      Control_OxidizedVar = Control_Oxidized,
      Control_UnoxidizedVar = Control_Unoxidized,
      Sample_OxidizedVar = Sample_Oxidized,
      Sample_UnoxidizedVar = Sample_Unoxidized
    )

  df_out$Control_OxidizedVar<- test$Control_OxidizedVar
  df_out$Control_UnoxidizedVar<- test$Control_UnoxidizedVar
  df_out$Sample_OxidizedVar<- test$Sample_OxidizedVar
  df_out$Sample_UnoxidizedVar<- test$Sample_UnoxidizedVar
  df_out$TotalSampleVar<- rowSums(test[, c("Sample_OxidizedVar", "Sample_UnoxidizedVar")], na.rm = TRUE)
  df_out$TotalControlVar<- rowSums(test[, c("Control_OxidizedVar", "Control_UnoxidizedVar")], na.rm = TRUE)

  #Calculating Total Variance
  test2 <- df_in %>%
    group_by(MasterProteinAccessions, Sequence, SampleControl, Condition) %>%
    reframe(TotalArea = var(`Precursor Abundance`, na.rm = TRUE)) %>%
    ungroup() %>%
    pivot_wider(
      id_cols = c("MasterProteinAccessions", "Sequence", "Condition"),
      names_from = c("SampleControl"),
      values_from = "TotalArea",
      values_fill = NA
    )

  # Rename columns
  test2 <- test2 %>%
    rename(
      Sample_TotalVar = Sample,
      Control_TotalVar = Control,

    )
  df_out$Sample_TotalVar<- test2$Sample_TotalVar
  df_out$Control_TotalVar<- test2$Control_TotalVar


  # Corrected formulas for Variance calculations in R
  df_out$VarianceSample <- df_out$EOMSample^2 * (((df_out$Sample_OxidizedVar) / (df_out$Sample_OxidizedArea^2)) + (df_out$TotalSampleVar) / (df_out$TotalSampleArea^2))
  df_out$VarianceControl <- df_out$EOMControl^2 * (((df_out$Control_OxidizedVar) / (df_out$Control_OxidizedArea^2)) + (df_out$TotalControlVar) / (df_out$TotalControlArea^2))

  # Sum of Variances for Total Variance
  df_out$TotalVariance <- df_out$VarianceSample + df_out$VarianceControl

  # Calculation of Standard Deviation
  df_out$SD <- sqrt(df_out$TotalVariance)

  # Return the final dataframe
  return(df_out)
}
