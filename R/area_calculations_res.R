#' Calculate the Extent of Modification at the Residue Level (Step 13)
#'
#' @param df_in The annotated data frame merged with the parsed FASTA file.
#'
#' @return a data frame containing the extent of modification calculations
#' at the residue level
#' @export
#'
#' @examples EOM<- area_calculations_res(pd_data)
#' @aliases area_calculations_res
area_calculations_res <- function(df_in) {

  # Filter for mod_count 0 or 1
  df_out <- df_in %>%
    filter(mod_count == 0 | mod_count == 1) %>%
    group_by(MasterProteinAccessions, Sequence, SampleControl, MOD, Condition) %>%
    summarize(TotalArea = sum(`Precursor Abundance`, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      id_cols = c("MasterProteinAccessions", "Sequence", "Condition"),
      names_from = c("SampleControl", "MOD"),
      values_from = "TotalArea",
      values_fill = NA
    )

  # Calculate total areas for Sample and Control
  df_out <- df_out %>%
    mutate(SampleTotalArea = coalesce(Sample_Oxidized, 0) + coalesce(Sample_Unoxidized, 0),
           ControlTotalArea = coalesce(Control_Oxidized, 0) + coalesce(Control_Unoxidized, 0))

  # Remove individual oxidized/unoxidized columns
  df_out <- df_out %>%
    select(-Sample_Oxidized, -Sample_Unoxidized, -Control_Unoxidized, -Control_Oxidized)

  # Filter and summarize oxidized areas
  df_out2 <- df_in %>%
    filter((mod_count == 0 | mod_count == 1) & MOD == "Oxidized") %>%
    group_by(MasterProteinAccessions, Sequence, Res, SampleControl, Condition) %>%
    summarize(OxidizedArea = sum(`Precursor Abundance`, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      id_cols = c("MasterProteinAccessions", "Sequence", "Res", "Condition"),
      names_from = c("SampleControl"),
      values_from = "OxidizedArea",
      values_fill = NA
    )

  # Merge the two dataframes
  df_out <- full_join(df_out, df_out2, by = c("MasterProteinAccessions", "Sequence", "Condition"))

  # Rename columns
  df_out <- df_out %>%
    rename(SampleOxidizedArea = Sample,
           ControlOxidizedArea = Control)

  # Calculate EOM values
  df_out <- df_out %>%
    mutate(EOMSample = SampleOxidizedArea / SampleTotalArea,
           EOMControl = ControlOxidizedArea / ControlTotalArea,
           EOM = EOMSample - EOMControl)

  # Calculate count (N)
  N_df <- df_in %>%
    filter(mod_count == 0 | mod_count == 1, !is.na(Res)) %>%
    group_by(MasterProteinAccessions, Sequence, Res) %>%
    summarize(N = n(), .groups = "drop")  # Count the occurrences

  df_out <- df_out %>%
    left_join(N_df, by = c("MasterProteinAccessions", "Sequence", "Res"))

  # Calculate standard deviation
  #Calculate the oxidized variation
  ox_var<- df_in %>%
    filter((mod_count == 0 | mod_count == 1) & MOD == "Oxidized") %>%
    group_by(MasterProteinAccessions, Sequence, Res, SampleControl, Condition) %>%
    summarize(OxidizedVar = var(`Precursor Abundance`, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      id_cols = c("MasterProteinAccessions", "Sequence", "Res", "Condition"),
      names_from = c("SampleControl"),
      values_from = "OxidizedVar",
      values_fill = NA
    )
  # Rename columns
  ox_var <- ox_var %>%
    rename(
      Control_OxidizedVar = Control,
      Sample_OxidizedVar = Sample,
    )

  #Calculating Total Vairance
  total_var <- df_in %>%
    filter((mod_count == 0 | mod_count == 1)) %>%
    group_by(MasterProteinAccessions, Sequence, SampleControl, Condition) %>%
    summarize(TotalVar = var(`Precursor Abundance`, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      id_cols = c("MasterProteinAccessions", "Sequence", "Condition"),
      names_from = c("SampleControl"),
      values_from = "TotalVar",
      values_fill = list(TotalVar = NA)
    ) %>%
    rename(
      SampleTotalVar = Sample,
      ControlTotalVar = Control
    )



  merged<- inner_join(total_var, ox_var , by = c("MasterProteinAccessions", "Sequence", "Condition"))
  #merge variance data with df_out
  df_out<- inner_join(df_out, merged, by = c("MasterProteinAccessions", "Sequence", "Res", "Condition"))
  # Corrected formulas for Variance calculations in R
  df_out$VarianceSample <- df_out$EOMSample^2 * (((df_out$Sample_OxidizedVar) / (df_out$SampleOxidizedArea^2)) + (df_out$SampleTotalVar) / (df_out$SampleTotalArea^2))
  df_out$VarianceControl <- df_out$EOMControl^2 * (((df_out$Control_OxidizedVar) / (df_out$ControlOxidizedArea^2)) + (df_out$ControlTotalVar) / (df_out$ControlTotalArea^2))

  # Sum of Variances for Total Variance
  df_out$TotalVariance <- df_out$VarianceSample + df_out$VarianceControl

  # Calculation of Standard Deviation
  df_out$SD <- sqrt(df_out$TotalVariance)
  # Filter out rows with missing Res or EOM values
  #df_out <- df_out[complete.cases(df_out[c("Res", "EOM")]), ]

  return(df_out)
}
