#' Calculate the Extent of Modification at the Residue Level (Step 13)
#' @param df_in data frame from modified raw data to calculate the EOM
#' @return a data frame containing the extent of modification calculations
#' at the residue level
#' @export
#'
#' @examples EOM<- area_calculations_res(pd_data)
#' @aliases area_calculations_res
area_calculations_res <- function(df_in) {

  df_out <- df_in %>%
    filter(mod_count == 0 | mod_count == 1) %>%
    group_by(MasterProteinAccessions, Sequence, SampleControl, MOD) %>%
    reframe(TotalArea = sum(`Precursor Abundance`)) %>%
    ungroup() %>%
    pivot_wider(
      id_cols = c("MasterProteinAccessions", "Sequence"),
      names_from = c("SampleControl", "MOD"),
      values_from = "TotalArea",
      values_fill = NA
    )


  df_out<- df_out %>%
    mutate(SampleTotalArea = Sample_Oxidized + Sample_Unoxidized,
           ControlTotalArea = Control_Oxidized + Control_Unoxidized)

  df_out <- df_out %>%
    select(-Sample_Oxidized, -Sample_Unoxidized, -Control_Unoxidized, -Control_Oxidized)


  df_out2 <- df_in %>%
    filter((mod_count == 0 | mod_count == 1) & MOD == "Oxidized")  %>%
    group_by(MasterProteinAccessions, Sequence, Res, SampleControl) %>%
    reframe(OxidizedArea = sum(`Precursor Abundance`)) %>%
    ungroup() %>%
    pivot_wider(
      id_cols = c("MasterProteinAccessions", "Sequence", "Res"),
      names_from = c("SampleControl"),
      values_from = "OxidizedArea",
      values_fill = NA
    )

  df_out<- full_join(df_out, df_out2, by = c("MasterProteinAccessions", "Sequence"))


  df_out <- df_out %>%
    rename(SampleOxidizedArea = Sample,
           ControlOxidizedArea = Control)




  df_out$EOMSample <- df_out$SampleOxidizedArea / df_out$SampleTotalArea
  df_out$EOMControl <- df_out$ControlOxidizedArea / df_out$ControlTotalArea
  df_out$EOM <- df_out$EOMSample - df_out$EOMControl

  ##################
  #########################

  N_df <- df_in %>%
    filter(mod_count == 0| mod_count == 1, !is.na(Res)) %>%
    group_by(MasterProteinAccessions, Sequence, Res) %>%
    summarize(N = n())  # Count the occurrences


  df_out$N <- df_out %>%
    left_join(N_df, by = c("MasterProteinAccessions", "Sequence", "Res")) %>%
    pull(N)  # Extract N column
   #colnames(df_out)[12] <- "N"  # Renaming the 12th column to "N"

  sd_df <- df_in %>%
    group_by(MasterProteinAccessions, Sequence, Res) %>%
    summarize(sdprep = sd(`Precursor Abundance`))

  df_out <- df_out %>%
    left_join(sd_df, by = c("MasterProteinAccessions", "Sequence", "Res"))
  df_out$SD <- df_out$sdprep/(df_out$SampleTotalArea+df_out$ControlTotalArea)

  df_out<- df_out[complete.cases(df_out[c("Res", "EOM")]), ]

  return(df_out)
}
