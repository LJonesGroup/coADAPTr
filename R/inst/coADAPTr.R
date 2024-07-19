##Before Analyzing your data be sure to run this function to install necessary Packages
# and resolve package conflicts.
before_beginning()

###########################################################################

# Step 1-Read in required inputs
FASTA<- FASTA_file()

# Step-2 Set output directory
file_output<- output_folder()

# Step-3 Read in data from Proteome Discoverer

raw_data<- import_data()

######TESTING SELECTION FUNCTION
column_selection <- function(df) {
  refined_data <- data.frame(matrix(ncol = 0, nrow = nrow(df)))  # Initialize an empty data frame with the same number of rows as df
  col_names <- colnames(df)

  # Function to prompt for column selection
  select_column <- function(prompt) {
    cat(prompt, "\n")
    selected_col <- select.list(
      col_names,
      multiple = FALSE,
      title = prompt,
      graphics = TRUE
    )
    if (!is.null(selected_col) && selected_col != "") {
      return(selected_col)
    } else {
      stop("No column selected, exiting the function.")
    }
  }

  # Prompt and select columns
  seq_col <- select_column("Please select the column containing the peptide sequences (unannotated):")
  acc_col <- select_column("Please select the column containing the Uniprot ID or master protein accessions:")
  mod_col <- select_column("Please select the column containing the protein modifications:")
  pre_col <- select_column("Please select the column containing precursor abundance/intensities:")
  spe_col <- select_column("Please select the column containing the spectrum file IDs:")
  cond_col <- select_column("Please select the column containing the experimental conditions:")

  # Add selected columns to refined_data
  refined_data <- cbind(refined_data, df[[seq_col]], df[[acc_col]], df[[mod_col]], df[[pre_col]], df[[spe_col]], df[[cond_col]])

  # Rename columns
  colnames(refined_data) <- c("Sequence", "Master Protein Accessions", "Modifications", "Precursor Abundance", "Spectrum File", "Condition")

  return(refined_data)
}


pd_data <- column_selection(raw_data)



# Step 4 Identify which spectrum files correlate to Sample and Control
pd_data<- SampleControl(pd_data)


# Step 5 Run this if MS files were analyzed via PD
##Skip if running Test Data Set
##Can remove duplicates in excel before using coADAPTr if desired

OG_pd_data<-pd_data


pd_data<- remove_dup(pd_data)

# Raw Data Annotations and Parsing of the FASTA File-----------------------------------------------------------
###########################################################################
# Step 6 Clean and parse data from PD output file
annotate_features <- function(raw_data) {
  raw_data <- raw_data %>%
    mutate(`Master Protein Accessions` = sapply(strsplit(`Master Protein Accessions`, ";"), `[`, 1),
           UniprotID = `Master Protein Accessions`)  # Adding this line to create the UniprotID column

  # Rename Modifications column to Mods
  raw_data <- raw_data %>%
    rename(Mods = Modifications)

  # Remove specific modifications and clean up semicolons
  raw_data <- raw_data %>%
    mutate(Modifications = Mods) %>%
    mutate(Modifications = gsub(";[A-Z]\\d+\\(Carbamidomethyl\\)", ";", Modifications)) %>%
    mutate(Modifications = gsub("[A-Z]\\d+\\(Carbamidomethyl\\);", ";", Modifications)) %>%
    mutate(Modifications = gsub("[A-Z]\\d+\\(Carbamidomethyl\\)", "", Modifications)) %>%
    mutate(Modifications = gsub(";[A-Z]\\d+\\(TMT[^)]*\\)", ";", Modifications)) %>%
    mutate(Modifications = gsub("[A-Z]\\d+\\(TMT[^)]*\\);", ";", Modifications)) %>%
    mutate(Modifications = gsub("[A-Z]\\d+\\(TMT[^)]*\\)", "", Modifications)) %>%
    mutate(Modifications = gsub("N-Term\\(Prot\\)\\(TMTpro\\);", "", Modifications)) %>%
    mutate(Modifications = gsub(";{2,}", ";", Modifications)) %>%  # Remove any double semicolons that might result from deletions
    mutate(Modifications = gsub("^;|;$", "", Modifications))      # Remove leading or trailing semicolons

  # Cleanup: Remove leading, trailing, and multiple consecutive semicolons
  raw_data <- raw_data %>%
    mutate(Modifications = gsub("\\s*;\\s*", ";", Modifications)) %>%  # Remove spaces around semicolons
    mutate(Modifications = gsub(";{2,}", ";", Modifications)) %>%  # Replace multiple semicolons with a single one
    mutate(Modifications = gsub("^;|;$", "", Modifications))  # Remove leading and trailing semicolons

  # Determine if the sequence is oxidized or unoxidized
  raw_data <- raw_data %>%
    mutate(MOD = ifelse(is.na(Modifications) | Modifications == "" | Modifications == "NA" | Modifications == " ", "Unoxidized", "Oxidized"))

  # Double check the MOD column for accuracy
  raw_data <- raw_data %>%
    mutate(MOD = ifelse(nchar(gsub("\\s+", "", Modifications)) >= 2, "Oxidized", MOD))

  # Ensure there are no empty values in the MOD column
  raw_data <- raw_data %>%
    mutate(MOD = ifelse(MOD == "" | is.na(MOD), "Unoxidized", MOD))

  # Create ModPositionL and ModPositionN columns
  raw_data <- raw_data %>%
    mutate(ModPositionL = sub("^\\s*([A-Z]).*", "\\1", Modifications)) %>%
    mutate(ModPositionN = gsub(".*?([0-9]+).*", "\\1", Modifications)) %>%
    mutate(ModPositionN = ifelse(ModPositionN == Modifications, NA, ModPositionN)) %>%
    mutate(ModPosition = ifelse(is.na(ModPositionL) | ModPositionL == "", NA,
                                paste(ModPositionL, ModPositionN, sep = "")))

  return(raw_data)
}

pd_data_annotated <- annotate_features(pd_data)

# Step 7 Parse the FASTA file for later manipulations
FASTA<- parse_fasta(FASTA)

# Step 8 Locate the residue number for the peptide termini and residues
locate_startend_res <- function(raw_data, FASTA){

  # Merge raw_data with FASTA by UniprotID
  raw_data <- merge(raw_data, FASTA, by = "UniprotID", all.x = TRUE)

  # Locate the start and end of the sequence within the protein sequence
  index <- str_locate(raw_data$protein_sequence, raw_data$Sequence)
  raw_data <- cbind(raw_data, index)

  # Create peptide column
  raw_data$peptide <- paste(raw_data$start, "-", raw_data$end, sep = "")

  # Count the number of modifications
  raw_data$mod_count <- str_count(raw_data$Modifications, "\\(.*?\\)")
  raw_data$mod_count <- ifelse(raw_data$MOD == "Unoxidized", 0, raw_data$mod_count)

  # Convert ModPositionN and start to numeric
  raw_data$ModPositionN <- as.numeric(raw_data$ModPositionN)
  raw_data$start <- as.numeric(raw_data$start)

  # Calculate modified residue positions
  raw_data$mod_res <- ifelse(!is.na(raw_data$ModPositionN) & raw_data$ModPositionN > 0, raw_data$start + raw_data$ModPositionN - 1, NA)

  # Create the Res column
  raw_data$Res <- paste(raw_data$ModPositionL, raw_data$mod_res, sep = "")

  return(raw_data)
}

pd_data_fasta_merged<- locate_startend_res(pd_data_annotated, FASTA)

# FPOP Calculations ---------------------------------------------------------------------------------------------
# Step 9 Calculating the total peptide areas and the extent of modification at the peptide level

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

  # Replace NA values with 0 where applicable
  #df_out$Control_OxidizedArea <- ifelse(df_out$Control_UnoxidizedArea > 0 & is.na(df_out$Control_OxidizedArea), 0, df_out$Control_OxidizedArea)
  #df_out$Sample_UnoxidizedArea <- ifelse(df_out$Sample_OxidizedArea > 0 & is.na(df_out$Sample_UnoxidizedArea), 0, df_out$Sample_UnoxidizedArea)

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

  #Calculating Total Valiance
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


Areas_pep<- area_calculations_pepp(pd_data_fasta_merged)




# Step 10 Subset sequence metadata like residue start/stop
##(grab_seq_metadata_pep SKIP-Used in graphing_df_pep)

# Step 11 Merge metadata with numeric graphing data

graphing_df_pep<- merge_metadata_pep(Areas_pep, pd_data_fasta_merged)


# Step 12 Filter Graphical Data to include data that is acceptable.
quant_graph_df_pep<- filtered_graphing_df_pep(graphing_df_pep)



#Residue Level-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Step 13 Calculating the total peptide areas and the extent of modification at the residue level

area_calculations_resss <- function(df_in) {

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
  sd_df <- df_in %>%
    group_by(MasterProteinAccessions, Sequence, Res) %>%
    summarize(sdprep = sd(`Precursor Abundance`, na.rm = TRUE), .groups = "drop")

  df_out <- df_out %>%
    left_join(sd_df, by = c("MasterProteinAccessions", "Sequence", "Res")) %>%
    mutate(SD = sdprep / (SampleTotalArea + ControlTotalArea + sdprep))

  # Filter out rows with missing Res or EOM values
  df_out <- df_out[complete.cases(df_out[c("Res", "EOM")]), ]

  return(df_out)
}

testresiduefunction<- area_calculations_resss(pd_data_fasta_merged)

#####TESTING RES LEVEL SD


#Calculate Variance and Standard Variation
# Group by the necessary columns including 'Condition'
RESSDTEST <- function(df_in) {
  df_out<- df_in %>%
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
  df_out <- df_out %>%
  rename(
    Control_OxidizedVar = Control,
    Sample_OxidizedVar = Sample,

    )


  # Return the final dataframe
  return(df_out)
}
sdtest<- RESSDTEST(pd_data_fasta_merged)

###UNOXIDIZED TEST
RESSDTESTall <- function(df_in) {
  df_out <- df_in %>%
    filter(mod_count == 0 | mod_count == 1) %>%
    group_by(MasterProteinAccessions, Sequence, SampleControl, MOD, Condition) %>%
    summarize(TotalArea = var(`Precursor Abundance`, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      id_cols = c("MasterProteinAccessions", "Sequence", "Condition"),
      names_from = c("SampleControl", "MOD"),
      values_from = "TotalArea",
      values_fill = NA
    )
  # Rename columns
  #df_out <- df_out %>%
    #rename(
     # Control_OxidizedVar = Control,
      #Sample_OxidizedVar = Sample,

    #)


  # Return the final dataframe
  return(df_out)
}
sdtestunox_all<- RESSDTEST(pd_data_fasta_merged)


#TESTING RES SD
RESSDTESTTOTALS <- function(df_in) {
  # Calculate Total Variance
  df_out <- df_in %>%
    filter((mod_count == 0 | mod_count == 1)) %>%
    group_by(MasterProteinAccessions, Sequence, Res, SampleControl, Condition) %>%
    summarize(TotalVar = var(`Precursor Abundance`, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      id_cols = c("MasterProteinAccessions", "Sequence", "Res", "Condition"),
      names_from = c("SampleControl"),
      values_from = "TotalVar",
      values_fill = list(TotalVar = NA)
    ) %>%
    rename(
      Sample_TotalVar = Sample,
      Control_TotalVar = Control
    )

  return(df_out)
}

# Assuming pd_data_fasta_merged is your dataset
sdtesttotals <- RESSDTESTTOTALS(pd_data_fasta_merged)










# Assuming pd_data_fasta_merged is your dataset
residue_data <- area_calculations_res(pd_data_fasta_merged)
# Assuming pd_data_fasta_merged is your dataset
residue_data <- area_calculations_res(pd_data_fasta_merged)

# Step 14  Subset sequence metadata like residue start/stop for residue level data
graphing_df_res<- graphing_data_res(Areas_res, pd_data_fasta_merged)

#Step 15 Filters for residue level graphing.

quant_graph_df_res<- filtered_graphing_df_res(graphing_df_res)



########Saving Tables and Plots------------------------------------------------------------------------------------------------------------


#Step 16A creating a table of totals
TotalsTable<-create_totals_tablelist(graphing_df_pep, graphing_df_res)


####Saving tables---------------------------------------------------------------------------

#Step 16B Save Data Frames as Excel Files

save_data_frames(file_output, TotalsTable = TotalsTable, quant_graph_df_pep = quant_graph_df_pep, quant_graph_df_res = quant_graph_df_res, graphing_df_pep = graphing_df_pep, graphing_df_res=graphing_df_res)

#### Saving Plots -----------------------------------------------------------------------

# Step 16C Saving Peptide Level Bar Graphs

generate_eom_plot_pep(df_in = quant_graph_df_pep, file_output = file_output)

#Step 16D Saving Residue Level Bar Graphs

generate_eom_plot_res(df_in = quant_graph_df_res, file_output = file_output)

#Step 17 Saving Venn Diagrams
venn_diagram()

#Step 18 Saving Peptide Level Grouped Bar Graphs

generate_grouped_bar_plot_pep()

generate_grouped_bar_plot_res()

dev.off()
