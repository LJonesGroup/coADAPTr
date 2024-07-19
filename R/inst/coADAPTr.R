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

rename_and_split_spectrum_files <- function(df_in) {
  # Check if the "Spectrum File" column exists
  if ("Spectrum File" %in% colnames(df_in)) {
    # Extract the part before the first '.' to find unique files
    df_in$FileIdentifier <- sapply(strsplit(as.character(df_in$`Spectrum File`), "\\."), `[`, 1)
    unique_files <- unique(df_in$FileIdentifier)
    cat("Unique File Identifiers detected:\n")
    print(unique_files)

    # Initialize a data frame to hold the mappings
    mappings <- data.frame(Original = character(), SampleType = character(), Condition = character(), stringsAsFactors = FALSE)

    # Loop through each unique file identifier
    for (file in unique_files) {
      response <- readline(prompt = paste("Enter new name for '", file, "' in the format 'Condition:SampleType': ", sep=""))
      parts <- strsplit(response, ":")[[1]]
      if (length(parts) == 2) {
        # Append the mappings
        mappings <- rbind(mappings, data.frame(Original = file, SampleType = parts[2], Condition = parts[1]))
      } else {
        cat("Invalid input. Skipping '", file, "'\n")
      }
    }

    # Rename and assign based on mappings
    if (nrow(mappings) > 0) {
      # Map FileIdentifier to SampleType
      df_in$SampleType <- df_in$FileIdentifier
      for (i in 1:nrow(mappings)) {
        df_in$SampleType[df_in$FileIdentifier == mappings$Original[i]] <- mappings$SampleType[i]
        df_in$Condition[df_in$FileIdentifier == mappings$Original[i]] <- mappings$Condition[i]
      }
      # Optionally remove the temporary identifier column
      df_in$FileIdentifier <- NULL
      # Rename column after all mappings are applied
      colnames(df_in)[colnames(df_in) == "Spectrum File"] <- "SampleType"
    }
  } else {
    cat("'Spectrum File' column not found in the data frame.\n")
  }

  return(df_in)
}

# Usage:
# cleaned_data should be your data frame variable
modified_data <- rename_and_split_spectrum_files(pd_data)



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


Areas_pep<- area_calculations_pep(pd_data_fasta_merged)




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

Areas_res<- area_calculations_resss(pd_data_fasta_merged)

###
##CORRECT SD






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

generate_grouped_bar_plot_pep <- function() {
  # Prompt the user to select the Excel file
  cat("Select the Excel file containing the data: ")
  filepath <- file.choose()

  # Load the data from the selected Excel file
  df_in <- readxl::read_excel(filepath)

  # Arrange the dataframe by start
  df_in <- df_in %>%
    arrange(start)

  # Auto-detect unique conditions
  unique_conditions <- unique(df_in$Condition)

  # Manually specify blue colors for filling the bars
  condition_colors <- c("#0570b0", "#3690c0", "#74a9cf", "#a6bddb", "#d0d1e6", "#084594", "#2171b5", "#4292c6", "#6baed6", "#9ecae1", "#c6dbef", "#08519c", "#3182bd", "#6baed6", "#9ecae1", "#c6dbef")

  # Create a factor variable to represent the sorted order of conditions
  df_in$Condition <- factor(df_in$Condition, levels = unique_conditions)

  # Prompt the user to specify the filename
  cat("Enter the filename (without extension): ")
  filename <- readline()

  # Remove any leading or trailing whitespace
  filename <- trimws(filename)

  # Iterate over each protein and make a grouped bar plot for it
  for (protein in unique(df_in$MasterProteinAccessions)) {
    # Subset the dataframe for this protein
    temp <- subset(df_in, MasterProteinAccessions == protein)

    # Generate a grouped bar plot of the extent of modification for each peptide
    # that maps to this protein, with different conditions represented by color
    fig <- ggplot(temp, aes(x = peptide, y = EOM, fill = Condition)) +
      geom_bar(position = position_dodge(width = 0.9), stat = "identity") +
      geom_errorbar(aes(ymin = EOM - SD, ymax = EOM + SD), position = position_dodge(width = 0.9), width = 0.25) +
      labs(title = paste(protein, "Peptide Level Analysis"),
           x = "Peptide",
           y = "Extent of Modification",
           fill = "Condition") +
      theme_classic() +
      theme(text = element_text(size = 18, family = "sans"),
            plot.title = element_text(hjust = 0.5, size = 20, family = "sans", face = "bold"),
            legend.text = element_text(size = 16, family = "sans"),
            legend.title = element_text(size = 18, family = "sans"),
            axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
            axis.title.x = element_text(margin = margin(t = 20))) +
      scale_fill_manual(values = condition_colors)

    # Create the output directory for bar graphs based on the file output and excel filename
    graph_output_directory <- file.path(file_output, paste0(filename, "_PeptideLevelGroupedBarGraphs"))
    dir.create(graph_output_directory, showWarnings = FALSE, recursive = TRUE)

    # Generate the full file path for this protein and save the figure
    file_out <- file.path(graph_output_directory, paste0(protein, ".png"))
    ggsave(filename = file_out, plot = fig, device = "png")

    # Print a message to indicate successful saving
    cat("Grouped bar graph for", protein, "saved as", file_out, "\n")
  }
}

generate_grouped_bar_plot_pep()




#Step 16D Saving Residue Level Bar Graphs


# Generating grouped bar plot at the residue level
generate_grouped_bar_plot_res <- function() {
  # Prompt the user to select the Excel file
  cat("Select the Excel file containing the data: ")
  filepath <- file.choose()

  # Load the data from the selected Excel file
  df_in <- readxl::read_excel(filepath)

  # Arrange the dataframe by start
  df_in <- df_in %>%
    arrange(start)

  # Auto-detect unique conditions
  unique_conditions <- unique(df_in$Condition)

  # Manually specify blue colors for filling the bars
  condition_colors <- c("#0570b0", "#3690c0", "#74a9cf", "#a6bddb", "#d0d1e6", "#084594", "#2171b5", "#4292c6", "#6baed6", "#9ecae1", "#c6dbef", "#08519c", "#3182bd", "#6baed6", "#9ecae1", "#c6dbef")

  # Create a factor variable to represent the sorted order of conditions
  df_in$Condition <- factor(df_in$Condition, levels = unique_conditions)

  # Prompt the user to specify the filename
  cat("Enter the filename (without extension): ")
  filename <- readline()

  # Remove any leading or trailing whitespace
  filename <- trimws(filename)

  # Iterate over each protein and make a grouped bar plot for it
  for (protein in unique(df_in$MasterProteinAccessions)) {
    # Subset the dataframe for this protein
    temp <- subset(df_in, MasterProteinAccessions == protein)

    # Generate a grouped bar plot of the extent of modification for each peptide
    # that maps to this protein, with different conditions represented by color
    fig <- ggplot(temp, aes(x = Res, y = EOM, fill = Condition)) +
      geom_bar(position = position_dodge(width = 0.9), stat = "identity") +
      geom_errorbar(aes(ymin = EOM - SD, ymax = EOM + SD), position = position_dodge(width = 0.9), width = 0.25) +
      labs(title = paste(protein, "Residue Level Analysis"),
           x = "Residue",
           y = "Extent of Modification",
           fill = "Condition") +
      theme_classic() +
      theme(text = element_text(size = 18, family = "sans"),
            plot.title = element_text(hjust = 0.5, size = 20, family = "sans", face = "bold"),
            legend.text = element_text(size = 16, family = "sans"),
            legend.title = element_text(size = 18, family = "sans"),
            axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
            axis.title.x = element_text(margin = margin(t = 20))) +
      scale_fill_manual(values = condition_colors)

    # Create the output directory for bar graphs based on the file output and excel filename
    graph_output_directory <- file.path(file_output, paste0(filename, "_ResidueLevelGroupedBarGraphs"))
    dir.create(graph_output_directory, showWarnings = FALSE, recursive = TRUE)

    # Generate the full file path for this protein and save the figure
    file_out <- file.path(graph_output_directory, paste0(protein, ".png"))
    ggsave(filename = file_out, plot = fig, device = "png")

    # Print a message to indicate successful saving
    cat("Grouped bar graph for", protein, "saved as", file_out, "\n")
  }
}

generate_grouped_bar_plot_res()

#Step 17 Saving Venn Diagrams
#Generate Venn Diagrams
venn_diagram <- function() {
  # Prompt user to select Excel file
  excel_file <- file.choose()

  # Read data from Excel file
  data <- read_excel(excel_file, col_names = TRUE)

  # Prompt user to select output folder
  output_folder <- choose.dir()

  # Get the number of conditions
  num_conditions <- ncol(data)

  # Get bubble names from the user
  venn_bubble_names <- readline(prompt = "Enter names for the Venn diagram bubbles separated by commas: ")
  venn_bubble_names <- strsplit(venn_bubble_names, ",")[[1]]

  # Ensure the number of bubble names matches the number of conditions
  if (length(venn_bubble_names) != num_conditions) {
    stop("Number of bubble names provided does not match the number of columns in the Excel file.")
  }

  # Process each column separately
  condition_lists <- lapply(seq_len(num_conditions), function(i) {
    column_data <- na.omit(data[[i]])
    if (length(column_data) == 0) {
      return(NULL)
    } else {
      return(column_data)
    }
  })

  # Generate a custom color palette (using shades of blue or gray)
  custom_palette <- c("#0570b0", "#3690c0", "#74a9cf", "#a6bddb", "#d0d1e6", "#084594", "#2171b5", "#4292c6", "#6baed6", "#9ecae1", "#c6dbef", "#08519c", "#3182bd", "#6baed6", "#9ecae1", "#c6dbef")
  custom_palette <- custom_palette[1:num_conditions] # ensuring the palette has the desired number of colors

  # Create Venn diagram plot
  venn_plot <- venn.diagram(
    x = condition_lists,
    category.names = venn_bubble_names,
    filename = NULL,
    col = custom_palette,
    fill = custom_palette,
    alpha = 0.5, # Adjust transparency for better visualization
    margin = 0.1, # Increase margin for better plot appearance
    fontfamily = "sans", # Set font family to Helvetica (or similar font)
    fontface = "bold",   # Set font face to bold
    cat.fontsize = 20,   # Set category font size to 20
    cex = 2.5,           # Adjust overall font size
    cat.cex = 2.5,         # Set category title font size to 24
    cat.fontfamily = "sans", # Set category title font family
    cat.fontface = "bold"  # Set category title font face to bold
  )

  # Prompt user to specify the name of the output PNG file
  png_file_name <- readline(prompt = "Enter the name of the output PNG file (without extension): ")
  png_file <- file.path(output_folder, paste0(png_file_name, ".png"))

  # Save plot as PNG
  png(filename = png_file, width = 1000, height = 1000) # Adjust width and height as needed
  grid.draw(venn_plot)
  dev.off()

  cat("Venn diagram plot saved at:", png_file, "\n")

  # Calculate overlap
  overlap <- calculate.overlap(condition_lists)

  # Save overlap to Excel
  overlap_file <- file.path(output_folder, paste0(png_file_name, "_Overlap.xlsx"))
  write.xlsx(overlap, overlap_file)

  cat("Overlap information saved at:", overlap_file, "\n")

  if (FALSE) {
    venn_diagram()
  }
}
venn_diagram()


#Step 18 Saving Modifications Per Peptide and Residue

#Counting modified peptides per protein and graphing

count_peptides_per_protein <- function() {
  data_list <- list()
  conditions <- c()

  repeat {
    file_path <- file.choose()
    condition <- readline(prompt = "What condition does this data correspond to? ")
    conditions <- c(conditions, condition)

    data <- read_excel(file_path)

    count_data <- data %>%
      group_by(MasterProteinAccessions) %>%
      summarise(SequenceCount = n()) %>%
      ungroup()

    data_list[[condition]] <- count_data

    more_files <- readline(prompt = "Do you want to input another file? (yes/no) ")
    if (tolower(more_files) != "yes") {
      break
    }
  }

  graph_title <- readline(prompt = "Enter the title for the Bar Graph: ")

  combined_data <- bind_rows(data_list, .id = "Condition")
  summary_data <- combined_data %>%
    group_by(Condition, SequenceCount) %>%
    summarise(ProteinCount = n(), .groups = 'drop')

  p <- ggplot(summary_data, aes(x = SequenceCount, y = ProteinCount, fill = Condition)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = graph_title, x = "Number of Peptides per Protein", y = "Number of Proteins") +
    theme_minimal(base_size = 15) +
    theme(panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA),
          legend.background = element_rect(fill = "white", colour = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 20))

  print(p)

  output_dir <- choose.dir(default = "", caption = "Select directory to save the count data and graph")
  file_name <- readline(prompt = "Enter the base name for the saved files (without extension): ")

  output_file <- file.path(output_dir, paste0(file_name, ".xlsx"))
  graph_file <- file.path(output_dir, paste0(file_name, ".png"))

  write_xlsx(summary_data, output_file)
  ggsave(graph_file, plot = p, width = 10, height = 6)

  print(paste("Data saved to", output_file))
  print(paste("Graph saved to", graph_file))
}

# To run the function
count_peptides_per_protein()
####


#Count Residues Per Protein
count_residue_entries_per_protein <- function() {
  data_list <- list()
  conditions <- c()

  repeat {
    file_path <- file.choose()
    condition <- readline(prompt = "What condition does this data correspond to? ")
    conditions <- c(conditions, condition)

    data <- read_excel(file_path)

    # Ensure "Res" column contains non-NA values and only valid entries with numbers followed by letters
    data <- data %>% filter(!is.na(Res) & grepl("^[A-Za-z][0-9]+", Res))

    count_data <- data %>%
      group_by(MasterProteinAccessions) %>%
      summarise(ResidueEntryCount = n()) %>%
      ungroup()

    data_list[[condition]] <- count_data

    more_files <- readline(prompt = "Do you want to input another file? (yes/no) ")
    if (tolower(more_files) != "yes") {
      break
    }
  }

  graph_title <- readline(prompt = "Enter the title for the Bar Graph: ")

  combined_data <- bind_rows(data_list, .id = "Condition")
  summary_data <- combined_data %>%
    group_by(Condition, ResidueEntryCount) %>%
    summarise(ProteinCount = n(), .groups = 'drop')

  p <- ggplot(summary_data, aes(x = ResidueEntryCount, y = ProteinCount, fill = Condition)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = graph_title, x = "Number of Residue Entries per Protein", y = "Number of Proteins") +
    theme_minimal(base_size = 15) +
    theme(panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA),
          legend.background = element_rect(fill = "white", colour = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 20))

  print(p)

  output_dir <- choose.dir(default = "", caption = "Select directory to save the count data and graph")
  file_name <- readline(prompt = "Enter the base name for the saved files (without extension): ")

  output_file <- file.path(output_dir, paste0(file_name, ".xlsx"))
  graph_file <- file.path(output_dir, paste0(file_name, ".png"))

  write_xlsx(summary_data, output_file)
  ggsave(graph_file, plot = p, width = 10, height = 6)

  print(paste("Data saved to", output_file))
  print(paste("Graph saved to", graph_file))
}

# To run the function
count_residue_entries_per_protein()


dev.off()
