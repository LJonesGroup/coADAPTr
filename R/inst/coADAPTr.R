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
  colnames(refined_data) <- c("Sequence", "Master Protein Accessions", "Modifications", "Precursor Abundance", "Spectrum File")

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
pd_data_annotated<- annotate_features(pd_data)
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

  # Correct the MOD column after cleaning up Carbamidomethyl modifications
  raw_data <- raw_data %>%
    mutate(MOD = ifelse(Mods == "" | Mods == "NA" | grepl("Carbamidomethyl", Mods), "Unoxidized", MOD))

  # Create ModPositionL and ModPositionN columns
  raw_data <- raw_data %>%
    mutate(ModPositionL = sub("^\\s*([A-Z]).*", "\\1", Modifications)) %>%
    mutate(ModPositionN = gsub(".*?([0-9]+).*", "\\1", Modifications)) %>%
    mutate(ModPositionN = ifelse(ModPositionN == Modifications, NA, ModPositionN)) %>%
    mutate(ModPosition = ifelse(is.na(ModPositionL) | ModPositionL == "", NA,
                                paste(ModPositionL, ModPositionN, sep = "")))

  # Extract letters before ":" in the "Spectrum File" column and add to "Condition" column
  raw_data <- raw_data %>%
    mutate(Condition = sub("^(.*):.*", "\\1", `Spectrum File`))

  return(raw_data)
}



pd_data_annotated<-annotate_features(pd_data)

# Step 7 Parse the FASTA file for later manipulations
FASTA<- parse_fasta(FASTA)

# Step 8 Locate the residue number for the peptide termini and residues

pd_data_fasta_merged <- locate_startend_res(pd_data_annotated)



# FPOP Calculations ---------------------------------------------------------------------------------------------
# Step 9 Calculating the total peptide areas and the extent of modification at the peptide level

Areas_pep<- area_calculations_pep(pd_data_fasta_merged)


# Step 10 Subset sequence metadata like residue start/stop
##(grab_seq_metadata_pep SKIP-Used in graphing_df_pep)

# Step 11 Merge metadata with numeric graphing data

graphing_df_pep<- merge_metadata_pep(Areas_pep, pd_data_fasta_merged)


# Step 12 Filter Graphical Data to include data that is acceptable.
quant_graph_df_pep<- filtered_graphing_df_pep(graphing_df_pep)



#Residue Level-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Step 13 Calculating the total peptide areas and the extent of modification at the residue level

Areas_res <- area_calculations_res(pd_data_fasta_merged)

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
