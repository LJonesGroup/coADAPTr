before_beginning <- function() {
  options(stringsAsFactors = FALSE)

  required_packages <- c("ExcelFunctionsR", "remotes", "data.table", "stringr", "protr",
                         "plyr", "extrafont", "readxl", "ggplot2", "eulerr",
                         "tidyverse", "EnvStats", "dplyr", "writexl", "conflicted",
                         "phylotools", "parallel", "rlist","argparser", "Cairo", "VennDiagram")


  package.check <- lapply(
    required_packages,
    FUN = function(x) {
      if (!require(x, character.only = TRUE)) {
        install.packages(x, dependencies = TRUE)
        library(x, character.only = TRUE)
      }
    }
  )


  setwd(getwd())


  conflict_prefer("arrange", winner = "dplyr")
  conflict_prefer("filter", winner = "dplyr")
  conflict_prefer("between", winner = "dplyr")
  conflict_prefer("compact", winner = "purrr")
  conflict_prefer("count", winner = "dplyr")
  conflict_prefer("desc", winner = "dplyr")
  conflict_prefer("failwith", winner = "dplyr")
  conflict_prefer("filter", winner = "dplyr")
  conflict_prefer("first", winner = "dplyr")
  conflict_prefer("hour", winner = "lubridate")
  conflict_prefer("id", winner = "dplyr")
  conflict_prefer("isoweek", winner = "lubridate")
  conflict_prefer("lag", winner = "dplyr")
  conflict_prefer("last", winner = "dplyr")
  conflict_prefer("mday", winner = "lubridate")
  conflict_prefer("minute", winner = "lubridate")
  conflict_prefer("month", winner = "lubridate")
  conflict_prefer("mutate", winner = "dplyr")
  conflict_prefer("quarter", winner = "lubridate")
  conflict_prefer("rename", winner = "dplyr")
  conflict_prefer("second", winner = "lubridate")
  conflict_prefer("summarise", winner = "dplyr")
  conflict_prefer("summarize", winner = "dplyr")
  conflict_prefer("transpose", winner = "purrr")
  conflict_prefer("wday", winner = "lubridate")
  conflict_prefer("week", winner = "lubridate")
  conflict_prefer("yday", winner = "lubridate")
  conflict_prefer("year", winner = "lubridate")
  conflict_prefer("read.fasta", winner = "phylotools")
  conflict_prefer("rename", winner = "dplyr")
}

before_beginning()
###########################################################################
#Homo sapien Reviewed 12062021
# Read in required inputs
FASTA_file <- function() {
  FASTA_path <- file.choose()
  FASTA <- read.fasta(FASTA_path)
  return(FASTA)
}

FASTA<-FASTA_file()

#USER INDICATES FILE PATH
import_data <- function() {
  library(openxlsx)

  file_path <- file.choose()
  df <- read.xlsx(file_path, check.names = FALSE)

  colnames(df) <- gsub("\\.", " ", colnames(df))

  return(df)
}

pd_data<-import_data()


# Set output directory
output_folder <- function() {
  file_output <- choose.dir()
  return(file_output)
}

file_output <- output_folder()

#Calculate Total Abundance (denominator)
sum_abundance <- function(data) {
  # Find columns containing "Abundance"
  abundance_cols <- grep("Abundance", colnames(data), value = TRUE)

  # Sum numerical content of abundance columns for each row
  data$AbundanceTotals <- rowSums(data[, abundance_cols, drop = FALSE], na.rm = TRUE)

  return(data)
}

# Apply function
pd_data <- sum_abundance(pd_data)

#Calculate True Precursor Abundance
calculate_abundances <- function(data) {
  # Find columns containing "Abundance"
  abundance_cols <- grep("Abundance:", colnames(data), value = TRUE)

  # Loop through each abundance column
  for (col in abundance_cols) {
    # Extract the number and letter from the original column name
    num_letter <- gsub("Abundance: ", "", col)

    # Calculate the new column
    new_col <- paste("TRUE Abundance:", num_letter, sep = " ")
    data[[new_col]] <- data[[col]] * data$Intensity / data$AbundanceTotals
  }

  return(data)
}

pd_data<- calculate_abundances(pd_data)

#Extract Necessary Columns
extract_columns <- function(data) {
  # Select the desired columns
  selected_cols <- c("Sequence", "Modifications", "Master Protein Accessions", "Protein Accessions")

  # Identify indices of columns containing "TRUE Abundance"
  abundance_indices <- grepl("TRUE Abundance", colnames(data))

  # Combine selected columns with abundance columns
  all_selected_cols <- c(selected_cols, colnames(data)[abundance_indices])

  # Subset the dataframe
  data_subset <- data[, all_selected_cols, drop = FALSE]

  return(data_subset)
}



cleaned_data<- extract_columns(pd_data)

#renaming columns
# Function to interactively rename columns
rename_columns_interactively <- function(data) {
  # Prompt user to select columns for renaming
  selected_columns <- select.list(
    colnames(data),
    multiple = TRUE,
    title = "Select columns to rename:",
    graphics = TRUE
  )

  # Loop through selected columns and rename them
  for (col in selected_columns) {
    new_name <- readline(prompt = paste("Enter new name for column", col, "in the format Condition:SampleType: "))
    # Rename the column
    colnames(data)[colnames(data) == col] <- new_name
  }

  return(data)
}


renamed_data <- rename_columns_interactively(cleaned_data)

# Transform the data in the column into rows


transform_data <- function(data) {
  # Pivot the data into long format, splitting columns containing ":"
  transformed_data <- pivot_longer(data, cols = contains(":"), names_to = "Spectrum File", values_to = "Precursor Abundance")

  return(transformed_data)
}

transformed_data <- transform_data(renamed_data)

#Identify Sample and Control Files
SampleControl <- function(pd_data) {
  pd_data$SampleControl <- ifelse(grepl("NL", pd_data$`Spectrum File`), "Control", "Sample")
  return(pd_data)
}
transformed_data <- SampleControl(transformed_data)

#Annotate features in transformed data
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
    mutate(MOD = ifelse(is.na(Modifications) | Modifications == "", "Unoxidized",
                        ifelse(grepl("Oxidation", Modifications), "Oxidized", "Unoxidized")))

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



annotated_data<-annotate_features(transformed_data)

#Parse the FASTA file

parse_fasta <- function(fasta_in){
  fasta_in <- fasta_in %>%
    rename(protein_sequence = seq.text)
  fasta_in<- fasta_in %>%
    rename(MasterProteinAccessions = seq.name)

  fasta_in$UniprotID <- gsub("^.+\\|(\\w+)\\|.*$", "\\1", fasta_in$MasterProteinAccessions)

  return(fasta_in)

}
FASTA<- parse_fasta(FASTA)

#Merge the FASTA file with the annotated data
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


pd_data_fasta_merged <- locate_startend_res(annotated_data, FASTA)

#Perform EOM Calculations at the peptide level

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
  df_out$Control_OxidizedArea <- ifelse(df_out$Control_UnoxidizedArea > 0 & is.na(df_out$Control_OxidizedArea), 0, df_out$Control_OxidizedArea)
  df_out$Sample_UnoxidizedArea <- ifelse(df_out$Sample_OxidizedArea > 0 & is.na(df_out$Sample_UnoxidizedArea), 0, df_out$Sample_UnoxidizedArea)

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

  # Calculate standard deviation, replacing NA with 0
  test <- df_in %>%
    group_by(MasterProteinAccessions, Sequence, Condition) %>%
    summarize(sdprep = sd(replace_na(`Precursor Abundance`, 0), na.rm = TRUE))
  df_out <- df_out %>%
    left_join(test, by = c("MasterProteinAccessions", "Sequence", "Condition")) %>%
    mutate(SD = sdprep / (TotalSampleArea + TotalControlArea + sdprep))

  # Return the final dataframe
  return(df_out)
}




Areas_pep<- area_calculations_pep(pd_data_fasta_merged)

#Merge metadata
grab_seq_metadata_pep <- function(df_in){
  df_out <- subset(df_in, select = c("MasterProteinAccessions", "Sequence",
                                     "peptide", "start"))
  df_out <- df_out[!duplicated(df_out), ]
  return(df_out)
}

merge_metadata_pep <- function(df_in, rawdatafastamerged) {
  df_out <- left_join(df_in, grab_seq_metadata_pep(rawdatafastamerged))

  df_out <- df_out[order(df_out$start), ]

  df_out$MasterProteinAccessions <- gsub(".*\\|(.*?)\\|.*", "\\1", df_out$MasterProteinAccessions)

  return(df_out)
}
graphing_df_pep<- merge_metadata_pep(Areas_pep, pd_data_fasta_merged)

#Filter out quantifiable modifications
filtered_graphing_df_pep <- function(df_in) {
  df_out = df_in[df_in$EOM > 0 & df_in$EOM > df_in$SD & df_in$N > 4, ]
  df_out<- df_out %>%
    arrange(start)
  df_out <- df_out[!is.na(df_out$MasterProteinAccessions), ]
  return(df_out)
}
quant_graph_df_pep<- filtered_graphing_df_pep(graphing_df_pep)

#Perform EOM Calculations at the residue level

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

Areas_res<-area_calculations_res(pd_data_fasta_merged)

#Merge metadata for residue level quantification
grab_seq_metadata_res <- function(df_in){
  df_out <- subset(df_in, select = c("MasterProteinAccessions", "Sequence", "Res", "start"))
  df_out <- df_out[!duplicated(df_out), ]
  return(df_out)
}
graphing_data_res <- function(df_in, pd_data_fasta_merged) {
  df_out <- df_in %>%
    left_join(grab_seq_metadata_res(pd_data_fasta_merged)) %>%
    filter(!(is.na(Res) | Res == "")) %>%
    arrange(start)%>%
    mutate(MasterProteinAccessions = gsub(".*\\|(.*?)\\|.*", "\\1", MasterProteinAccessions))

  return(df_out)
}

graphing_df_res<- graphing_data_res(Areas_res, pd_data_fasta_merged)

#Filter out quantifiable modifications at the residue level
filtered_graphing_df_res <- function(df_in) {
  df_out = df_in[df_in$EOM > 0 & df_in$EOM > df_in$SD & df_in$N > 4, ]
  df_out<- df_out %>%
    arrange(start)
  df_out <- df_out[!is.na(df_out$MasterProteinAccessions), ]
  return(df_out)
}

quant_graph_df_res<- filtered_graphing_df_res(graphing_df_res)

#Createing Totals Table for saving

create_totals_tablelist <- function(df_in, df_res) {
  all_conditions <- unique(df_in$Condition)
  result_list <- list()

  for (condition in all_conditions) {
    df_condition_in <- df_in %>%
      filter(Condition == condition)

    df_condition_res <- df_res %>%
      filter(Condition == condition)

    df_out <- data.frame(
      Condition = condition,
      UniqueProteinDet = NA,
      UniqueSeqDet = NA,
      UniqueResDet = NA,
      QuantifiableModProtein = NA,
      QuantifiableModSeq = NA,
      QuantifiableModRes = NA
    )

    unique_protein_det <- unique(df_condition_in$MasterProteinAccessions)
    df_out$UniqueProteinDet <- length(unique_protein_det)

    unique_seq_det <- unique(df_condition_in$Sequence)
    df_out$UniqueSeqDet <- length(unique_seq_det)

    unique_seq_res_det <- unique(df_condition_res %>%
                                   select(Sequence, Res))
    df_out$UniqueResDet <- nrow(unique_seq_res_det)

    df_filtered <- df_condition_in %>%
      filter(EOM > 0 & EOM > SD & !is.na(SD))

    df_filtered_res <- df_condition_res %>%
      filter(EOM > 0 & EOM > SD)

    quantifiable_protein_mod <- unique(df_filtered$MasterProteinAccessions)
    df_out$QuantifiableModProtein <- length(quantifiable_protein_mod)

    quantifiable_seq_mod <- unique(df_filtered$Sequence)
    df_out$QuantifiableModSeq <- length(quantifiable_seq_mod)

    quantifiable_seq_res_mod <- unique(df_filtered_res %>%
                                         select(Sequence, Res))
    df_out$QuantifiableModRes <- nrow(quantifiable_seq_res_mod)

    result_list[[condition]] <- df_out
  }

  final_df <- bind_rows(result_list)

  return(final_df)
}



TotalsTable<-create_totals_tablelist(graphing_df_pep, graphing_df_res)


####Saving tables---------------------------------------------------------------------------

#Step 18 Save Data Frames as Excel Files
save_data_frames <- function(output_directory, ...) {
  # Prompt the user to input the file name
  cat("Enter the file name (without extension): ")
  file_name <- readline()

  # Remove any leading or trailing whitespace
  file_name <- trimws(file_name)

  # Check if the file name is empty, if so, use a default name
  if (nchar(file_name) == 0) {
    file_name <- "output"
  }

  # Create a list of data frames
  data_frames <- list(...)

  # Iterate over each data frame and save as separate Excel files
  for (i in seq_along(data_frames)) {
    # Get the current data frame
    df <- data_frames[[i]]

    # Create the output file name
    output_file_name <- paste0(file_name, "_", names(data_frames)[i])

    # Create the output file path in the output directory
    output_file_path <- file.path(output_directory, paste0(output_file_name, ".xlsx"))

    # Save the data frame as an Excel file
    writexl::write_xlsx(df, output_file_path)

    # Print a message to indicate successful saving
    cat(names(data_frames)[i], "saved as", output_file_path, "\n")
  }
}

save_data_frames(file_output, TotalsTable = TotalsTable, Areas_pep = Areas_pep, Areas_res = Areas_res, quant_graph_df_pep = quant_graph_df_pep, quant_graph_df_res = quant_graph_df_res, graphing_df_pep = graphing_df_pep, graphing_df_res=graphing_df_res)

#### Saving Plots -----------------------------------------------------------------------
# Creating a bar graph for each condition at the peptide level
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

#####
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
