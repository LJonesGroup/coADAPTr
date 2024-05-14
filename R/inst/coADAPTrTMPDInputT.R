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

  raw_data$Modifications <- gsub(";C\\d+\\(Carbamidomethyl\\)", ";", raw_data$Modifications)
  raw_data$Modifications <- gsub("C\\d+\\(Carbamidomethyl\\);", ";", raw_data$Modifications)
  raw_data$Modifications <- gsub("C\\d+\\(Carbamidomethyl\\)", "", raw_data$Modifications)
  raw_data$Modifications <- gsub(";(?=[^A-Za-z0-9]*$)", "", raw_data$Modifications, perl = TRUE)
  raw_data$Modifications <- sub("^;\\s", "", raw_data$Modifications)
  raw_data$Modifications <- ifelse(grepl("^; [A-Z]", raw_data$Modifications),
                                   sub("^;\\s", "", raw_data$Modifications),
                                   raw_data$Modifications)
  raw_data$MOD <- ifelse(is.na(raw_data$Modifications) | raw_data$Modifications == "", "Unoxidized",
                         ifelse(grepl("()", raw_data$Modifications), "Oxidized", "Unoxidized"))
  raw_data$MOD <- ifelse(is.na(raw_data$Modifications) | raw_data$Modifications == "" | grepl("^\\s*$", raw_data$Modifications) | !grepl("[A-Za-z]", raw_data$Modifications), "Unoxidized", raw_data$MOD)

  raw_data$ModPositionL <- sub("^([[:alpha:]]*).*", "\\1",
                               raw_data$Modifications)
  raw_data$ModPositionN <- as.numeric(gsub(".*?([0-9]+).*", "\\1",
                                           raw_data$Modifications))

  raw_data$ModPosition <- ifelse(raw_data$ModPositionL == "NA", "NA",
                                 paste(raw_data$ModPositionL,
                                       raw_data$ModPositionN))
  raw_data$ModPosition <- gsub("^; ", "", raw_data$ModPosition)

  # Extract letters before ":" in the "Spectrum File" column and add to "Condition" column
  raw_data$Condition <- sub("^(.*):.*", "\\1", raw_data$`Spectrum File`)

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
locate_startend_res <- function(raw_data){

  uniqueMPA <- unique(raw_data[, c("Master Protein Accessions",'Sequence')])
  uniqueMPA <- as.data.frame(uniqueMPA)

  raw_data <- merge(raw_data, FASTA, by = "UniprotID")

  index <- str_locate(raw_data$protein_sequence, raw_data$Sequence)
  raw_data <- cbind(raw_data, index)  # TODO: memory efficiency?

  raw_data$peptide<- paste(raw_data$start,"-", raw_data$end)

  raw_data$mod_count <- str_count(raw_data$Modifications, "\\(.*?\\)")
  raw_data$mod_count <- ifelse(raw_data$MOD == "Unoxidized", 0, raw_data$mod_count)

  raw_data$mod_res <- ifelse(raw_data$ModPositionN>0, raw_data$start + raw_data$ModPositionN - 1, NA)

  raw_data$Res<- paste(raw_data$ModPositionL, raw_data$mod_res)

  return(raw_data)
}

pd_data_fasta_merged <- locate_startend_res(annotated_data)

#Perform EOM Calculations

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
