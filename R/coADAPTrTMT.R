options(stringsAsFactors = FALSE)

# list of required packages
required_packages <- c("ExcelFunctionsR", "remotes", "data.table", "stringr", "protr",
                       "plyr", "extrafont", "readxl", "ggplot2", "eulerr",
                       "tidyverse", "EnvStats", "dplyr", "writexl","stringr",
                       "phylotools", "parallel", "rlist","argparser", "Cairo")

# load or install & load all packages
package.check <- lapply(
  required_packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

# set working directory to current directory
setwd(getwd())


###########################################################################
#Homo sapien Reviewed 12062021
# Read in required inputs
FASTA <- read.fasta(file = ("./data/SWISSPROT Homo Sapian 9606 + cRAP Reviewed.fasta"),
                    clean_name = FALSE)
#USER INDICATES FILE PATH
file_path<- "C:/Users/Raqie/Desktop/UMBPSC/RCoding/LJonesGroup/coADAPTr/coADAPTr/data/EECTMT.xlsx"

#UMBPSC/RCoding/JonesLabScripts/data/100 mM H2O2 1 Min BR1 Core.xlsx"
#"C:/Users/Raqie/Desktop/UMBPSC/RCoding/JonesLabScripts/data/100 mM H2O2 1 Min BR1 Outer.xlsx"

# Set output directory
file_output= "C:/Users/Raqie/Desktop/AutoDataAnalysis/TMT"



#Read in data from Proteome Discoverer
pd_data <- read_excel(file_path)

##If your file input is from fragpipe start here

#extracting proetin sequence
#adding L or NL to files
#extracting modification
# Troubleshooting
# List of column names containing "Abundance"
abundance_cols <- grep("Abundance", colnames(pd_data), value = TRUE)

# Define a function to convert to numeric with NAs as 0
convert_to_numeric <- function(x) {
  ifelse(is.na(x), 0, as.numeric(x))
}

# Apply the function to abundance columns
pd_data[, abundance_cols] <- lapply(pd_data[, abundance_cols], convert_to_numeric)

# Now you can perform the calculations
# Example: Sum of Abundance: 126 and Abundance: 127N columns
pd_data$SumRow_Abundances <- rowSums(pd_data[, c("Abundance: 126", "Abundance: 127N", "Abundance: 127C", "Abundance: 128N", "Abundance: 128C", "Abundance: 129N", "Abundance: 129C", "Abundance: 130C", "Abundance: 131N", "Abundance: 131C")], na.rm = TRUE)
#Assign new column names here.
#Label which TMT channel corresponds to which sample by changing the column name (pd_data$*****)
pd_data$VehicleL1 <- (pd_data$"Abundance: 126" * pd_data$Intensity )/ pd_data$SumRow_Abundances
pd_data$VehicleL2 <- (pd_data$"Abundance: 127N" * pd_data$Intensity ) / pd_data$SumRow_Abundances
pd_data$VehicleL3 <- (pd_data$"Abundance: 127C" * pd_data$Intensity ) / pd_data$SumRow_Abundances
pd_data$VehicleNL1 <- (pd_data$"Abundance: 128N" * pd_data$Intensity ) / pd_data$SumRow_Abundances
pd_data$VehicleNL1 <- (pd_data$"Abundance: 128C" * pd_data$Intensity ) / pd_data$SumRow_Abundances
#pd_data$VehicleNL3 <- (pd_data$"Abundance: 129" * pd_data$Intensity ) / pd_data$SumRow_Abundances      #Not used in test data
pd_data$DrugL1 <- (pd_data$"Abundance: 129N" * pd_data$Intensity ) / pd_data$SumRow_Abundances
pd_data$DrugL2 <- (pd_data$"Abundance: 129C" * pd_data$Intensity ) / pd_data$SumRow_Abundances
pd_data$DrugL3 <- (pd_data$"Abundance: 130C" * pd_data$Intensity ) / pd_data$SumRow_Abundances
pd_data$DrugNL1 <- (pd_data$"Abundance: 131N" * pd_data$Intensity ) / pd_data$SumRow_Abundances
pd_data$DrugNL2 <- (pd_data$"Abundance: 131C" * pd_data$Intensity ) / pd_data$SumRow_Abundances
#pd_data$DrugnNL3 <- (pd_data$"Abundance: 130N" * pd_data$Intensity ) / pd_data$SumRow_Abundances

#"Abundance: 126"  "Abundance: 127N" "Abundance: 127C" "Abundance: 128N" "Abundance: 128C" "Abundance: 129N" "Abundance: 129C" "Abundance: 130C" "Abundance: 131N" "Abundance: 131C"

# Calculate the result for each column, handling NAs


# Rename the columns
new_colnames <- c(  "Vehicle NL2",
                  "Drug L1", "Drug L2", "Drug L3", "Drug NL1", "Drug NL2")
colnames(pd_data)[1:length(abundance_cols)] <- new_colnames

return(pd_data)






calculate_tmt_fpop <- function(pd_data) {
  # List of column names containing "Abundance"
  abundance_cols <- grep("Abundance", colnames(pd_data), value = TRUE)

  # Convert abundance columns to numeric, handling NAs as 0
  pd_data[, abundance_cols] <- lapply(pd_data[, abundance_cols], function(x) as.numeric(ifelse(x == "NA", 0, x)))

  # Calculate the result for each column, handling NAs
  for (i in 1:length(abundance_cols)) {
    pd_data[, i] <- ifelse(pd_data[, i] == 0, NA, (pd_data[, i] * pd_data$Intensity) / rowSums(pd_data[, abundance_cols]))
  }

  # Rename the columns
  new_colnames <- c("Vehicle L1", "Vehicle L2", "Vehicle L3", "Vehicle NL1", "Vehicle NL2",
                    "Drug L1", "Drug L2", "Drug L3", "Drug NL1", "Drug NL2")
  colnames(pd_data)[1:length(abundance_cols)] <- new_colnames

  return(pd_data)
}

# Call the function with your pd_data dataframe
pd_data <- calculate_tmt_fpop(pd_data)



#Does this need to have a separate function or can it be merged with Annotate Features.
# Annotate Sample/Control based on MS acquisitoin filename convention
pd_data$SampleControl <- ifelse(pd_data$'Spectrum File' %like%
                                  "NL", "Control", "Sample")

#####Run this if MS files were analyzed via PD

remove_dup <- function(df) {

  identical_cols <- setdiff(names(df), c("Identifying Node", "DeltaScore"))

  df <- df %>%
    group_by(across(all_of(identical_cols))) %>%
    mutate(MultipleID = ifelse(n() > 1, "Yes", "No")) %>%
    ungroup()

  # Group the rows by the columns specified above and find the row with the highest DeltaScore in each group
  df <- df %>%
    group_by(across(all_of(identical_cols))) %>%
    mutate(DMax = ifelse(MultipleID == "Yes" & DeltaScore == max(DeltaScore), "Max", NA_character_),
           Multiple = ifelse(DMax == "Max", "Max", "Min")) %>%
    ungroup()

  # Filter the rows based on the specified criteria (MultipleID == "Yes" and DMax == "Max" or MultipleID == "No")
  df <- df %>%
    filter(MultipleID == "Yes" & DMax == "Max" | MultipleID == "No")

  return(df)
}



OG_pd_data<-pd_data


pd_data<- remove_dup(pd_data)


#Map acquisition file names
###########################################################################
# function to clean and parse data from PD output file
annotate_features <- function(raw_data){

  # Changing column header to match FASTA file
  raw_data <- raw_data %>%
    rename('UniprotID' = 'Protein Accessions')


  #Remove Conserved Proteins
  raw_data<- raw_data%>%
    filter(!grepl(";", raw_data$UniprotID)
    )


  #Remove Carbamidomethyls

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


  # Mod Position*
  raw_data$ModPositionL <- sub("^([[:alpha:]]*).*", "\\1",
                               raw_data$Modifications)
  raw_data$ModPositionN <- as.numeric(gsub(".*?([0-9]+).*", "\\1",
                                           raw_data$Modifications))

  raw_data$ModPosition <- ifelse(raw_data$ModPositionL == "NA", "NA",
                                 paste(raw_data$ModPositionL,
                                       raw_data$ModPositionN))
  raw_data$ModPosition <- gsub("^; ", "", raw_data$ModPosition)


  return(raw_data)
}

# parse the PD data
pd_data_annotated <- annotate_features(pd_data)
