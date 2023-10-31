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
file_path<- "C:/Users/Raqie/Desktop/Data/BR1_EEC_MTX_5th.xlsx"
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

####TMT data

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
