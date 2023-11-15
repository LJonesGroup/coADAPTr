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
FASTA <- read.fasta(file = ("./data/TNFa_Adalimum.fasta"),
                    clean_name = FALSE)
#SWISSPROT Homo Sapian 9606 + cRAP Reviewed.fasta
#USER INDICATES FILE PATH
file_path<- "C:/Users/Raqie/Desktop/Data/TNF-alphaGN.xlsx"
#UMBPSC/RCoding/JonesLabScripts/data/100 mM H2O2 1 Min BR1 Core.xlsx"
#"C:/Users/Raqie/Desktop/UMBPSC/RCoding/JonesLabScripts/data/100 mM H2O2 1 Min BR1 Outer.xlsx"

# Set output directory
file_output= "C:/Users/Raqie/Desktop/AutoDataAnalysis/Building"



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

# function to parse the FASTA file for later manipulations

parse_fasta <- function(fasta_in){

  # Rename Columns
  fasta_in <- fasta_in %>%
    rename(protein_sequence = seq.text)
  fasta_in<- fasta_in %>%
    rename(MasterProteinAccessions = seq.name)

  # Isolate MPA
  fasta_in$UniprotID <- gsub("^.+\\|(\\w+)\\|.*$", "\\1", fasta_in$MasterProteinAccessions)

  return(fasta_in)

}

# rename columns of the FASTA file for merging dataframes later
FASTA <- parse_fasta(FASTA)


#####################################################################
# function to locate the residue number for the peptide termini and residues

locate_startend_res <- function(raw_data){

  uniqueMPA <- unique(raw_data[, c("Master Protein Accessions",'Sequence')])
  uniqueMPA <- as.data.frame(uniqueMPA)

  # merge MPA and fasta dataframes together
  raw_data <- merge(raw_data, FASTA, by = "UniprotID")

  # get index of peptide sequence from protein sequence
  index <- str_locate(raw_data$protein_sequence, raw_data$Sequence)
  raw_data <- cbind(raw_data, index)  # TODO: memory efficiency?

  #annotate peptide
  raw_data$peptide<- paste(raw_data$start,"-", raw_data$end)

  raw_data$mod_count <- str_count(raw_data$Modifications, "\\(.*?\\)")
  raw_data$mod_count <- ifelse(raw_data$MOD == "Unoxidized", 0, raw_data$mod_count)

  # get modified residue and number
  raw_data$mod_res <- ifelse(raw_data$ModPositionN>0, raw_data$start + raw_data$ModPositionN - 1, NA)

  #combine the mod AA and res
  raw_data$Res<- paste(raw_data$ModPositionL, raw_data$mod_res)


  return(raw_data)
}



# annotate the PD output file with peptide and modification residue positions ERROR HERE. LOTS OF DATA LOST DURING MERGE
pd_data_fasta_merged <- locate_startend_res(pd_data_annotated)
# Assuming pd_data_fasta_merged is your data frame


# garbage cleanup
rm(pd_data_annotated);gc()
#rm(FASTA);gc()


# FPOP Calculations ---------------------------------------------------------------------------------------------
#Calculating the total peptide areas (Oxidized and unoxidized)
area_calculations_pep <- function(df_in) {
  df_out <- df_in %>%
    group_by(MasterProteinAccessions, Sequence, SampleControl, MOD) %>%
    reframe(TotalArea = sum(`Precursor Abundance`) ) %>%
    ungroup() %>%
    pivot_wider(
      id_cols = c("MasterProteinAccessions", "Sequence"),
      names_from = c("SampleControl", "MOD"),
      values_from = "TotalArea",
      values_fill = NA
    )

  df_out <- df_out %>%
    rename(
      Control_OxidizedArea = Control_Oxidized,    # Rename OldColumn1 to NewColumn1
      Control_UnoxidizedArea = Control_Unoxidized,    # Rename OldColumn2 to NewColumn2
      Sample_OxidizedArea = Sample_Oxidized,     # Rename OldColumn3 to NewColumn3
      Sample_UnoxidizedArea = Sample_Unoxidized     # Rename OldColumn3 to NewColumn3
    )

  #create conditions so that the bare minimum criteria for FPOP can be calculated
  #specifically cases where sample oxidized area is detected and control oxidized area re detected.
  #Cases where that peptide was not detected in the sample unoxidized area and control oxidized are will be turned to 0
  ###*****MAKE SURE THIS IS OK WITH LISA WANT TO AVOID IMPLICIT BIAS
  #The peptide at least has to be detected by the mass spec in order to calculate control total area?
  df_out$Control_OxidizedArea <- ifelse(df_out$Control_UnoxidizedArea > 0 & is.na(df_out$Control_OxidizedArea), 0, df_out$Control_OxidizedArea)

  df_out$Sample_UnoxidizedArea <- ifelse(df_out$Sample_OxidizedArea > 0 & is.na(df_out$Sample_UnoxidizedArea), 0, df_out$Sample_UnoxidizedArea)

  #Calculating the extent of modification (EOM) for the sample and control data for each peptide

  df_out$TotalSampleArea <- rowSums(df_out[, c("Sample_OxidizedArea", "Sample_UnoxidizedArea")])
  df_out$TotalControlArea <- rowSums(df_out[, c("Control_OxidizedArea", "Control_UnoxidizedArea")])


  #Calculating the extent of modification (EOM) by subtracting background oxidation (control)
  df_out$EOMSample <- df_out$Sample_OxidizedArea / df_out$TotalSampleArea
  df_out$EOMControl <- df_out$Control_OxidizedArea / df_out$TotalControlArea
  df_out$EOM <- df_out$EOMSample - df_out$EOMControl

  #Calculating N or the number of times the petpdide was detected by the MS
  df_out$N <- df_in %>%
    group_by(MasterProteinAccessions, Sequence) %>%
    count()

  df_out$N$Sequence <- NULL
  df_out$N$MasterProteinAccessions <- NULL
  colnames(df_out)[12] <- "N"

  #Calcculating standard deviation
  test <- pd_data_fasta_merged %>%
    group_by(MasterProteinAccessions, Sequence) %>%
    reframe(sdprep = sd(`Precursor Abundance`))
  df_out$sdprep<- test$sdprep
  df_out$SD<- df_out$SD <- df_out$sdprep/(df_out$TotalSampleArea+df_out$TotalControlArea)

  # Renaming the 12th column to "N"


  return(df_out)
}

#Make sure that if the data detected at least and not modified calculations still occur or is it ok

Areas_pep <- area_calculations_pep(pd_data_fasta_merged)


##############################################################################
# Data Grabbing ---------------------------------------------------------------
###############################################################################

# function to subset sequence metadata like residue start/stop
grab_seq_metadata_pep <- function(df_in){

  df_out <- subset(df_in, select = c("MasterProteinAccessions", "Sequence",
                                     "peptide", "start"))
  df_out <- df_out[!duplicated(df_out), ]

  return(df_out)
}
#####Merge these variable with a function
# merge metadata with numeric graphing data
graphing_df_pep <- Areas_pep %>%
  left_join(grab_seq_metadata_pep(pd_data_fasta_merged))


graphing_df_pep <- graphing_df_pep[order(graphing_df_pep$start), ]
graphing_df_pep$MasterProteinAccessions <- gsub(".*\\|(.*?)\\|.*", "\\1", graphing_df_pep$MasterProteinAccessions)


#Create function for filtering graphical data
#The above data is complete but, not all of the data meets standards for graphic. ie error bars greater than EOM

# Filter Graphical Data to include data that is acceptable.
filtered_graphing_df_pep <- function(df_in) {
  df_in <- df_in %>% filter(EOM > 0)
  df_in <- df_in %>% filter(EOM > SD)
  df_in <- df_in %>% filter(df_in$N >= 4)  # Filter for the 12th column > 4
  df_in <- df_in %>% arrange(start)
  return(df_in)
}



quant_graph_df_pep <- filtered_graphing_df_pep(graphing_df_pep)


##
#Residue Level-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Remove multiply modified speceis
#RefactoredCorrect
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

  #create conditions so that the bare minimum criteria for FPOP can be calculated
  #specifically cases where sample oxidized area is detected and control oxidized area re detected.
  #Cases where that peptide was not detected in the sample unoxidized area and control oxidized are will be turned to 0

  #The peptide at least has to be detected by the mass spec in order to calculate control total area?
  df_out$Control_Oxidized <- ifelse(df_out$Control_Unoxidized > 0 & df_out$Control_Oxidized == "NA", 0, df_out$Control_Oxidized)

  #For extent of modification calculations the peptide needs to be modified
  df_out$Sample_Unoxidized <- ifelse(df_out$Sample_Oxidized > 0 & df_out$Sample_Unoxidized == "NA", 0, df_out$Control_Oxidized)



  df_out<- df_out %>%
    mutate(SampleTotalArea = Sample_Oxidized + Sample_Unoxidized,
           ControlTotalArea = Control_Oxidized + Control_Unoxidized)


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

  #FIX CALCULATE N FUNCTION FOR RESIDUE LEVEL N > $
  # Calculate N values and store in a separate data frame
  N_df <- df_in %>%
    filter(mod_count == 0| mod_count == 1, !is.na(Res)) %>%
    group_by(MasterProteinAccessions, Sequence, Res) %>%
    summarize(N = n())  # Count the occurrences

  # Merge the N column into df_out
  df_out$N <- df_out %>%
    left_join(N_df, by = c("MasterProteinAccessions", "Sequence", "Res")) %>%
    pull(N)  # Extract N column
  # colnames(df_out)[12] <- "N"  # Renaming the 12th column to "N"

  # Calculate standard deviation and store in a separate data frame
  sd_df <- df_in %>%
    group_by(MasterProteinAccessions, Sequence, Res) %>%
    summarize(sdprep = sd(`Precursor Abundance`))

  # Join the calculated sdprep values to df_out
  df_out <- df_out %>%
    left_join(sd_df, by = c("MasterProteinAccessions", "Sequence"))
  df_out$SD <- df_out$sdprep/(df_out$SampleTotalArea+df_out$ControlTotalArea)

  return(df_out)
}

Areas_res <- area_calculations_res(pd_data_fasta_merged)


##############################################################################
# Data Grabbing ---------------------------------------------------------------
###############################################################################

# function to subset sequence metadata like residue start/stop
grab_seq_metadata_res <- function(df_in){
  df_out <- subset(df_in, select = c("MasterProteinAccessions", "Sequence", "Res", "start"))
  df_out <- df_out[!duplicated(df_out), ]
  return(df_out)
}
# merge metadata with numeric graphing data
graphing_df_res <- Areas_res%>%
  left_join(grab_seq_metadata_res(pd_data_fasta_merged))

graphing_df_res <- graphing_df_res[!(is.na(graphing_df_res$Res) | graphing_df_res$Res == ""), ]

# Ascending order
graphing_df_res <- graphing_df_res[order(graphing_df_res$start), ]
graphing_df_res$MasterProteinAccessions <- gsub(".*\\|(.*?)\\|.*", "\\1", graphing_df_res$MasterProteinAccessions)

#Add filters for residue level graphing. Similar to above.
#filter for graphing data
filtered_graphing_df_res <- function(df_in) {
  df_in <- df_in %>%filter(EOM >0)
  df_in <- df_in %>% filter(EOM > SD)
  df_in <- df_in %>% filter(df_in[[12]] > 3)  # Filter for the 12th column > 4
  df_in <- df_in %>% arrange(start)
  return(df_in)
}
quant_graph_df_res<-filtered_graphing_df_res(graphing_df_res)








#creating a table of totals
create_totals_tablelist <- function(df_in, df_res) {
  # Initialize an empty data frame
  df_out <- data.frame(
    UniqueProteinMod = NA,
    UniqueSeqMod = NA,
    UniqueResMod = NA,
    QuantifiableModProtein = NA,
    QuantifiableModSeq = NA,
    QuantifiableModRes = NA
  )

  # Count the number of unique MasterProteinAccessions values
  unique_protein_mod <- unique(df_in$MasterProteinAccessions)
  df_out$UniqueProteinMod <- length(unique_protein_mod)

  # Count the number of unique Sequence values
  unique_seq_mod <- unique(df_in$Sequence)
  df_out$UniqueSeqMod <- length(unique_seq_mod)

  # Count the number of unique Residue values
  unique_res_mod <- unique(df_res$Res)
  df_out$UniqueResMod <- length(unique_res_mod)



  # Select rows where EOM is greater than 0 and greater than SD and SD is not NA
  df_filtered <- df_in %>%
    filter(EOM > 0 & EOM > SD & !is.na(SD))
  df_filteredres <- df_res %>%
    filter(EOM > 0 & EOM > SD & SD > 0)

  # Count the number of unique MasterProteinAccessions values in the filtered data frame
  quantifiable_protein_mod <- unique(df_filtered$MasterProteinAccessions)
  df_out$QuantifiableModProtein <- length(quantifiable_protein_mod)


  # Count the number of unique Sequence values in the filtered data frame
  quantifiable_seq_mod <- unique(df_filtered$Sequence)
  df_out$QuantifiableModSeq <- length(quantifiable_seq_mod)


  # Count the number of unique Residue values in the filtered data frame
  quantifiable_res_mod <- unique(df_filteredres$Res)
  df_out$QuantifiableModRes <- length(quantifiable_res_mod)


  # Return the resulting data frame
  return(df_out)
}

TotalsTable<-create_totals_tablelist(graphing_df_pep, graphing_df_res)







####Saving tables---------------------------------------------------------------------------
save_data_frames <- function(file_path, output_directory, ...) {
  # Extract the file name from the file path
  file_name <- basename(file_path)

  # Remove the extension from the file name
  file_name <- tools::file_path_sans_ext(file_name)

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
# Save the data frames as separate Excel files
save_data_frames(file_path, file_output, TotalsTable = TotalsTable, quant_graph_df_pep = quant_graph_df_pep, quant_graph_df_res = quant_graph_df_res, graphing_df_pep = graphing_df_pep, graphing_df_res=graphing_df_res)




# Saving Plots -----------------------------------------------------------------------

# function to generate extent of modification figures

generate_eom_plot_pep <- function(df_in, file_output, excel_filename) {
  df_in <- df_in %>%
    arrange(start)

  # Create a factor variable to represent the sorted order
  df_in$peptide <- factor(df_in$peptide, levels = df_in$peptide)
  # Grab the protein that is being plotted
  protein <- unique(df_in$MasterProteinAccessions)
  # Generate a bargraph of the extent of modification for each peptide
  # that maps to this protein
  fig <- ggplot(df_in, aes(x = peptide, y = EOM)) +
    geom_bar(position = "dodge", stat = "identity") +
    geom_errorbar(aes(ymin = EOM - SD, ymax = EOM + SD), width= .5, linewidth= 1) +
    labs(title = paste(protein," Peptide Level Analysis"),
         x = "Peptide",
         y = "Extent of Modification") +
    theme_classic() +
    theme(text = element_text(size = 18, family = "Arial")) +
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 20,
                                    family = "Arial",
                                    face = "bold")) +
    scale_fill_manual(values = c("grey42"))

  # Create the output directory for bar graphs based on the file output and excel filename
  graph_output_directory <- file.path(file_output, paste0(excel_filename, "_PeptideLevelBarGraphs"))
  dir.create(graph_output_directory, showWarnings = FALSE, recursive = TRUE)

  # Generate the full file path for this protein and save the figure
  file_out <- file.path(graph_output_directory, paste0(protein, ".png"))
  ggsave(filename = file_out, plot = fig, device = "png")

  # Print a message to indicate successful saving
  cat("Bar graph for", protein, "saved as", file_out, "\n")
}

# Get the excel file name from the file path
excel_filename <- tools::file_path_sans_ext(basename(file_path))

# Iterate over each protein and make an extent of mod plot for it
for (protein in quant_graph_df_pep$MasterProteinAccessions) {
  # subset the dataframe for this protein
  temp <- subset(quant_graph_df_pep, MasterProteinAccessions == protein)
  generate_eom_plot_pep(temp, file_output, excel_filename)
}

####################################################################################################


generate_eom_plot_res <- function(df_in, file_output, excel_filename) {
  df_in <- df_in %>%
    arrange(start)

  # Grab the protein that is being plotted
  protein <- unique(df_in$MasterProteinAccessions)
  # Generate a bargraph of the extent of modification for each peptide
  # that maps to this protein
  fig <- ggplot(df_in, aes(x = Res, y = EOM)) +
    geom_bar(position = "dodge", stat = "identity") +
    geom_errorbar(aes(ymin = EOM - SD, ymax = EOM + SD), width= .5, linewidth= 1) +
    labs(title = paste(protein," Residue Level Analysis"),
         x = "Residue",
         y = "Extent of Modification") +
    theme_classic() +
    theme(text = element_text(size = 18, family = "Arial")) +
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 20,
                                    family = "Arial",
                                    face = "bold")) +
    scale_fill_manual(values = c("grey42"))

  # Create the output directory for bar graphs based on the file output and excel filename
  graph_output_directory <- file.path(file_output, paste0(excel_filename, "_ResidueLevelBarGraphs"))
  dir.create(graph_output_directory, showWarnings = FALSE, recursive = TRUE)

  # Generate the full file path for this protein and save the figure
  file_out <- file.path(graph_output_directory, paste0(protein, ".png"))
  ggsave(filename = file_out, plot = fig, device = "png")

  # Print a message to indicate successful saving
  cat("Bar graph for", protein, "saved as", file_out, "\n")
}

# Get the excel file name from the file path
excel_filename <- tools::file_path_sans_ext(basename(file_path))

# Iterate over each protein and make an extent of mod plot for it
for (protein in quant_graph_df_res$MasterProteinAccessions) {
  # subset the dataframe for this protein
  temp <- subset(quant_graph_df_res, MasterProteinAccessions == protein)
  generate_eom_plot_res(temp, file_output, excel_filename)
}



dev.off()
