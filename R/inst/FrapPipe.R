####add functions for grouped bar grpah and venn diagrams
options(stringsAsFactors = FALSE)

# list of required packages
required_packages <- c("ExcelFunctionsR", "remotes", "data.table", "stringr", "protr",
                       "plyr", "extrafont", "readxl", "ggplot2", "eulerr",
                       "tidyverse", "EnvStats", "dplyr", "writexl", "conflicted",
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
#Resolve library conflicts if necessary

conflict_prefer("arrange", winner = "dplyr")
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

###########################################################################
#Homo sapien Reviewed 12062021
# Read in required inputs
#FASTA <- read.fasta(file = ("./data/Homo sapien Reviewed 12062021.fasta"),
                    #clean_name = FALSE)
#USER INDICATES FILE PATH
#file_path<- "C:/Users/rasho/Desktop/coADAPTR/MTXfpop_psm.tsv"


# Set output directory
file_output= "C:/Users/rasho/Desktop/AutoDataAnalysis"

#Read in data from FragPipe
#fragpipedata<- read.table(file_path)
fragpipedata<- read.delim("C:/Users/rasho/Desktop/coADAPTR/MTXfpop_psm.tsv", header = TRUE)



####TMT data Removign Extra Columns
columns_to_remove <- c(
  "Spectrum", "Spectrum.File", "Modified.Peptide",
  "Extended.Peptide", "Prev.AA", "Next.AA", "Peptide.Length",
  "Charge", "Retention", "Observed.Mass", "Calibrated.Observed.Mass",
  "Observed.M.Z", "Calibrated.Observed.M.Z", "Calculated.Peptide.Mass",
  "Calculated.M.Z", "Delta.Mass", "Expectation", "Hyperscore",
  "Nextscore", "PeptideProphet.Probability", "Number.of.Enzymatic.Termini",
  "Number.of.Missed.Cleavages", "Intensity", "Assigned.Modifications",
  "Observed.Modifications", "Protein.Description", "Gene", "MSFragger.Localization",
  "Best.Score.with.Delta.Mass", "Best.Score.without.Delta.Mass",
  "Purity", "Is.Unique", "Mapped.Genes", "Mapped.Proteins", "Quan.Usage"
)

remove_columns <- function(raw_data, columns_to_remove) {
  # Create a duplicate of the original dataframe
  new_data <- raw_data

  # Remove specified columns
  new_data <- new_data[, !names(new_data) %in% columns_to_remove]

  return(new_data)
}

# Create a new dataframe without the specified columns
fpdata <- remove_columns(fragpipedata, columns_to_remove)


#####################################################################
# function to locate the residue number for the peptide termini and residues
fpannotate <- function(raw_data) {
  # Extract everything before the first parenthesis in the "FPOP Modifications" column
  #USE THIS FOR RESIDUE LEVEL CALCULATIONS!!!!!!!!!!!
  #raw_data$Residue <- gsub("\\(.*", "", raw_data$FPOP.Modifications)

  # Create a new column to store the modified peptide
  raw_data$peptidelocation <- paste(raw_data$Protein.Start, "-", raw_data$Protein.End)
  raw_data$MOD<- ifelse(is.na(raw_data$FPOP.Modifications) | raw_data$FPOP.Modifications == "", "Unoxidized",
                                          ifelse(grepl("()", raw_data$FPOP.Modifications), "Oxidized", "Unoxidized"))
  raw_data$MOD <- ifelse(is.na(raw_data$FPOP.Modifications) | raw_data$FPOP.Modifications == "" | grepl("^\\s*$", raw_data$FPOP.Modifications) | !grepl("[A-Za-z]", raw_data$FPOP.Modifications), "Unoxidized", raw_data$MOD)


  return(raw_data)
}

fpdata <- fpannotate(fpdata)



# Data Cleaning Calculations ---------------------------------------------------------------------------------------------
#Summing multiple occurrences
summarize_duplicates <- function(data) {
  # Step 1: Remove the FPOP.Modifications column
  data <- select(data, -FPOP.Modifications)

  # Step 2: Auto-detect columns containing VL, VC, DL, and DC
  cols_to_sum <- grep("(VL|VC|DL|DC)\\d*$", names(data), value = TRUE)

  # Step 3: Group the data by specified columns and calculate the sum of detected columns
  summarized_data <- data %>%
    group_by(Peptide, Protein.Start, Protein.End, Protein, Protein.ID, Entry.Name, MOD, peptidelocation) %>%
    summarize(across(all_of(cols_to_sum), sum, na.rm = TRUE), .groups = "drop")

  summarized_data<- summarized_data %>%
    group_by(Peptide, Protein.Start, Protein.End, Protein, Protein.ID, Entry.Name, peptidelocation) %>%
    filter(all(c("Oxidized", "Unoxidized") %in% MOD))

  return(summarized_data)
}


# Usage example:
fpdata_sums <- summarize_duplicates(fpdata)

#Spread the data out before math
spread_unoxidized_data <- function(data) {
  spread_data <- data %>%
    mutate(across(starts_with(c("VL", "VC", "DL", "DC")), ~if_else(MOD == "Unoxidized", as.numeric(.), .))) %>%
    pivot_wider(names_from = MOD, values_from = starts_with(c("VL", "VC", "DL", "DC")), names_sep = "_")

  return(spread_data)
}
spreadfpdata<- spread_unoxidized_data(fpdata_sums)

####FPOP Calculations (Peptide Level)---------------------------------------------------------------------------------------------
##############################################################################

EOM_pep <- function(df_in) {
  df_in$VCOxidized<- df_in$VC1_Oxidized + df_in$VC2_Oxidized #+ df_in$VC3_Oxidized
  df_in$DCOxidized<- df_in$DC1_Oxidized + df_in$DC2_Oxidized #+ df_in$DC3_Oxidized
  df_in$VLOxidized<- df_in$VL1_Oxidized + df_in$VL2_Oxidized + df_in$VL3_Oxidized
  df_in$DLOxidized<- df_in$DL1_Oxidized + df_in$DL2_Oxidized + df_in$DL3_Oxidized

  df_in$VCUnoxidized<- df_in$VC1_Unoxidized + df_in$VC2_Unoxidized #+ df_in$VC3_Unoxidized
  df_in$DCUnoxidized<- df_in$DC1_Unoxidized + df_in$DC2_Unoxidized #+ df_in$DC3_Unoxidized
  df_in$VLUnoxidized<- df_in$VL1_Unoxidized + df_in$VL2_Unoxidized + df_in$VL3_Unoxidized
  df_in$DLUnoxidized<- df_in$DL1_Unoxidized + df_in$DL2_Unoxidized + df_in$DL3_Unoxidized

  df_in$VLEOM <- (df_in$VLOxidized) / (df_in$VLOxidized + df_in$VLUnoxidized)
  df_in$VCEOM <- (df_in$VCOxidized) / (df_in$VCOxidized + df_in$VCUnoxidized)

  df_in$DLEOM <- (df_in$DLOxidized) / (df_in$DLOxidized + df_in$DLUnoxidized)
  df_in$DCEOM <- (df_in$DCOxidized) / (df_in$DCOxidized + df_in$DCUnoxidized)

  df_in$VEOM <- df_in$VLEOM - df_in$VCEOM
  df_in$DEOM <- df_in$DLEOM - df_in$DCEOM
#add back for full dat set #2nd ine "+ df_in$DC3_Oxidized)"
  #first line "+ df_in$VC3_Oxidized)"
  df_in$VCSD<- sd(df_in$VC1_Oxidized + df_in$VC2_Oxidized) /(df_in$VCOxidized + df_in$VCUnoxidized)
  df_in$DCSD<- sd(df_in$DC1_Oxidized + df_in$DC2_Oxidized) / (df_in$DCOxidized + df_in$DCUnoxidized)
  df_in$VLSD<- sd(df_in$VL1_Oxidized + df_in$VL2_Oxidized + df_in$VL3_Oxidized)/ (df_in$VLOxidized + df_in$VLUnoxidized)
  df_in$DLSD<- sd(df_in$DL1_Oxidized + df_in$DL2_Oxidized + df_in$DL3_Oxidized)/ (df_in$DLOxidized + df_in$DLUnoxidized)
  df_in$VSD<- df_in$VLSD - df_in$VCSD
  df_in$DSD<- df_in$DLSD - df_in$DCSD


  df_out<- df_in%>%
    subset(select = c("Protein", "Protein.ID", "Peptide", "peptidelocation", "VLEOM", "VCEOM", "DLEOM", "DCEOM", "VEOM", "DEOM", "VLSD", "VCSD", "DLSD", "DCSD",
                      "VSD", "DSD"))

  return(df_out)
}


EOM <- EOM_pep(spreadfpdata)




###########################################################################################################

# Filter Graphical Data to include data that is acceptable.
filtered_graphing_df_pep <- function(df_in) {
  df_in <- df_in %>% filter(VEOM > 0)
  df_in <- df_in %>% filter(DEOM > 0)
  df_in <- df_in %>% filter(VEOM > SD)
  df_in <- df_in %>% filter(DEOM > SD)
  df_in <- df_in %>% arrange(Protein)
  return(df_in)
}



quant_graph_df_pep <- filtered_graphing_df_pep(EOM)


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
    left_join(sd_df, by = c("MasterProteinAccessions", "Sequence", "Res"))
  df_out$SD <- df_out$sdprep/(df_out$SampleTotalArea+df_out$ControlTotalArea)

  df_out<- df_out[complete.cases(df_out[c("Res", "EOM")]), ]

  return(df_out)
}

Areas_res <- area_calculations_res(pd_data_fasta_merged)


#rebuilding areas_res function##########################################################





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
graphing_data_res <- function(df_in) {
  df_out <- df_in %>%
    left_join(grab_seq_metadata_res(pd_data_fasta_merged)) %>%
    filter(!(is.na(Res) | Res == "")) %>%
    arrange(start) %>%
    mutate(MasterProteinAccessions = gsub(".*\\|(.*?)\\|.*", "\\1", MasterProteinAccessions))

  return(df_out)
}

graphing_df_res<- graphing_data_res(Areas_res)
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
