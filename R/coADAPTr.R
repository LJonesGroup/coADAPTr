# Setup -------------------------------------------------------------------

options(stringsAsFactors = FALSE)

# list of required packages
required_packages <- c("ExcelFunctionsR", "remotes", "data.table", "stringr", "protr",
                       "plyr", "extrafont", "readxl", "ggplot2", "eulerr", 
                       "tidyverse", "EnvStats", "dplyr", "writexl","stringr", 
                       "phylotools", "parallel", "rlist","argparser")

# load or install&load all packages 
# source: 
# https://vbaliga.github.io/verify-that-r-packages-are-installed-and-loaded/
package.check <- lapply(
  required_packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, type= source, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

# set working directory to current directory
setwd(getwd())


###########################################################################

# Read in required inputs
FASTA <- read.fasta(file = ("./data/Homo sapien Reviewed 12062021.fasta"), 
                    clean_name = FALSE)  
#USER INDICATES FILE PATH
file_path<- "C:/Users/Raqie/Desktop/UMBPSC/RCoding/JonesLabScripts/data/coADAPTrTestData.xlsx"

pd_data <- read_excel(file_path)

# Set output directory
fileoutput= "./results/OPHEKDigest.xlsx"  

#Does this need to have a separate function or can it be merged with Annotate Features. 
# Annotate Sample/Control based on MS acquisitoin filename convention
pd_data$SampleControl <- ifelse(pd_data$'Spectrum File' %like% 
                                  "NL", "Control", "Sample") 

#Keep lowest delta score from seq
#Remove unused identifying node column   KEEP COLUMN TO SEE THAT DIFFERENT NODES PRODUCE DIFFERENT DELSTA SCORE FOR OTHERWISE IDENTICAL PEPTIDE

remove_pd_deplicates<- function (df_in){
  df_in = subset(df_in, select = c("Identifying Node",  "Sequence", "Modifications", "Master Protein Accessions",
  "Protein Accessions", "Spectrum File", "Precursor Abundance",  "SampleControl"))
  df_in = unique(df_in)
  return (df_in)
}

pd_data<- remove_pd_deplicates(pd_data)

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
  
  #THE CELLS WITH OTHER MODS ARE REMOVED!!!!!FIX HERE
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
  
  
  
  #raw_data$Modifications <- ifelse(grepl("Carbamidomethyl",
                                         #raw_data$Modifications),
                                   #NA, raw_data$Modifications)
  
  #raw_data$MOD <- ifelse(is.na(raw_data$Modifications) | raw_data$Modifications == "", "Unoxidized", 
  #                       ifelse(grepl("Oxidation", raw_data$Modifications), "Oxidized", raw_data$Modifications))
  
  
  
  #raw_data$MOD <- ifelse(raw_data$Modifications == is.na(raw_data$Modifications),"Unoxidized","Oxidized")
  #raw_data$MOD <- ifelse(is.na(raw_data$Modifications), "Unoxidized", ifelse(grepl("Oxidation", raw_data$MOD), "Oxidized", raw_data$MOD))
  #raw_data$MOD[is.na(raw_data$MOD)] <- "Unoxidized"
  
  
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
  
  # get number of modified residues
  raw_data$mod_count <- str_count(raw_data$Modifications, ";")
  
  # get modified residue and number
  raw_data$mod_res <- ifelse(raw_data$ModPositionN>0, raw_data$start + raw_data$ModPositionN - 1, NA)
  
  #combine the mod AA and res
  raw_data$Res<- paste(raw_data$ModPositionL, raw_data$mod_res)
  
  return(raw_data)
}



# annotate the PD output file with peptide and modification residue positions ERROR HERE. LOTS OF DATA LOST DURING MERGE
pd_data_fasta_merged <- locate_startend_res(pd_data_annotated)

# garbage cleanup
rm(pd_data_annotated);gc()
rm(FASTA);gc()


# FPOP Calculations -----------------------------------------------------------

# functions for summing areas
calculate_area_pep <- function(df_in){
  
  # Calculate the sum of precursor areas for each Protein-Peptide pair
  df_out <- df_in %>%
    group_by(MasterProteinAccessions, Sequence) %>%
    summarise(Area = sum(`Precursor Abundance`))
  
  return(df_out)
}

calculate_area_ox_pep <- function(df_in){
  
  # Calculate the precusor area for each unique combination of 
  # Protein-Peptide-SampleGroup-Modification 
  df_out <- df_in %>%
    group_by(MasterProteinAccessions, Sequence, SampleControl, MOD) %>%
    summarise(Oxidized.Area = sum(`Precursor Abundance`)[MOD == "Oxidized"])
  
  # Reformat the dataframe to match conventions
  df_out <- df_out[complete.cases(df_out), ]
  df_out$MOD <- NULL
  df_out <- spread(df_out, SampleControl, Oxidized.Area)
  df_out <- rename(df_out, OxidizedSampleArea = Sample)
  df_out <- rename(df_out, OxidizedControlArea = Control)
  
  return(df_out)
}

# function to calculate number of sequences in original dataframe 
#CAN BE USED AS PREFILTERING MECHANISM TOO
calculate_n_pep <- function(df_in) {
  
  # Calculate the number of unique Protein-Peptide pairs
  NCalc <- df_in %>%
    group_by(MasterProteinAccessions, Sequence) %>%
    count(df_in$Sequence)
  
  return(NCalc)
}
#Check error here, Sample col name
# function for calculating the extent of extent of modification, EOM
calculate_eom_pep <- function(areas, ox_areas){
  
  # Merge the unmodified and modified areas dataframes
  mod_data <- merge(areas, ox_areas, 
                    by = c("MasterProteinAccessions", "Sequence"), all = TRUE)
  mod_data <- mod_data[complete.cases(mod_data), ]
  
  # Calculate extent of modification = oxidized area / total area
  mod_data$EOMSample <- mod_data$OxidizedSampleArea / mod_data$Area
  mod_data$EOMControl <- mod_data$OxidizedControlArea / mod_data$Area
  mod_data$EOM <- mod_data$EOMSample-mod_data$EOMControl

  # replace missing data with zero
  #mod_data[is.na(mod_data)] <- 0
  
  return(mod_data)
}

area_df_pep <- calculate_area_pep(pd_data_fasta_merged)
area_ox_df_pep <- calculate_area_ox_pep(pd_data_fasta_merged)  

#Can the 2nd line be merged into the function? What about the other lines?
NCalcPep <- calculate_n_pep(pd_data_fasta_merged)
NCalcPep$`df_in$Sequence` <- NULL

area_ox_df_pep<- area_ox_df_pep%>%
  left_join(NCalcPep)

mod_df_pep <- calculate_eom_pep(area_df_pep, area_ox_df_pep)


# garbage cleanup
rm(area_df_pep, area_ox_df_pep);gc()

######Calculating standard deviation (SD)
calculateSD <- function(df_in){
  
  # Calculate oxidized area sd() for each Protein-Peptide in the sample
  df_out <- df_in %>%
    group_by(MasterProteinAccessions, Sequence, SampleControl, MOD) %>%
    filter(SampleControl == "Sample", MOD %like% "Oxidized") %>%
    summarise(SDOxidized.Area = sd(`Precursor Abundance`))
  
  df_out$SampleControl<- NULL
  df_out$MOD<- NULL
  return(df_out)
}

SDcalc_pep <- calculateSD(pd_data_fasta_merged)

mod_df_pep <- mod_df_pep %>%
  left_join(SDcalc_pep)

mod_df_pep$SD<- mod_df_pep$SDOxidized.Area/mod_df_pep$Area

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
graphing_df_pep <- mod_df_pep %>%
  left_join(grab_seq_metadata_pep(pd_data_fasta_merged))

# Sort start amino acid position by ascending order
graphing_df_pep <- graphing_df_pep[order(graphing_df_pep$start), ]

# Save graphing df as "'File Name' Peptide Level Data Analysis" 
write_xlsx(graphing_df_pep, path = fileoutput)

#Totals table for saving as well.
create_totals_table <- function(df_in) {
  # Initialize an empty data frame
  df_out <- data.frame(UniqueProteinMod = NA,
                       UniqueSeqMod = NA,
                       QuantifiableModProtein = NA,
                       QuantifiableModPep = NA)
  
  # Count the number of unique MasterProteinAccessions values
  df_out$UniqueProteinMod <- n_distinct(unique(df_in$MasterProteinAccessions))
  
  # Count the number of unique Sequence values
  df_out$UniqueSeqMod <- n_distinct(unique(df_in$Sequence))
  
  # Select rows where EOM is greater than 0 and greater than SD and SD is not NA
  df_filtered <- df_in %>%
    filter(EOM > 0 & EOM > SD & !is.na(SD))
  
  # Count the number of unique MasterProteinAccessions values in the filtered data frame
  df_out$QuantifiableModProtein <- n_distinct(unique(df_filtered$MasterProteinAccessions))
  
  # Count the number of unique Sequence values in the filtered data frame
  df_out$QuantifiableModPep <- n_distinct(unique(df_filtered$Sequence))
  
  # Return the resulting
  return(df_out)
}

TotalsTable<- (create_totals_table(graphing_df_pep))

#Create function for filtering graphical data 
#The above data is complete but, not all of the data meets standards for graphic. ie error bars greater than EOM

# Filter Graphical Data to include data that is acceptable. 
filtered_graphing_df_pep <- function(df_in) {
  df_in <- df_in %>% filter(EOM > 0)
  df_in <- df_in %>% filter(EOM > SD)
  df_in <- df_in %>%
    mutate(UniProtID = gsub("\\|", "", str_extract(MasterProteinAccessions, "\\|([A-Z0-9]+)\\|")))
  df_in <- df_in %>% arrange(start)
  return(df_in)
}


quant_graph_df_pep <- filtered_graphing_df_pep(graphing_df_pep)


# Plots -----------------------------------------------------------------------
#CHECK FUNCTION HERE. PLOTS NOT SAVED
# function to generate extent of modification figure
generate_eom_plot_pep <- function(df_in) {
  
  # Sort the data frame by the "start" column
  df_in <- df_in %>% arrange(start)
  
  # Grab the protein that is being plotted
  protein <- unique(df_in$UniProtID)
  
  # Generate a bargraph of the extent of modification for each peptide 
  # that maps to this protein
  fig <- ggplot(df_in, aes(x = peptide, y = EOM)) +
    geom_bar(position = "dodge", stat = "identity") +
    geom_errorbar(aes(ymin = EOM - SD, ymax = EOM + SD), width= .5, size= 1) +
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
  
  # Generate the full file path for this protein and save the figure
  fileout <- paste("results", protein, paste(protein,".png"), sep = "/") 
  ggsave(filename = fileout, plot = fig, device = "png")
}


# Iterate over each protein  and make an extent of mod plot for it
for (protein in quant_graph_df_pep$UniProtID) { 
  
  # subset the dataframe for this protein
  temp <- subset(quant_graph_df_pep, UniProtID == protein)
  plotout <- generate_eom_plot_pep(temp)
}






##
#Residue Level------------------------------------------------------------
#Remove multiply modified speceis
#RefactoredCorrect
calculate_area_res <- function(df_in){
  df_out <- df_in %>%
    group_by(MasterProteinAccessions, Sequence) %>%
    summarise(Area = sum(`Precursor Abundance`))
  
  print(df_out)
  return(df_out)
}
#miscalulation here!
#Ox area Correct 
calculate_area_ox_res <- function(df_in){
  df_out <- df_in %>%
    group_by(MasterProteinAccessions, Sequence, Res, `SampleControl`, mod_count, MOD) %>%
    summarise(Oxidized.Area = sum(`Precursor Abundance`)[MOD == "Oxidized" & mod_count == 0 ])
  df_out <- df_out[complete.cases(df_out), ]
  df_out$MOD <- NULL
  df_out$mod_count<- NULL
  df_out<- spread(df_out, SampleControl, Oxidized.Area)
  df_out<- rename(df_out, OxidizedSampleArea = Sample)
  df_out<- rename(df_out, OxidizedControlArea = Control)
  return(df_out)
}
#Functions to calculuate N
calculate_n_res <- function(df_in) {
  NCalc <- df_in %>%
    group_by(MasterProteinAccessions, Res, mod_count) %>%
    count(df_in$Res, mod_count = 0)
  
  return(NCalc)
}

#Refactored functions for clean up and EOM calc for saple and control
#Sike THe error happens with the merge. data gets scrambled.
calculate_eom_res <- function(areas, ox_areas){
  mod_data <- merge(area_df_res, area_ox_df_res, by = c("MasterProteinAccessions", "Sequence"), all = TRUE)
  print(mod_data)
  
  
  # extent of modification = oxidized area / total area
  mod_data$EOMSample <- mod_data$OxidizedSampleArea / mod_data$Area
  mod_data$EOMControl <- mod_data$OxidizedControlArea / mod_data$Area
  mod_data$EOM<- mod_data$EOMSample-mod_data$EOMControl
  print(mod_data)
  return(mod_data)
}

#Running Factored Functions

area_df_res <- calculate_area_res(pd_data_fasta_merged)
area_ox_df_res <- calculate_area_ox_res(pd_data_fasta_merged)
NCalcRes<- calculate_n_res(pd_data_fasta_merged)
NCalcRes$`df_in$Res`<- NULL
#Does this need to go in a function?
area_ox_df_res<- area_ox_df_res%>%
  left_join(NCalcRes)
mod_df_res <- calculate_eom_res(area_df_res, area_ox_df_res)
mod_df_res$`df_in$Res`<- NULL
mod_df_res$mod_count<- NULL
# garbage cleanup 
rm(area_df, area_ox_df);gc()

#Calculating SD
calculateSD <- function(df_in){
  df_out <- df_in %>%
    group_by(MasterProteinAccessions, Sequence, Res, SampleControl, MOD) %>%
    filter(SampleControl == "Sample", MOD %like% "Oxidized")%>%
    summarise(SDOxidized.Area = sd(`Precursor Abundance`))
  df_out$SampleControl<- NULL
  df_out$MOD<- NULL
  return(df_out)
}

sDcalc_res<- calculateSD(pd_data_fasta_merged)
mod_df_res<- mod_df_res%>%
  left_join(sDcalc_res)
mod_df_res$SD<- mod_df_res$SDOxidized.Area/mod_df_res$Area


# functions for generating data for making plots
# function subset sequence metadata like residue start/stop
grab_seq_metadata_res <- function(df_in){
  print(df_in)
  df_out <- subset(df_in, select = c("MasterProteinAccessions", "Sequence", "Res", "ModPositionN"))
  print(df_out)
  df_out <- df_out[!duplicated(df_out), ]
  print(df_out)
  return(df_out)
}
# merge metadata with numeric graphing data
graphing_df_res <- mod_df_res%>%
  left_join(grab_seq_metadata_res(pd_data_fasta_merged))


# Ascending order
graphing_df_res <- graphing_df_res[order(graphing_df_res$ModPositionN), ]

#Add filters for residue level graphing. Similar to above. 

#save graphing df as "'File Name' Peptide Level Data Analysis" 

write_xlsx(graphing_df_res, path = fileoutput)
#filter for graphing data
filtered_graphing_df_res <- function(df_in) {
  df_in <- df_in %>%filter(EOM >0)
  df_in <- df_in %>% filter(EOM > SD)
  df_in<- df_in%>%
    mutate(UniProtID = gsub("\\|", "", str_extract(MasterProteinAccessions, "\\|([A-Z0-9]+)\\|")))
  return(df_in)
}
quant_graph_df_res<-filtered_graphing_df_res(graphing_df_res)





#Plot
# function to generate extent of modification figure
# function to generate extent of modification figure
generate_eom_plot_res <- function(df_in) {
  
  # Grab the protein that is being plotted
  protein <- unique(df_in$UniProtID)
  residue<- unique(df_in$Res)
  # Generate a bargraph of the extent of modification for each peptide 
  # that maps to this protein
  fig <- ggplot(df_in, aes(x = Res, y = EOM)) +
    geom_bar(position = "dodge", stat = "identity") +
    geom_errorbar(aes(ymin = EOM - SD, ymax = EOM + SD), width= .5, size= 1) +
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
  
  # Generate the full file path for this protein and save the figure
  fileout <- paste("results", protein, paste(protein,".png"), sep = "/") 
  ggsave(filename = fileout, plot = fig, device = "png")
}

# Iterate over each protein  and make an extent of mod plot for it
for (protein in quant_graph_df_res$UniProtID) { 
  
  # subset the dataframe for this protein
  temp <- subset(quant_graph_df_res, UniProtID == protein)
  plotout <- generate_eom_plot_res(temp)
}

dev.off()

