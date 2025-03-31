# Setting Up ----------------------------------------------------------------------------------------------------
#Before Analyzing your data be sure to run this function to install necessary
#Packages, resolve package conflicts, import required data, and select result
#output folders. Follow the console prompts and anticipate a File Explorer pop up.

setup()

# Removing PD Generated Duplicates ---------------------------------------------------------------------------------------------------
#SKIP if PD 3.0>, Frag Pipe or FOXWare was used to search the data
#Run this if MS files were analyzed via PD <3.0
#PD Versions before 3.0 relied on a multi-level sequence searching algorithm
#that created duplicate identifications for the same peptide since multiple Sequest HT nodes
#were used to search the myriad of HRPF modifications

#Save the original raw data frame as a reference
OG_raw_data<-raw_data

raw_data<- remove_dup(raw_data)

#Check to ensure duplicated ID represented by different Sequest Nodes picking up the same MOD
#are removed

#Preparation ----------------------------------------------------------------------------------------------------
#Prepare the input data frame by selecting the relevant data, renaming the columns
#appropriately, and searching the FASTA file against peptide spectral matches to
#identify the start and end residues.

LFQ_Prep()


# EOM Calculations ---------------------------------------------------------------------------------------------
#Calculate the Extent of Modification (EOM) for each peptide and residue based on the LFQ data
#quant_graph_df containing all of the data that is acceptable for graphing
#Areas_pep and Areas_res contain the data that is used to calculate the EOM and will
#contain all data. This includes cases where the EOM could be negative (high background oxidation)
#or the SD is greater than the EOM (data has a high variance-likely due to experimental conditions)

EOM_Calculations()

#Saving Tables and Plots------------------------------------------------------------------------------------------------------------
#Save data frames as Excel files and save grouped bar from the plots as PNG files

Tables_and_Graphs()

#Creating Venn Diagrams ------------------------------------------------------------------------------------------------------------
#For the data you want in a Venn Diagram, create an excel file with headers of the Master Protein
#accessions that were modified in your data set. Save it as an excel file and run this function
venn_diagram()

dev.off()
