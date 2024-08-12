# Setting Up ----------------------------------------------------------------------------------------------------
#Before Analyzing your data be sure to run this function to install necessary
#Packages, resolve package conflicts, import required data, and select result
#output folders. Follow the console prompts and anticipate a File Explorer up.

setup()


# Normalize TMT Data and Select Required Columns --------------------------------------------------------------------
#This function will normalize the TMT abundances and allow you to rename the
#required columns accordingly to prapre for annotation

TMT_Quant()

# Annotate the TMT data --------------------------------------------------------------------
#This function will locate the start and end residues to the TMT data.

TMT_Prep()

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

