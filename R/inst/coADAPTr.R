##Before Analyzing your data be sure to run this function to install necessary Packages
# and resolve package conflicts.
before_beginning()

###########################################################################

# Step 1-Read in required inputs
FASTA<- FASTA_file()

# Step-2 Set output directory
file_output<- output_folder()

# Step-3 Read in data from Proteome Discoverer

pd_data<- import_data()

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


# Step 7 Parse the FASTA file for later manipulations
FASTA<- parse_fasta(FASTA)

# Step 8 Locate the residue number for the peptide termini and residues

pd_data_fasta_merged <- locate_startend_res(pd_data_annotated)



# FPOP Calculations ---------------------------------------------------------------------------------------------
# Step 9 Calculating the total peptide areas and the extent of modification at the peptide level

Areas_pep<- area_calculations_pep(pd_data_fasta_merged)

library(dplyr)
library(tidyr)

testarea_calculations_pep <- function(df_in) {
  df_out <- df_in %>%
    group_by(MasterProteinAccessions, Sequence, SampleControl, MOD) %>%
    summarise(TotalArea = sum(`Precursor Abundance`), .groups = 'drop') %>%
    pivot_wider(
      id_cols = c("MasterProteinAccessions", "Sequence"),
      names_from = c("SampleControl", "MOD"),
      values_from = "TotalArea",
      values_fill = list(TotalArea = 0)
    )

  df_out <- df_out %>%
    rename(
      Control_OxidizedArea = Control_Oxidized,
      Control_UnoxidizedArea = Control_Unoxidized,
      Sample_OxidizedArea = Sample_Oxidized,
      Sample_UnoxidizedArea = Sample_Unoxidized
    )

  df_out$Control_OxidizedArea <- ifelse(df_out$Control_UnoxidizedArea > 0 & is.na(df_out$Control_OxidizedArea), 0, df_out$Control_OxidizedArea)
  df_out$Sample_UnoxidizedArea <- ifelse(df_out$Sample_OxidizedArea > 0 & is.na(df_out$Sample_UnoxidizedArea), 0, df_out$Sample_UnoxidizedArea)

  df_out$TotalSampleArea <- rowSums(df_out[, c("Sample_OxidizedArea", "Sample_UnoxidizedArea")], na.rm = TRUE)
  df_out$TotalControlArea <- rowSums(df_out[, c("Control_OxidizedArea", "Control_UnoxidizedArea")], na.rm = TRUE)

  df_out$EOMSample <- df_out$Sample_OxidizedArea / df_out$TotalSampleArea
  df_out$EOMControl <- df_out$Control_OxidizedArea / df_out$TotalControlArea
  df_out$EOM <- df_out$EOMSample - df_out$EOMControl

  df_out$N <- df_in %>%
    group_by(MasterProteinAccessions, Sequence) %>%
    count() %>%
    ungroup() %>%
    pull(n)

  # Calculate SD based on the specified formula
  sd_data <- df_in %>%
    group_by(MasterProteinAccessions, Sequence, MOD, SampleControl) %>%
    summarise(
      SumPrecAbund = sum(`Precursor Abundance`),
      VarPrecAbund = var(`Precursor Abundance`),
      .groups = 'drop'
    )

  # Debugging: Print the first few rows of sd_data
  print(head(sd_data))

  sd_calculation <- sd_data %>%
    pivot_wider(
      names_from = c(SampleControl, MOD),
      values_from = c(SumPrecAbund, VarPrecAbund),
      names_glue = "{.value}_{SampleControl}_{MOD}",
      values_fill = list(SumPrecAbund = 0, VarPrecAbund = 0)
    ) %>%
    rowwise() %>%
    mutate(
      SD_Oxidized = sqrt(
        (SumPrecAbund_Sample_Oxidized^2 * VarPrecAbund_Sample_Oxidized / SumPrecAbund_Sample_Oxidized^2) +
          (VarPrecAbund_Sample / SumPrecAbund_Sample^2) +
          (SumPrecAbund_Control_Oxidized^2 * VarPrecAbund_Control_Oxidized / SumPrecAbund_Control_Oxidized^2) +
          (VarPrecAbund_Control / SumPrecAbund_Control^2)
      ),
      SD = SD_Oxidized / (SumPrecAbund_Sample + SumPrecAbund_Control)
    ) %>%
    select(MasterProteinAccessions, Sequence, SD)

  # Debugging: Print the first few rows of sd_calculation
  print(head(sd_calculation))

  df_out <- df_out %>%
    left_join(sd_calculation, by = c("MasterProteinAccessions", "Sequence"))

  return(df_out)
}


testarea_calculations_pep(pd_data_fasta_merged)

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
