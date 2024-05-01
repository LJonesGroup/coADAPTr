pkgname <- "coADAPTr"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('coADAPTr')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("FASTA_file")
### * FASTA_file

flush(stderr()); flush(stdout())

### Name: FASTA_file
### Title: Import FASTA File
### Aliases: FASTA

### ** Examples

  # Call the function to read a FASTA file
  FASTA <- FASTA_file()




cleanEx()
nameEx("annotate_features")
### * annotate_features

flush(stderr()); flush(stdout())

### Name: annotate_features
### Title: Annnotate
### Aliases: annotate

### ** Examples

  pd_data_annotated <- annotate_features(pd_data)



cleanEx()
nameEx("area_calculations_pep")
### * area_calculations_pep

flush(stderr()); flush(stdout())

### Name: area_calculations_pep
### Title: calculate the peptide total areas
### Aliases: area_pep

### ** Examples

  Areas_pep <- area_calculations_pep(pd_data_fasta_merged)



cleanEx()
nameEx("area_calculations_res")
### * area_calculations_res

flush(stderr()); flush(stdout())

### Name: area_calculations_res
### Title: Calculate Area and Associated Metrics
### Aliases: area_res

### ** Examples

  # Call the function with a data frame pd_data_fasta_merged
  Areas_res <- area_calculations_res(pd_data_fasta_merged)



cleanEx()
nameEx("filtered_graphing_df_pep")
### * filtered_graphing_df_pep

flush(stderr()); flush(stdout())

### Name: filtered_graphing_df_pep
### Title: filter the quantifiable peptide extent of modifications
### Aliases: filter_pep

### ** Examples

 quant_graph_df_pep <- filtered_graphing_df_pep(graphing_df_pep)



cleanEx()
nameEx("filtered_graphing_df_res")
### * filtered_graphing_df_res

flush(stderr()); flush(stdout())

### Name: filtered_graphing_df_res
### Title: Filter Graphing Data
### Aliases: filtered_graphing_df_res

### ** Examples

  # Call the function to filter graphing data
  quant_graph_df_res <- filtered_graphing_df_res(graphing_df_res)
  print(quant_graph_df_res)



cleanEx()
nameEx("generate_eom_plot_pep")
### * generate_eom_plot_pep

flush(stderr()); flush(stdout())

### Name: generate_peptide_bar_graph
### Title: Generate Peptide Bar Graph
### Aliases: generate_peptide_bar_graph

### ** Examples

  # Call the function to generate peptide bar graphs
  generate_peptide_bar_graph(df_in, file_output, excel_filename)



cleanEx()
nameEx("generate_eom_plot_res")
### * generate_eom_plot_res

flush(stderr()); flush(stdout())

### Name: generate_protein_eom_plots
### Title: Generate Extent of Modification (EOM) Plots for Proteins
### Aliases: generate_protein_eom_plots

### ** Examples

  # Define the excel_filename variable
  excel_filename <- tools::file_path_sans_ext(basename(file_path))

  # Call the function to generate EOM plots for proteins
  generate_protein_eom_plots(quant_graph_df_pep, file_output, excel_filename)



cleanEx()
nameEx("grab_seq_metadata_pep")
### * grab_seq_metadata_pep

flush(stderr()); flush(stdout())

### Name: grab_seq_metadata_pep
### Title: grab the quantifiable peptide extent of modifications
### Aliases: grab_pep

### ** Examples

  graphing_df_pep <- Areas_pep 
  left_join(grab_seq_metadata_pep(pd_data_fasta_merged))


graphing_df_pep <- graphing_df_pep[order(graphing_df_pep$start), ]
graphing_df_pep$MasterProteinAccessions <- gsub(".*\|(.*?)\|.*", "\1", graphing_df_pep$MasterProteinAccessions)



cleanEx()
nameEx("grab_seq_metadata_res")
### * grab_seq_metadata_res

flush(stderr()); flush(stdout())

### Name: grab_seq_metadata_res
### Title: Grab the quantifiable residue extent of modifications
### Aliases: grab_res

### ** Examples

  graphing_df_pep <- Areas_pep 
  left_join(grab_seq_metadata_pep(pd_data_fasta_merged))


graphing_df_pep <- graphing_df_pep[order(graphing_df_pep$start), ]
graphing_df_pep$MasterProteinAccessions <- gsub(".*\|(.*?)\|.*", "\1",
graphing_df_pep$MasterProteinAccessions)



cleanEx()
nameEx("graphing_data_res")
### * graphing_data_res

flush(stderr()); flush(stdout())

### Name: graphing_data_res
### Title: Prepare Data for Graphing
### Aliases: graphing_df_res

### ** Examples

  # Call the function to prepare data for graphing
  graphing_df <- graphing_data_res(Areas_res)
  print(graphing_df)



cleanEx()
nameEx("locate_startend_res")
### * locate_startend_res

flush(stderr()); flush(stdout())

### Name: locate_startend_res
### Title: Locate start and end residues
### Aliases: locate

### ** Examples

  pd_data_fasta_merged <- locate_startend_res(pd_data_annotated)



cleanEx()
nameEx("output_folder")
### * output_folder

flush(stderr()); flush(stdout())

### Name: output_folder
### Title: Select a folder for your tables and figures to be saved
### Aliases: output

### ** Examples

  # Call the function to select a folder
  selected_path <- output_folder()




cleanEx()
nameEx("parse_fasta")
### * parse_fasta

flush(stderr()); flush(stdout())

### Name: parse_fasta
### Title: Parse fasta files
### Aliases: parse

### ** Examples

  FASTA <- parse_fasta(FASTA)



cleanEx()
nameEx("process_data")
### * process_data

flush(stderr()); flush(stdout())

### Name: process_data
### Title: Select the file location of sequence searched MS sata
### Aliases: file_path

### ** Examples

  # Call the function to process sequence searched MS data
  processed_data <- process_data()
  file<- process_data()



cleanEx()
nameEx("save_data_frames")
### * save_data_frames

flush(stderr()); flush(stdout())

### Name: save_data_frames
### Title: Save Data Frames as Separate Excel Files
### Aliases: save_df

### ** Examples

  # Example usage:
  # Define data frames to be saved
  TotalsTable <- data.frame()
  quant_graph_df_pep <- data.frame()
  quant_graph_df_res <- data.frame()
  graphing_df_pep <- data.frame()
  graphing_df_res <- data.frame()

  # Save the data frames as separate Excel files
  save_data_frames(file_path, file_output, TotalsTable = TotalsTable, quant_graph_df_pep = quant_graph_df_pep, quant_graph_df_res = quant_graph_df_res, graphing_df_pep = graphing_df_pep, graphing_df_res = graphing_df_res)



cleanEx()
nameEx("save_graphs_pep")
### * save_graphs_pep

flush(stderr()); flush(stdout())

### Name: save_graphs_pep
### Title: Save Extent of Modification (EOM) Plots for modified residues
### Aliases: save_pep_graphs

### ** Examples

 save_graphs_pep(df_in_pep, file_output, excel_filename)



cleanEx()
nameEx("save_graphs_res")
### * save_graphs_res

flush(stderr()); flush(stdout())

### Name: save_graphs_res
### Title: Save Graphs for Each Protein in Data Frame
### Aliases: save_res_graphs

### ** Examples

  ## Not run: 
##D     # Example usage
##D     save_graphs_res(my_data, "plots/", "protein_data.xlsx")
##D   
## End(Not run)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
