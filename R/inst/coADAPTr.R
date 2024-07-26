# Setting Up ----------------------------------------------------------------------------------------------------
#Before Analyzing your data be sure to run this function to install necessary
#Packages, resolve package conflicts, import required data, and select result
#output folders. Follow the console prompts and anticipate a File Explorer up.

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


# FPOP Calculations ---------------------------------------------------------------------------------------------
#Calculate the Extent of Modification (EOM) for each peptide and residue based on the LFQ data
#quant_graph_df containing all of the data that is acceptable for graphing
#Areas_pep and Areas_res contain the data that is used to calculate the EOM and will
#contain all data. This includes cases where the EOM could be negative (high background oxidation)
#or the SD is greater than the EOM (data has a high variance-likely due to experimental conditions)

FPOP_Calculations()

#Saving Tables and Plots------------------------------------------------------------------------------------------------------------
#Save data frames as Excel files and save grouped bar from the plots as PNG files

Tables_and_Graphs()
library(ggplot2)
library(readxl)
library(dplyr)

generate_grouped_bar_plot_res <- function() {
  # Prompt the user to select the Excel file
  cat("Select the Excel file containing the quantifiable residue level data: ")
  filepath <- file.choose()

  # Load the data from the selected Excel file
  df_in <- readxl::read_excel(filepath)

  # Arrange the dataframe by start
  df_in <- df_in %>%
    arrange(start)

  # Auto-detect unique conditions
  unique_conditions <- unique(df_in$Condition)

  # Manually specify blue colors for filling the bars
  condition_colors <- c("#0570b0", "#3690c0", "#74a9cf", "#a6bddb", "#d0d1e6", "#084594", "#2171b5", "#4292c6", "#6baed6", "#9ecae1", "#c6dbef", "#08519c", "#3182bd", "#6baed6", "#9ecae1", "#c6dbef")

  # Create a factor variable to represent the sorted order of conditions
  df_in$Condition <- factor(df_in$Condition, levels = unique_conditions)

  # Prompt the user to specify the filename
  cat("Enter the filename for the residue level data (without extension): ")
  filename <- readline()

  # Remove any leading or trailing whitespace
  filename <- trimws(filename)

  # Iterate over each protein and make a grouped bar plot for it
  for (protein in unique(df_in$MasterProteinAccessions)) {
    # Subset the dataframe for this protein
    temp <- subset(df_in, MasterProteinAccessions == protein)

    # Generate a grouped bar plot of the extent of modification for each residue
    # that maps to this protein, with different conditions represented by color
    fig <- ggplot(temp, aes(x = Res, y = EOM, fill = Condition)) +
      geom_bar(position = position_dodge(width = 0.9), stat = "identity") +
      geom_errorbar(aes(ymin = EOM - SD, ymax = EOM + SD), position = position_dodge(width = 0.9), width = 0.25) +
      labs(title = paste(protein, "Residue Level Analysis"),
           x = "Residue",
           y = "Extent of Modification",
           fill = "Condition") +
      theme_classic() +
      theme(text = element_text(size = 18, family = "sans"),
            plot.title = element_text(hjust = 0.5, size = 20, family = "sans", face = "bold"),
            legend.text = element_text(size = 16, family = "sans"),
            legend.title = element_text(size = 18, family = "sans"),
            axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
            axis.title.x = element_text(margin = margin(t = 20))) +
      scale_fill_manual(values = condition_colors)

    # Create the output directory for bar graphs based on the file output and excel filename
    graph_output_directory <- file.path(dirname(filepath), paste0(filename, "_ResidueLevelGroupedBarGraphs"))
    dir.create(graph_output_directory, showWarnings = FALSE, recursive = TRUE)

    # Generate the full file path for this protein and save the figure
    file_out <- file.path(graph_output_directory, paste0(protein, ".png"))
    ggsave(filename = file_out, plot = fig, device = "png", width = 10, height = 8, dpi = 300)

    # Print a message to indicate successful saving
    cat("Grouped bar graph for", protein, "saved as", file_out, "\n")
  }
}

generate_grouped_bar_plot_res()

#Step 17 Saving Venn Diagrams
#Generate Venn Diagrams
venn_diagram <- function() {
  # Prompt user to select Excel file
  excel_file <- file.choose()

  # Read data from Excel file
  data <- read_excel(excel_file, col_names = TRUE)

  # Prompt user to select output folder
  output_folder <- choose.dir()

  # Get the number of conditions
  num_conditions <- ncol(data)

  # Get bubble names from the user
  venn_bubble_names <- readline(prompt = "Enter names for the Venn diagram bubbles separated by commas: ")
  venn_bubble_names <- strsplit(venn_bubble_names, ",")[[1]]

  # Ensure the number of bubble names matches the number of conditions
  if (length(venn_bubble_names) != num_conditions) {
    stop("Number of bubble names provided does not match the number of columns in the Excel file.")
  }

  # Process each column separately
  condition_lists <- lapply(seq_len(num_conditions), function(i) {
    column_data <- na.omit(data[[i]])
    if (length(column_data) == 0) {
      return(NULL)
    } else {
      return(column_data)
    }
  })

  # Generate a custom color palette (using shades of blue or gray)
  custom_palette <- c("#0570b0", "#3690c0", "#74a9cf", "#a6bddb", "#d0d1e6", "#084594", "#2171b5", "#4292c6", "#6baed6", "#9ecae1", "#c6dbef", "#08519c", "#3182bd", "#6baed6", "#9ecae1", "#c6dbef")
  custom_palette <- custom_palette[1:num_conditions] # ensuring the palette has the desired number of colors

  # Create Venn diagram plot
  venn_plot <- venn.diagram(
    x = condition_lists,
    category.names = venn_bubble_names,
    filename = NULL,
    col = custom_palette,
    fill = custom_palette,
    alpha = 0.5, # Adjust transparency for better visualization
    margin = 0.1, # Increase margin for better plot appearance
    fontfamily = "sans", # Set font family to Helvetica (or similar font)
    fontface = "bold",   # Set font face to bold
    cat.fontsize = 20,   # Set category font size to 20
    cex = 2.5,           # Adjust overall font size
    cat.cex = 2.5,         # Set category title font size to 24
    cat.fontfamily = "sans", # Set category title font family
    cat.fontface = "bold"  # Set category title font face to bold
  )

  # Prompt user to specify the name of the output PNG file
  png_file_name <- readline(prompt = "Enter the name of the output PNG file (without extension): ")
  png_file <- file.path(output_folder, paste0(png_file_name, ".png"))

  # Save plot as PNG
  png(filename = png_file, width = 1000, height = 1000) # Adjust width and height as needed
  grid.draw(venn_plot)
  dev.off()

  cat("Venn diagram plot saved at:", png_file, "\n")

  # Calculate overlap
  overlap <- calculate.overlap(condition_lists)

  # Save overlap to Excel
  overlap_file <- file.path(output_folder, paste0(png_file_name, "_Overlap.xlsx"))
  write.xlsx(overlap, overlap_file)

  cat("Overlap information saved at:", overlap_file, "\n")

  if (FALSE) {
    venn_diagram()
  }
}
venn_diagram()


#Step 18 Saving Modifications Per Peptide and Residue

#Counting modified peptides per protein and graphing

count_peptides_per_protein <- function() {
  data_list <- list()
  conditions <- c()

  repeat {
    file_path <- file.choose()
    condition <- readline(prompt = "What condition does this data correspond to? ")
    conditions <- c(conditions, condition)

    data <- read_excel(file_path)

    count_data <- data %>%
      group_by(MasterProteinAccessions) %>%
      summarise(SequenceCount = n()) %>%
      ungroup()

    data_list[[condition]] <- count_data

    more_files <- readline(prompt = "Do you want to input another file? (yes/no) ")
    if (tolower(more_files) != "yes") {
      break
    }
  }

  graph_title <- readline(prompt = "Enter the title for the Modified Peptides Per Protein Bar Graph: ")

  combined_data <- bind_rows(data_list, .id = "Condition")
  summary_data <- combined_data %>%
    group_by(Condition, SequenceCount) %>%
    summarise(ProteinCount = n(), .groups = 'drop')

  p <- ggplot(summary_data, aes(x = SequenceCount, y = ProteinCount, fill = Condition)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = graph_title, x = "Number of Peptides per Protein", y = "Number of Proteins") +
    theme_minimal(base_size = 15) +
    theme(panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA),
          legend.background = element_rect(fill = "white", colour = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 20))

  print(p)

  output_dir <- choose.dir(default = "", caption = "Select directory to save the count data and graph")
  file_name <- readline(prompt = "Enter the base name for the modified peptides per protein files (without extension): ")

  output_file <- file.path(output_dir, paste0(file_name, ".xlsx"))
  graph_file <- file.path(output_dir, paste0(file_name, ".png"))

  write_xlsx(summary_data, output_file)
  ggsave(graph_file, plot = p, width = 10, height = 6)

  print(paste("Data saved to", output_file))
  print(paste("Graph saved to", graph_file))
}

# To run the function
count_peptides_per_protein()
####


#Count Residues Per Protein
count_residue_entries_per_protein <- function() {
  data_list <- list()
  conditions <- c()

  repeat {
    file_path <- file.choose()
    condition <- readline(prompt = "What condition does this data correspond to? ")
    conditions <- c(conditions, condition)

    data <- read_excel(file_path)

    # Ensure "Res" column contains non-NA values and only valid entries with numbers followed by letters
    data <- data %>% filter(!is.na(Res) & grepl("^[A-Za-z][0-9]+", Res))

    count_data <- data %>%
      group_by(MasterProteinAccessions) %>%
      summarise(ResidueEntryCount = n()) %>%
      ungroup()

    data_list[[condition]] <- count_data

    more_files <- readline(prompt = "Do you want to input another file? (yes/no) ")
    if (tolower(more_files) != "yes") {
      break
    }
  }

  graph_title <- readline(prompt = "Enter the title for the Modified Residues Per Protein Bar Graph: ")

  combined_data <- bind_rows(data_list, .id = "Condition")
  summary_data <- combined_data %>%
    group_by(Condition, ResidueEntryCount) %>%
    summarise(ProteinCount = n(), .groups = 'drop')

  p <- ggplot(summary_data, aes(x = ResidueEntryCount, y = ProteinCount, fill = Condition)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = graph_title, x = "Number of Residue Entries per Protein", y = "Number of Proteins") +
    theme_minimal(base_size = 15) +
    theme(panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA),
          legend.background = element_rect(fill = "white", colour = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 20))

  print(p)

  output_dir <- choose.dir(default = "", caption = "Select directory to save the count data and graph")
  file_name <- readline(prompt = "Enter the base name for modified residues per protein files (without extension): ")

  output_file <- file.path(output_dir, paste0(file_name, ".xlsx"))
  graph_file <- file.path(output_dir, paste0(file_name, ".png"))

  write_xlsx(summary_data, output_file)
  ggsave(graph_file, plot = p, width = 10, height = 6)

  print(paste("Data saved to", output_file))
  print(paste("Graph saved to", graph_file))
}

# To run the function
count_residue_entries_per_protein()


dev.off()
