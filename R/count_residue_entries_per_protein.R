
#' Count the Modified Residues Per Protein
#'
#' @return An excel file of modified residues per protein and a bar graph of the data
#' @export
#'
#' @examples count_residue_entries_per_protein()
#' @aliases count_residue_entries_per_protein
count_residue_entries_per_protein <- function() {
  data_list <- list()

  repeat {
    message("Please select the file corresponding to quantifiable residue level data:")
    file_path <- file.choose()
    data <- read_excel(file_path)

    # Filter data and count residues per protein
    data <- data %>%
      filter(!is.na(Res) & grepl("^[A-Za-z][0-9]+", Res)) %>%
      group_by(MasterProteinAccessions, Condition) %>%
      summarise(ResidueEntryCount = n(), .groups = 'drop')

    data_list <- append(data_list, list(data))

    more_files <- readline(prompt = "Do you want to input another file? (yes/no) ")
    if (tolower(more_files) != "yes") {
      break
    }
  }

  combined_data <- bind_rows(data_list)
  summary_data <- combined_data %>%
    group_by(Condition, ResidueEntryCount) %>%
    summarise(ProteinCount = n(), .groups = 'drop')

  graph_title <- readline(prompt = "Enter the title for the Modified Residues Per Protein Bar Graph: ")

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

  output_dir <- choose.dir(default = "", caption = "Select directory/folder to save the count data and graph")
  file_name <- readline(prompt = "Enter the base name for modified residues per protein files (without extension): ")

  output_file <- file.path(output_dir, paste0(file_name, ".xlsx"))
  graph_file <- file.path(output_dir, paste0(file_name, ".png"))

  write_xlsx(summary_data, output_file)
  ggsave(graph_file, plot = p, width = 10, height = 6)

  print(paste("Data saved to", output_file))
  print(paste("Graph saved to", graph_file))
}
