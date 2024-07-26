#' Count the Number of Modified Peptides Per Protein
#'
#' @return A Table and Graph of the Number of Modified Peptides Per Protein
#' @export
#'
#' @examples count_peptides_per_protein()
#' @aliases count_peptides_per_protein
count_peptides_per_protein <- function() {
  data_list <- list()

  repeat {
    file_path <- file.choose()
    data <- read_excel(file_path)

    # Aggregate count data by protein accessions
    count_data <- data %>%
      group_by(MasterProteinAccessions, Condition) %>%
      summarise(SequenceCount = n(), .groups = 'drop')

    data_list <- append(data_list, list(count_data))

    more_files <- readline(prompt = "Do you want to input another file? (yes/no) ")
    if (tolower(more_files) != "yes") {
      break
    }
  }

  combined_data <- bind_rows(data_list)
  summary_data <- combined_data %>%
    group_by(Condition, SequenceCount) %>%
    summarise(ProteinCount = n(), .groups = 'drop')

  graph_title <- readline(prompt = "Enter the title for the Modified Peptides Per Protein Bar Graph: ")

  p <- ggplot(summary_data, aes(x = SequenceCount, y = ProteinCount, fill = Condition)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = graph_title, x = "Number of Peptides per Protein", y = "Number of Proteins") +
    theme_minimal()

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



