#' Generate and Save Grouped Bar Graphs for Each Modified Residue Per Condition
#'
#' @return Grouped bar graphs for each modified residue in the data frame per condition
#' @export
#'
#' @examples generate_grouped_bar_plot_res()
#' @aliases generate_grouped_bar_plot_res
generate_grouped_bar_plot_res <- function() {
  # Select the Excel file
  cat("Select the Excel file containing the quantifiable residue level data: ")
  filepath <- file.choose()
  df_in <- readxl::read_excel(filepath)

  # Ask for color palette
  palettes <- list(
    "rainbow" = rainbow, "heat.colors" = heat.colors,
    "terrain.colors" = terrain.colors, "topo.colors" = topo.colors,
    "cm.colors" = cm.colors, "grayscale" = gray.colors,
    "blue_hues" = function(n) colorRampPalette(c("#084594", "#2171b5", "#6baed6", "#c6dbef"))(n),
    "viridis" = viridis::viridis, "plasma" = viridis::plasma, "inferno" = viridis::inferno
  )

  cat("Available color palettes:\n")
  for (i in seq_along(palettes)) cat(paste0(i, ": ", names(palettes)[i], "\n"))
  choice <- as.numeric(readline(prompt = "Select a color palette by number: "))
  if (is.na(choice) || choice < 1 || choice > length(palettes)) stop("Invalid palette choice.")
  condition_colors <- palettes[[choice]](length(unique(df_in$Condition)))

  # Ask for output filename
  cat("Enter the new filename for the residue level data that will be saved (without extension): ")
  filename <- trimws(readline())

  # --- AGGREGATE data by Protein, Condition, and Residue ---
  df_agg <- df_in %>%
    group_by(MasterProteinAccessions, Res, Condition) %>%
    summarise(
      EOM = mean(EOM, na.rm = TRUE),
      SD = sd(EOM, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(MasterProteinAccessions, Res)

  # Make Res a factor so ggplot keeps order
  df_agg$Res <- factor(df_agg$Res, levels = unique(df_agg$Res))

  # Create output dir
  graph_output_directory <- file.path(dirname(filepath), paste0(filename, "_ProteinResidueBarGraphs"))
  dir.create(graph_output_directory, showWarnings = FALSE, recursive = TRUE)

  # --- PLOT each protein ---
  for (protein in unique(df_agg$MasterProteinAccessions)) {
    temp <- subset(df_agg, MasterProteinAccessions == protein)

    fig <- ggplot(temp, aes(x = Res, y = EOM, fill = Condition)) +
      geom_bar(stat = "identity", position = position_dodge(width = 1)) +
      geom_errorbar(aes(ymin = EOM - SD, ymax = EOM + SD), width = 0.3,
                    position = position_dodge(width = 1)) +
      labs(title = paste(protein, "- Residue Level Modification"),
           x = "Residue", y = "Extent of Modification (EOM)", fill = "Condition") +
      theme_classic() +
      theme(
        text = element_text(size = 18),
        axis.text.x = element_text(angle = 45, hjust = 0.5),
        plot.title = element_text(size = 22, face = "bold")
      ) +
      scale_fill_manual(values = condition_colors)

    file_out <- file.path(graph_output_directory, paste0(protein, "_ResiduePlot.png"))
    ggsave(filename = file_out, plot = fig, device = "png", width = 12, height = 8, dpi = 300)

    cat("Saved plot for", protein, "at", file_out, "\n")
  }
}

