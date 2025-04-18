#' Generate and Save Grouped Bar Graphs for Each Modified Peptide Per Condition
#'
#' @return Grouped bar graphs for each modified peptide in the data frame per condition
#' @export
#'
#' @examples generate_grouped_bar_plot_pep()
#' @aliases generate_grouped_bar_plot_pep
generate_grouped_bar_plot_pep <- function() {
  cat("Select the Excel file containing the quantifiable peptide level data: ")
  filepath <- file.choose()
  df_in <- readxl::read_excel(filepath)

  palettes <- list(
    "rainbow" = rainbow,
    "heat.colors" = heat.colors,
    "terrain.colors" = terrain.colors,
    "topo.colors" = topo.colors,
    "cm.colors" = cm.colors,
    "grayscale" = gray.colors,
    "blue_hues" = function(n) colorRampPalette(c("#084594", "#2171b5", "#6baed6", "#c6dbef"))(n),
    "viridis" = viridis::viridis,
    "plasma" = viridis::plasma,
    "inferno" = viridis::inferno
  )

  show_palette_preview <- function(palettes, num_colors = 5) {
    num_palettes <- length(palettes)
    par(mfrow = c(ceiling(num_palettes / 2), 2), mar = c(1, 1, 2, 1))
    for (i in seq_along(palettes)) {
      colors <- palettes[[i]](num_colors)
      plot(1:num_colors, pch = 15, cex = 3, col = colors, xlab = "", ylab = "",
           xaxt = 'n', yaxt = 'n', bty = 'n', main = paste("Palette:", names(palettes)[i]))
    }
    par(mfrow = c(1, 1))
  }

  cat("Would you like to see a preview of the color palettes before building your plot? (yes/no): ")
  show_preview <- tolower(readline())
  if (show_preview == "yes") {
    show_palette_preview(palettes)
  }

  cat("Available color palettes:\n")
  for (i in seq_along(palettes)) {
    cat(paste0(i, ": ", names(palettes)[i], "\n"))
  }

  choice <- as.numeric(readline(prompt = "Select a color palette by number: "))
  if (is.na(choice) || choice < 1 || choice > length(palettes)) {
    stop("Invalid selection. Please run the function again and select a valid option.")
  }

  condition_colors <- palettes[[choice]](length(unique(df_in$Condition)))

  cat("Enter the new filename for the peptide level data that will be saved (without extension): ")
  filename <- trimws(readline())

  df_in <- df_in %>%
    arrange(start)

  df_in$peptide_plot <- make.unique(as.character(df_in$peptide))
  df_in$peptide_plot <- factor(df_in$peptide_plot, levels = df_in$peptide_plot)

  for (protein in unique(df_in$MasterProteinAccessions)) {
    temp <- subset(df_in, MasterProteinAccessions == protein)

    fig <- ggplot(temp, aes(x = peptide_plot, y = EOM, fill = Condition)) +
      geom_bar(position = position_dodge(width = 1), stat = "identity") +
      geom_errorbar(aes(ymin = EOM - SD, ymax = EOM + SD), position = position_dodge(width = 1), width = 0.4) +
      labs(title = paste(protein, "Peptide Level Analysis"),
           x = "Peptide",
           y = "Extent of Modification",
           fill = "Condition") +
      theme_classic() +
      theme(
        text = element_text(size = 20, family = "sans"),
        plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
        legend.text = element_text(size = 18, family = "sans"),
        legend.title = element_text(size = 20, family = "sans"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
        axis.title.x = element_text(margin = margin(t = 22))
      ) +
      scale_fill_manual(values = condition_colors)

    graph_output_directory <- file.path(dirname(filepath), paste0(filename, "_PeptideLevelGroupedBarGraphs"))
    dir.create(graph_output_directory, showWarnings = FALSE, recursive = TRUE)

    file_out <- file.path(graph_output_directory, paste0(protein, ".png"))
    ggsave(filename = file_out, plot = fig, device = "png", width = 12, height = 10, dpi = 1200)

    cat("Grouped bar graph for", protein, "saved as", file_out, "\n")
  }
}

