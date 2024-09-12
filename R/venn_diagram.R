#' Create and Save a Venn Diagram Plot (Step 17)
#'
#' @return A Venn diagram plot of the modified protein lists per condition saved as a PNG file
#' @export
#'
#' @examples venn_diagram()
#' @aliases venn_diagram
venn_diagram <- function() {
  # Prompt user to select Excel file
  readline(prompt = "Please select the file containing the lists of modified proteins per condition (one condition per column):")
  excel_file <- file.choose()

  # Read data from Excel file
  readline(prompt = "Reading data from Excel file...")
  data <- read_excel(excel_file, col_names = TRUE)

  # Prompt user to select output folder
  readline(prompt = "Please select the output folder for the Venn diagram plot and overlap information to be saved:")
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

  # Calculate an appropriate cat.dist value
  cat.dist.value <- 0.1 # Adjust this value if necessary for proper distance

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
    cat.fontsize = 30,   # Set category font size to 20
    cex = 3,           # Adjust overall font size
    cat.cex = 3,         # Set category title font size to 24
    cat.fontfamily = "sans", # Set category title font family
    cat.fontface = "bold",  # Set category title font face to bold
    cat.dist = rep(cat.dist.value, num_conditions) # Set consistent distance for category titles
  )

  # Prompt user to specify the name of the output PNG file
  png_file_name <- readline(prompt = "Enter the name of the output PNG file (without extension): ")
  png_file <- file.path(output_folder, paste0(png_file_name, ".png"))

  # Save plot as PNG
  png(filename = png_file, width = 1500, height = 1500) # Adjust width and height as needed
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
