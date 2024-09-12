venn_diagram <- function() {

  # List of 10 available color palettes, including grayscale and blue-hued palette
  palettes <- list(
    "rainbow" = rainbow,
    "heat.colors" = heat.colors,
    "terrain.colors" = terrain.colors,
    "topo.colors" = topo.colors,
    "cm.colors" = cm.colors,
    "grayscale" = gray.colors,  # Grayscale palette
    "blue_hues" = function(n) colorRampPalette(c("#084594", "#2171b5", "#6baed6", "#c6dbef"))(n),  # Blue hues
    "viridis" = viridis::viridis,  # Viridis color palette
    "plasma" = viridis::plasma,    # Plasma color palette
    "inferno" = viridis::inferno   # Inferno color palette
  )

  # Define the function to display the color palettes in a pop-up
  show_palette_popup <- function() {
    # Create a window
    win <- tktoplevel()
    tkwm.title(win, "Select a Color Palette")

    # Add instructions
    tkgrid(tklabel(win, text = "Select a Color Palette for your Venn Diagram:"))

    # Loop over each palette and display a preview of the colors
    for (i in seq_along(palettes)) {
      palette_name <- names(palettes)[i]
      palette_colors <- palettes   # Generate 5 colors for preview

      # Create a row with the palette name
      color_row <- tkframe(win)
      tkgrid(tklabel(color_row, text = palette_name))

      # Display color swatches as individual labels with colored background
      for (color in palette_colors) {
        color_box <- tklabel(color_row, width = 2, relief = "raised", background = color)
        tkgrid(color_box, padx = 2)
      }

      # Add the row to the window
      tkgrid(color_row)

      # Add a button to select this palette
      select_button <- tkbutton(win, text = "Select", command = function() {
        assign("selected_palette", palettes[[i]], envir = .GlobalEnv)
        assign("palette_name", palette_name, envir = .GlobalEnv)
        tkdestroy(win)  # Close the window when selected
      })
      tkgrid(select_button)
    }

    # Wait for the user to select a palette
    tkwait.window(win)
  }

  # Display the pop-up for palette selection
  show_palette_popup()

  # If the user did not select a palette, stop the function
  if (!exists("selected_palette", envir = .GlobalEnv)) {
    stop("No color palette selected. Please run the function again.")
  }

  # Prompt user to select Excel file
  cat("Please select the file containing the lists of modified proteins per condition (one condition per column):\n")
  excel_file <- file.choose()

  # Read data from Excel file
  cat("Reading data from Excel file...\n")
  data <- read_excel(excel_file, col_names = TRUE)

  # Prompt user to select output folder
  cat("Please select the output folder for the Venn diagram plot and overlap information to be saved:\n")
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

  # Generate a custom color palette based on the user's selection
  custom_palette <- selected_palette(num_conditions)

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
    cat.fontsize = 40,   # Set category font size to 20
    cex = 4,             # Adjust overall font size
    cat.cex = 5,         # Set category title font size to 24
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

  # Clean up the global environment
  rm(selected_palette, palette_name, envir = .GlobalEnv)

  if (FALSE) {
    venn_diagram()
  }
}
