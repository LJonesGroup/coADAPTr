library(readxl)
library(openxlsx)
library(VennDiagram)
library(RColorBrewer)

venn_diagram <- function() {
  # Prompt user to select Excel file
  excel_file <- file.choose()

  # Read data from Excel file
  data <- read_excel(excel_file, col_names = FALSE)

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
  custom_palette <- brewer.pal(num_conditions, "Blues") # You can replace "Blues" with "Greys" for gray shades

  # Create Venn diagram plot
  venn_plot <- venn.diagram(
    x = condition_lists,
    category.names = venn_bubble_names,
    filename = NULL,
    col = custom_palette,
    fill = custom_palette,
    alpha = 0.5, # Adjust transparency for better visualization
    margin = 0.1, # Increase margin for better plot appearance
    fontfamily = "sans", # Set font family to Arial
    fontface = "bold",   # Set font face to bold
    cat.fontsize = 35,   # Set category font size to 35 (bubble titles)
    cex = 2.5,
    cat.fontfamily = "sans", # Set category title font family to sans
    cat.fontface = "bold"  # Set category title font face to bold
  )

  # Increase the font size of bubble titles (category names)
  venn_plot$vpList$category[[1]]$fontsize <- 35  # Adjust the font size as needed

  # Prompt user to specify the name of the output Excel file for data
  excel_file_name <- readline(prompt = "Enter the name of the output Excel file for data (without extension): ")
  excel_file_path <- file.path(output_folder, paste0(excel_file_name, ".xlsx"))

  # Create a new Excel workbook
  wb <- createWorkbook()

  # Add a worksheet for all data
  addWorksheet(wb, "Venn_Data")

  # Write data to different columns within the same worksheet
  start_col <- 1  # Start writing data from column 1
  for (i in 1:length(venn_bubble_names)) {
    # Calculate the end column
    end_col <- start_col + length(condition_lists[[i]]) - 1

    # Write data to the specified range of columns
    col_header <- venn_bubble_names[i]
    writeData(wb, sheet = "Venn_Data", x = data.frame(Condition = col_header, Data = condition_lists[[i]]),
              startCol = start_col, startRow = 1, colNames = TRUE)

    # Update the start column for the next set of data
    start_col <- end_col + 2  # Add some padding between sets of data
  }

  # Save the Excel workbook
  saveWorkbook(wb, excel_file_path)

  cat("Data from each Venn bubble and intersection saved to:", excel_file_path, "\n")

  # Prompt user to specify the name of the output PNG file
  png_file_name <- readline(prompt = "Enter the name of the output PNG file (without extension): ")
  png_file <- file.path(output_folder, paste0(png_file_name, ".png"))

  # Save plot as PNG
  png(filename = png_file, width = 800, height = 800) # Adjust width and height as needed
  grid.draw(venn_plot)
  dev.off()

  cat("Venn diagram plot saved at:", png_file, "\n")
}

# Call the function to create the Venn diagram plot and save data to Excel
venn_diagram()





