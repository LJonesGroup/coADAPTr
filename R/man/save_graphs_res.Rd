\name{save_graphs_res}
\alias{save_graphs_res}
\title{
  Save Graphs for Each Protein in Data Frame
}
\usage{
  save_graphs_res(df, file_output, excel_filename)
}
\arguments{
  \item{df}{
    A data frame containing protein information.
  }
  \item{file_output}{
    Character string specifying the directory to save the graphs.
  }
  \item{excel_filename}{
    Character string specifying the Excel filename to save the data.
  }
}
\description{
  This function saves graphs for each unique protein in the provided data frame.
}
\details{
  The function loops through each unique protein in the \code{MasterProteinAccessions} column of the data frame \code{df}, subsets the data frame for each protein, and generates a plot using the \code{generate_eom_plot_res} function. The plots are saved in the specified directory (\code{file_output}) and the data is saved in the specified Excel file (\code{excel_filename}).
}
\examples{
  \dontrun{
    # Example usage
    save_graphs_res(my_data, "plots/", "protein_data.xlsx")
  }
}
