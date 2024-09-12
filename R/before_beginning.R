#' Before Analyzing Your Data be Sure to Run this Function
#'
#' @return Required  Packaged and Conflict Preferences
#' @export
#'
#' @examples before_beginning()
#' @aliases before_beginning
before_beginning <- function() {
  options(stringsAsFactors = FALSE)

  required_packages <- c("ExcelFunctionsR", "remotes", "data.table", "stringr", "protr",
                         "plyr","viridis", "extrafont", "readxl", "ggplot2", "eulerr", "grid",
                         "tidyverse", "EnvStats", "dplyr", "writexl", "conflicted", "tcltk",
                         "phylotools", "parallel", "rlist","argparser", "Cairo", "VennDiagram")


  package.check <- lapply(
    required_packages,
    FUN = function(x) {
      if (!require(x, character.only = TRUE)) {
        install.packages(x, dependencies = TRUE)
        library(x, character.only = TRUE)
      }
    }
  )


  setwd(getwd())


  conflict_prefer("arrange", winner = "dplyr")
  conflict_prefer("filter", winner = "dplyr")
  conflict_prefer("between", winner = "dplyr")
  conflict_prefer("compact", winner = "purrr")
  conflict_prefer("count", winner = "dplyr")
  conflict_prefer("desc", winner = "dplyr")
  conflict_prefer("failwith", winner = "dplyr")
  conflict_prefer("filter", winner = "dplyr")
  conflict_prefer("first", winner = "dplyr")
  conflict_prefer("hour", winner = "lubridate")
  conflict_prefer("id", winner = "dplyr")
  conflict_prefer("isoweek", winner = "lubridate")
  conflict_prefer("lag", winner = "dplyr")
  conflict_prefer("last", winner = "dplyr")
  conflict_prefer("mday", winner = "lubridate")
  conflict_prefer("minute", winner = "lubridate")
  conflict_prefer("month", winner = "lubridate")
  conflict_prefer("mutate", winner = "dplyr")
  conflict_prefer("quarter", winner = "lubridate")
  conflict_prefer("rename", winner = "dplyr")
  conflict_prefer("second", winner = "lubridate")
  conflict_prefer("summarise", winner = "dplyr")
  conflict_prefer("summarize", winner = "dplyr")
  conflict_prefer("transpose", winner = "purrr")
  conflict_prefer("wday", winner = "lubridate")
  conflict_prefer("week", winner = "lubridate")
  conflict_prefer("yday", winner = "lubridate")
  conflict_prefer("year", winner = "lubridate")
  conflict_prefer("read.fasta", winner = "phylotools")
  conflict_prefer("rename", winner = "dplyr")
}
