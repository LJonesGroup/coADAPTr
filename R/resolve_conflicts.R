resolve_conflicts <- function() {
  preferred_packages <- list(
    arrange = "dplyr",
    between = "dplyr",
    compact = "purrr",
    count = "dplyr",
    desc = "dplyr",
    failwith = "dplyr",
    filter = "dplyr",
    first = "dplyr",
    hour = "lubridate",
    id = "dplyr",
    isoweek = "lubridate",
    lag = "dplyr",
    last = "dplyr",
    mday = "lubridate",
    minute = "lubridate",
    month = "lubridate",
    mutate = "dplyr",
    quarter = "lubridate",
    rename = "dplyr",
    second = "lubridate",
    summarise = "dplyr",
    summarize = "dplyr",
    transpose = "purrr",
    wday = "lubridate",
    week = "lubridate",
    yday = "lubridate",
    year = "lubridate",
    read.fasta = "phylotools"
  )

  resolve <- function(conflict_list) {
    resolved_conflicts <- list()

    for (conflict in conflict_list) {
      if (conflict %in% names(preferred_packages)) {
        resolved_conflicts[[conflict]] <- preferred_packages[[conflict]]
      } else {
        resolved_conflicts[[conflict]] <- "No preference"
      }
    }

    return(resolved_conflicts)
  }

  return(resolve)
}
