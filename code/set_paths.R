# paths to data on local machine


set_paths <- function(info = Sys.info()) {
  if (info["user"] == "adighe") {
    # paths for Amy's laptop
    paths = list(
      local = "/Users/ADIGHE/Documents/HEV IFR/data"
    )
  } else {
    # generic paths, if the machine is not above
    paths = list(
      local = "~/Users/"
    )
    warning("User and computer not recognized, please edit R/set_paths.R.")
  }
  return(paths)
}

