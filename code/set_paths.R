# paths to data on local machine

set_paths <- function(info = Sys.info()) {
  if (info["user"] == "adighe") {
    # paths for Amy's JHU laptop
    paths = list(
      local = "/Users/ADIGHE/Documents/HEV IFR/data"
    )
  } else if (info["user"] == "ad3716") {
    # paths for Amy's current laptop
    paths = list(
      local = "C:/Users/ad3716/OneDrive - Imperial College London/Documents/IFR/data_private"
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

