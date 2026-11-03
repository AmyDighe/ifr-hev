# loading libraries and internal functions needed

###############
## libraries ##
###############

library(tidyverse)
library(cowplot)
library(ggExtra)
library(purrr)
library(ggridges)
library(INLA)
library(sf)
library(raster)
library(here)
library(purrr)
library(rstan)
library(survival) # for poisson CIs
library(posterior)

source("code/set_paths.R")
my_path <- set_paths()

########################
## internal functions ##
########################

## projects sf object into transverse mercator coordinates reference system
## used for Bangladesh

transform_to_btm <- function(my_obj){
  
  if ("sf" %in% class(my_obj)){
    rc <- st_transform(my_obj,
                       crs="+proj=tmerc +lat_0=0 +lon_0=90 +k=0.9996 +x_0=500000 +y_0=0 +a=6377276.345 +b=6356075.41314024 +towgs84=283.7,735.9,261.1,0,0,0,0 +units=m +no_defs")
  } else if(class(my_obj) == "RasterLayer"){
    rc <- raster::projectRaster(my_obj,
                                crs="+proj=tmerc +lat_0=0 +lon_0=90 +k=0.9996 +x_0=500000 +y_0=0 +a=6377276.345 +b=6356075.41314024 +towgs84=283.7,735.9,261.1,0,0,0,0 +units=m +no_defs")
  } else {
    stop("can only transform rasterlayers and sf objects at the moment")
  }
  
}

## loads GADM data for Bangladesh from saved shapefiles
##' could have used raster::getData('GADM',country="BGD",level=level)
##' but this seems more robust to any changes in GADM
##' @title
##' @return
##' @author Andrew A

load_bgd_shp <- function(level=2){
  
  if(level==4){
    rc <- readRDS(here("public_data", "shapefiles", "2018", "BGD_adm4.rds")) %>% sf::st_as_sf() %>% mutate(adm4name = NAME_4)
  } else if(level == 3){
    rc <- readRDS(here("public_data", "shapefiles", "2018", "BGD_adm3.rds")) %>% sf::st_as_sf() %>% mutate(adm3name = NAME_3)
  } else if(level == 2){
    rc <- readRDS(here("public_data", "shapefiles", "2018", "BGD_adm2.rds")) %>% sf::st_as_sf() %>% mutate(adm2name = NAME_2)
  } else if(level == 1){
    rc <- readRDS(here("public_data", "shapefiles", "2018", "BGD_adm1.rds"))%>% sf::st_as_sf() %>% mutate(adm1name = NAME_1)
  } else if(level == 0){
    rc <- readRDS(here("public_data", "shapefiles", "2018", "BGD_adm0.rds")) %>% sf::st_as_sf()
  }
  return(rc)
}

