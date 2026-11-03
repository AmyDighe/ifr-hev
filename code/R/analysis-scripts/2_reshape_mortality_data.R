# This script takes the clean hospital mortality data and the clean community 
# mortality data and merges them, reshaping for model fitting

# it depends on:
##  - "code/dependencies.R"
##  - the clean community data at 
##    my_path,"/community_survey/clean/clean_community_decedant_level.rds"
##  - the clean hospital data at:
##    my_path "/hospital_surveillance/clean/minimal_deaths_data.csv"

# it produces:
## - the reformatted mortality data ready to inut into the IFR estimation 
##   framework, at:
##   my_path, "/model_input/deaths.rds"

# load dependencies
source("code/dependencies.R")

#-------------------------------------------------------------------------------

## read in hospitalisation data
hosp <- read.csv(file = paste0(my_path, "/hospital_surveillance/minimal_deaths_data.csv"))

## read in community data
comm <- readRDS(file = paste0(my_path,"/community_survey/clean_community_decedant_level.rds"))

## summarise each by MCH, age class and pregnancy status
hosp <- hosp %>% 
  mutate(MCH = recode(hospital,
                      CMCH = "Chittagong",
                      JIMCH = "Kishoreganj",
                      MAGOMCH = "Sylhet",
                      SBMCH = "Barisal",
                      SSMCH = "Dhaka",
                      SZRMCH = "Bogra"),
         age_class = recode(age_class, `(59,90]` = "(59,101]"))%>%
  group_by(MCH, age_class, maternal)%>%
  summarise(DAJ = sum(died == "YES", na.rm = T),
            DHEV = sum(died == "YES" & BGD_IgM_Result == "Positive", na.rm = T))%>%
  drop_na(MCH)

comm <- comm %>% 
  drop_na(dataid)%>%
  mutate(age_class = cut(age, breaks = c(13, 44, 59, 101)),
         MCH = recode(MCH,
                      `Chittagong (CMCH)` = "Chittagong",
                      `Kishoreganj (JIMCH)` = "Kishoreganj",
                      `Sylhet (SOMCH)` = "Sylhet",
                      `Barisal (SBMCH)` = "Barisal",
                      `Mitford (SSMCH)` = "Dhaka",
                      `Bogra (SZRMCH)` = "Bogra"),
         maternal = recode(maternal, `0` = "NO", `1` = "YES"))%>%
  group_by(MCH, age_class, maternal)%>%
  summarise(DcAJ = n())
  

## merge
deaths <- full_join(hosp, comm, by = c("MCH", "age_class", "maternal"))


## save
saveRDS(deaths, file = paste0(my_path, "/data_private/model_input/deaths.rds"))

# 
# ## scale up by population
# scalers <- readRDS("generated_data/popScalers.rds")
# 
# deaths <- left_join(deaths, scalers, by = "MCH")%>%
#   mutate(DHEV = if_else(DHEV == 0, 0.5, DHEV),
#          DcAJ = round(DcAJ*scaler, 0))
# saveRDS(deaths, file = paste0(my_path, 
#                               "/model_input/deaths_scaled_to_full_catchmentpop.rds"))

