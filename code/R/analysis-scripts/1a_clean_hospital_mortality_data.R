# This script takes the hospital mortality data and ensures the place names are 
# formatted the same as the shapefiles, and summarises the data for the purpose 
# of fitting p, the probability of being HEV IgM +ve given fatal acute jaundice.

# it depends on:
##  - "code/dependencies.R"
##  - the hospital data at 
##    my_path "/hospital_surveillance/raw/hosp_cleanHEV_complete.csv"
##    my_path "/hospital_surveillance/raw/
##    HEV_Bangladesh_Hos surveillance_deidentified_30May2024.csv""

# it produces:
## - the reshaped set of data to fit p at 
##   my_path "/hospital_surveillance/clean/minimal_deaths_data.csv"
##   my_path "/hospital_surveillance/clean/minimal_deaths_data_coinfections.csv"

# load dependencies
source(here::here(
  "code/dependencies.R")
)

#--------------------------------------------------------------------------------

# read in hospital level data and update name formatting
deaths_hosp <- read.csv(here(
  paste0(my_path, "/hospital_surveillance/raw/hosp_cleanHEV_complete.csv"))
  )%>% # match upazila names to shape files for mapping later
  mutate(address_thana = case_when(grepl(address_thana, pattern = "Thana Total") ~ str_remove(address_thana, "\\s*Thana Total.*$"),
                                   address_thana == "Barguna Sadar Upazila Tota" ~ "Barguna Sadar",
                                   ((address_thana == "Brahmanpara")&(address_dist == "Comilla")) ~ "Brahman Para",
                                   address_thana == "Barisal Sadar" ~ "Barisal Sadar (Kotwali)",
                                   ((address_thana == "Daulat Khan") & (address_dist == "Bhola")) ~ "Daulatkhan",
                                   ((address_thana == "Chakoria") & (address_dist == "Cox'S Bazar")) ~ "Chakaria",
                                   address_thana == "Coxs Bazar Sadar" ~ "Cox'S Bazar Sadar",
                                   ((address_thana == "Saghata") & (address_dist == "Gaibandha")) ~ "Saghatta",
                                   ((address_thana == "Kotalipara") & (address_dist == "Gopalganj")) ~ "Kotali Para",
                                   ((address_thana == "Bagherpara") & (address_dist == "Jessore")) ~ "Bagher Para",
                                   ((address_thana == "Daulat Khan") & (address_dist == "Bhola")) ~ "Daulatkhan",
                                   ((address_thana == "Noakhali Sadar") & (address_dist == "Noakhali"))~ "Noakhali Sadar (Sudharam)",
                                   ((address_thana == "Nesarabad") & (address_dist == "Pirojpur"))~"Nesarabad (Swarupkati)",
                                   ((address_thana == "Kawkhali") & (address_dist == "Rangamati")) ~ "Kawkhali (Betbunia)",
                                   ((address_thana == "Golapganj") & (address_dist == "Sylhet"))~ "Golabganj",
                                   TRUE~ address_thana))

# add column for whether patient died 
deaths_hosp <- deaths_hosp %>%
  mutate(died = case_when((D_outcome == "death" | 
                             f1_pt_healthStatus == "deceased" | 
                             f2_pt_healthStatus == "deceased"| 
                             f3_pt_healthStatus == "deceased"|
                             f4_pt_healthStatus == "deceased")~"YES",
                          (D_outcome != "death" & 
                             f1_pt_healthStatus != "deceased" & 
                             f2_pt_healthStatus != "deceased"&
                             f3_pt_healthStatus != "deceased"&
                             f4_pt_healthStatus != "deceased")~"NO",
                          TRUE~ as.character(NA)))
table(deaths_hosp$died, useNA = "always")

# add column for maternal
deaths_hosp <- deaths_hosp %>%
  mutate(maternal = case_when(
    (preg_status == "pregnant" | preg_status == "recent_pregnant")~"YES",
    ((preg_status == "not_pregnant") | is.na(preg_status)) ~"NO",
    TRUE~as.character(NA)))

# add age_class
deaths_hosp <- deaths_hosp %>%
  mutate(age_class = cut(age, breaks = c(13, 44, 59, 90)))

# save minimal data set
d <- deaths_hosp %>% dplyr::select(BGD_ID, BGD_IgM_Result, 
                                   hospital, hospitalID, year_adm, 
                                   died, maternal, age_class, age, sex)

write.csv(d, paste0(my_path, "/hospital_surveillance/clean/minimal_deaths_data.csv"))


# what about coinfections?

## merge with HBV and HAV results

### read in HAV/HBV data
coinf <- read.csv(here(
  paste0(my_path, "/hospital_surveillance/raw/HEV_Bangladesh_Hos surveillance_deidentified_30May2024.csv"))
  )%>%
  dplyr::select(dataid, HAV, HBcIgM, HBsAg)%>%
  mutate(dataid = as.character(dataid))

d_co <- full_join(d, coinf, by = c("BGD_ID" = "dataid")) %>%
  filter(!is.na(BGD_IgM_Result))

### save - eventually use this instead of version without coinfections 
### but for now save both so it doesn't break everything
write.csv(d_co, here(
  paste0(my_path, "/hospital_surveillance/clean/minimal_deaths_data_coinfections.csv"))
)
