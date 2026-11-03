# This script takes the community mortality data and fixes cluster/hospital 
# mismatches, merges with age and date of death variables, and cluster info.

# it depends on:
##  - "code/dependencies.R"
##  - the community data at 
##    my_path, "/community_survey/raw/
##    Jaundice deaths_adults_VA_Bangladesh_selected variables.xls"
##    my_path,"/community_survey/raw/
##    Deaths AJS by cluster and age_catchment area_v2.csv"
##    my_path,"/community_survey/raw/
##    HEV mortality Bangladesh_dataid description.xlsx"

# it produces:
## - the clean community survey data describing number of AJ deaths at 
##   my_path,"/community_survey/clean/clean_community_decedant_level.rds"
##   my_path, "/community_survey/clean/clean_community_cluster_level.rds"

# load dependencies

source("code/dependencies.R")

#-------------------------------------------------------------------------------

#############################################################################
# deaths with acute jaundice recalled by the community - 1 row per decedent #
#############################################################################
# read in the data
## original data from Repon
## filter on just acute jaundice cases
## delete duplicated row
deaths_comm <- readxl::read_xls(
  paste0(
    my_path, 
    "/community_survey/Jaundice deaths_adults_VA_Bangladesh_selected variables.xls")
  )%>% filter(
    acute_chronic == "Acute")%>%
  distinct()

# correct mismatches of clusterid and hospital id (see quarto for process of identification)
# according to Repon
deaths_comm[which(deaths_comm$dataid == "112233"), "clusterid"] <- "24"
deaths_comm[which(deaths_comm$dataid == "224244"), "dataid"] <- "124244"
deaths_comm[which(deaths_comm$dataid == "124244"), "hospitalid"] <- "1"
deaths_comm[which(deaths_comm$dataid == "260205"), "dataid"] <- "360205"
deaths_comm[which(deaths_comm$dataid == "360205"), "hospitalid"] <- "3"
deaths_comm[which(deaths_comm$dataid == "310241"), "dataid"] <- "358241"
deaths_comm[which(deaths_comm$dataid == "358241"), "clusterid"] <- "58"

## age and date of death from Repon
add_comm <- read.csv(
  paste0(
    my_path,"/community_survey/Deaths AJS by cluster and age_catchment area_v2.csv"
    ), colClasses = c(clusterid = "character")
  )%>% mutate(
    dataid = as.character(dataid))

## merge age and date data into the original data
deaths_comm <- full_join(
  deaths_comm %>% filter(
    acute_chronic == "Acute"),
  add_comm, by = c("dataid", "maternal", "clusterid")
  )

# We have a .csv file with 7 columns and 611 rows, where each row represents an
# acute jaundice related death recalled in the community mortality survey. 
# Jaundice was defined as having new onset of either yellow sclera or skin 
# during the illness prior to death.

# The first 6 rows are shown here:
head(deaths_comm)

# Columns should be read as follows:
#   
#   -   dataid\
# -   1st digit = hospital catchment area
# -   2nd and 3rd digit = cluster id
# -   4th digit = type of decedent (1 = maternal, 2 = adult \>14 yrs)
# -   5th and 6th digit = decedent id

#   -   hospitalid = id of hospital catchment area where:
# -   1 = Shaheed Ziaur Rahman MCH, Bogra
# -   2 = Sher-e-Bangla MCH, Barisal
# -   3 = Jahurul Islam MCH, Kishoregonj
# -   4 = Chittagong MCH, Chittagong
# -   5 = Sylhet MCH, Sylhet
# -   6 = Sir Salimullah MCH, Dhaka

table(deaths_comm$hospitalid, useNA = "always")

# could be 58 or 60, have gone for 58
deaths_comm[which(deaths_comm$dataid == "310240"), "dataid"] <- "358240"
deaths_comm[which(deaths_comm$dataid == "358240"), "clusterid"] <- "58"

# -   cluster and cluster ID (from 01-91)
# -   there are only 84 clusters in the data set - 
# -   note cluster 87 is missing as there were no jaundice associated deaths and
# -   clusters 08, 44, 47, 53, 64, 82 and 87
missing <- c("08", "44", "47", "53", "64", "82", "87")
# -   previous data set had both `cluster` and `clusterid` with two discrepancies. 
# -   I had assumed `cluster` was correct as it matched the 2nd and 3rd digits 
# -   of the dataid. The updated dataset only has `clusterid`. 
# -   I still correct the two that mismatch the 2nd and 3rd digits of the dataid

# check length and range
length_c <- length(unique(deaths_comm$clusterid))

# create new cluster variable based on dataid
deaths_comm$cluster <- str_sub(deaths_comm$dataid, start = 2L, end = 3L)

# discrepancies between cluster and cluster id
deaths_comm[deaths_comm$cluster != deaths_comm$clusterid, ] # two incident
# one is 112233 which I can't change to match Repon's corrected cluster 
# as it already is an ID for another decedant
# the second is 201211 which the cluster should be 01:
deaths_comm$clusterid[which(deaths_comm$dataid == "201211")] <- sprintf(
  "%02d", as.numeric(01)
  )

# -   `decedanttype`
# -   1 = maternal death
# -   2 = adult death (\>= 14yrs)

# -   `decedentid` (from 01-99) - this is per cluster and the starting digit was 
# -   assigned to an interviewer so that is the reason for the structure seen 
#     in the histogram below:

ggplot(deaths_comm)+
  geom_histogram(aes(x = as.numeric(decedantid)), binwidth = 1)+
  theme_minimal()+
  xlab("decedentid")

#   -   `maternal`
# -   1 = maternal death
# -   0 = adult death (\>= 14yrs)
# -   there are 5 discrepancies between `maternal` and `decedanttype`. 
# -   Three of them do not have acute jaundice, so just `358203` and `360201` 
# -   needed to be checked with Repon - he confirmed they are both maternal and 
# -   should be classified as 1

# correct maternal classification:
deaths_comm %>% filter(
  (deaths_comm$decedanttype == "2" & deaths_comm$maternal != 0) |
    (deaths_comm$decedanttype == "1" & deaths_comm$maternal != 1))
deaths_comm[which(deaths_comm$dataid == "358203"),"decedanttype"] <- "1"
deaths_comm[which(deaths_comm$dataid == "360201"),"decedanttype"] <- "1"

#   -   `acute_chronic` which has options `Acute`, `No Jaundice` and `Chronic`.

#######################
# cluster information #
#######################

# We also have a look-up table with 8 columns describing each cluster in terms of 
# admin units, rural/urban designation and population size, and 91 rows, 
# where each row is a cluster. Example of the first 6 rows shown below:

# names and location of each cluster (clusters are unions or wards)
clusters <- readxl::read_xlsx(
  paste0(
    my_path,"/community_survey/HEV mortality Bangladesh_dataid description.xlsx"
         ), sheet = "91 Cluster"
  )%>% rename(
    "cluster_name" = `Cluster Name\r\n(Union/Ward)`,
         "clusterid" = `Cluster ID`,
         "pop11" = `Population in 2011 Census`,
         "rural_urban" = `Rural/ Urban`)
head(clusters)

### Merging the look-up table with the list of deaths:

# The 91 clusters are spread over 64 Upazilas (90 clusters spread over 63 if 
# we don't count cluster 87 where there were no jaundice deaths - 
# but we want to keep all).

# merge 
df <- merge(deaths_comm, clusters, by = "clusterid", all.x = T, all.y = T) %>%
  dplyr::select(-cluster) # for clarity

### Sanity checks:

# 1. Any unmatched rows?
sum(is.na(df$MCH)) # check a column from look-up table
sum(is.na(df$hospitalid)) # check a column from the data - 
# its just the clusters with no acute jaundice deaths
# There are no rows in the data that could not be matched up with a cluster ID 
# in the look-up table.

# 2.  Does the hospital catchment area pair one-to-one with the MCH after 
# Repon's corrections for acute jaundice cases?
df %>% filter(
  acute_chronic == "Acute"
  ) %>% group_by(
    hospitalid, MCH
    ) %>% summarise(
      n()
      )
# yes

# Save the clean data
### by decedant
saveRDS(
  df, file = paste0(
    my_path,"/community_survey/clean/clean_community_decedant_level.rds")
  )

# Summarise acute jaundice deaths by cluster
### split into adult and maternal deaths
# tmp <- df %>% filter(acute_chronic == "Acute")%>% 
#   group_by(clusterid, hospitalid, MCH, Division, District, Upazila, 
#   cluster_name, rural_urban, pop11) %>%
#   summarise(maternal = sum(decedanttype == 1),
#             adult = sum(decedanttype == 2))


tmp <- df %>% group_by(
  clusterid, hospitalid, MCH, Division, District, Upazila,
  cluster_name, rural_urban, pop11, pop_upazila_tot,
  pop_upazila_rural, pop_upazila_urban, pop_upazila_urban_other
  ) %>% summarise(
    maternal = sum(decedanttype == 1 & acute_chronic == "Acute"),
    adult = sum(decedanttype == 2 & acute_chronic == "Acute")
    )

# clean clusters with no deaths - replace NA with 0 as reason not included 
# was because no jaundice deaths
tmp[tmp$clusterid %in% missing, c("maternal", "adult")] <- 0
tmp[tmp$clusterid == "08", "hospitalid"] <- as.character(2)
tmp[tmp$clusterid == "44", "hospitalid"] <- as.character(4)
tmp[tmp$clusterid == "47", "hospitalid"] <- as.character(4)
tmp[tmp$clusterid == "53", "hospitalid"] <- as.character(4)
tmp[tmp$clusterid == "64", "hospitalid"] <- as.character(6)
tmp[tmp$clusterid == "82", "hospitalid"] <- as.character(5)
tmp[tmp$clusterid == "87", "hospitalid"] <- as.character(5)


# project population
tmp <- tmp %>%
  mutate(pop12 = 1.0137 * pop11)%>%
  mutate(pop13 = 1.0137 * pop12)%>%
  mutate(pop14 = 1.0137 * pop13)%>%
  mutate(pop14_14orOver = 0.677 * pop14,
         live_births = 0.0222 * pop14,
         nonpreg_pop14_14orOver = round(pop14_14orOver - live_births, 0),
         preg_pop14 = round(live_births, 0))

saveRDS(
  tmp, file = paste0(
    my_path, "/community_survey/clean/clean_community_cluster_level.rds")
  )
