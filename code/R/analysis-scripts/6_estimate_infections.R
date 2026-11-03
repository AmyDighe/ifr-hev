# This script estimates the HEV infections per catchment area during the mortality
# survey period using the estimated FoI, seroprevalence, and total population size

# It depends on:
## - the community mortality data at the cluster level:
##   at my_path/community_survey/clean_community_cluster_level.rds"
## - age distributions from the 2011 census
## - seroprevalence estimates at "outputs/sp_for_catalytic_model_fitting.rds"
## - FoI estimates (the chains) at outputs/FoI/ (either core estimates, or nosr)


# It produces:
## - estimates numbers of infections per group and catchment as input for IFR
##   model, at outputs/infections.rds and 
##   outputs/infections_nosr.rds for sensitivity analysis with no seroreversion
## - estimates summarised for ease of quoting in the manuscript text
##   at outputs/infections_summary_for_text(/_nosr).rds


# load libraries
source(here::here("code/dependencies.R"))

#-------------------------------------------------------------------------------

# read in pregnant and non-pregnant population of each cluster
co <- readRDS(
  here(
    paste0(my_path, "/community_survey/clean_community_cluster_level.rds"))
)

# add column tying hospital to sp division to pull age distribution
co <- co %>%
  mutate(eff_div = case_when(MCH == "Barisal (SBMCH)" ~ "Barisal.Division",
                             MCH == "Bogra (SZRMCH)" ~ "Rajshahi.Division",
                             MCH == "Chittagong (CMCH)" ~ "Chittagong.Division",
                             MCH == "Sylhet (SOMCH)" ~ "Sylhet.Division",
                             MCH == "Kishoreganj (JIMCH)" ~ "Dhaka.Division",
                             MCH == "Mitford (SSMCH)"~ "Dhaka.Division")) %>%
  group_by(eff_div, MCH)%>%
  summarise(pop14_14orOver = sum(pop14_14orOver),
            live_births = sum(live_births),
            nonpreg_pop14_14orOver = sum(nonpreg_pop14_14orOver),
            preg_pop14 = sum(preg_pop14))

# read in age distribution for each division from: 
# http://203.112.218.65:8008/WebTestApplication/userfiles/Image/PopMonographs/Volume-9_Age-Sex.pdf 
age_dist <- data.frame(
  "age" = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", 
            "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70+"),
  "Barisal.Division" = c(10.36, 13.31, 12.70, 8.38, 7.83, 8.31, 6.92, 6.38, 5.53, 4.49, 4.09, 2.70, 3.26, 1.77, 3.96),
  "Chittagong.Division" = c(11.47, 13.79, 12.79, 10.07, 9.26, 8.66, 6.44, 5.84, 4.96, 3.88, 3.47, 2.18, 2.60, 1.35, 3.24),
  "Dhaka.Division" = c(10.15, 12.19, 11.02, 8.83, 9.96, 10.04, 7.62, 6.76, 5.80, 4.36, 3.84, 2.31, 2.65, 1.37, 3.08),
  "Rajshahi.Division" = c(9.59, 11.56, 10.71, 8.36, 9.24, 9.76, 7.82, 7.28, 6.36, 4.91, 4.10, 2.66, 2.80, 1.58, 3.27),
  "Sylhet.Division" = c(12.82, 14.79, 12.39, 9.66, 8.56, 8.16, 6.33, 5.82, 5.13, 3.76, 3.51, 2.00, 2.60, 1.22, 3.24))
age_dist_match <- age_dist %>% mutate(
  age_class = case_when(
    (age == "15-19" | age == "20-24" | age == "25-29" | age == "30-34" | 
       age == "35-39" | age == "40-44") ~ "14-44",
    (age == "45-49" | age == "50-54" | age == "55-59") ~ "45-59",
    (age == "60-64" | age == "65-69" | age == "70+") ~ "60+",
    TRUE ~ as.character(NA))) %>%
  filter(!is.na(age_class)) # get rid of ages not sampled

# get per over 14 to normalise other % since we only have pop>14
total_perc_over_14 <- colSums(age_dist_match[2:6]) 
age_dist_match[,2:6] <- mapply('/', age_dist_match[,2:6], total_perc_over_14)

age_dist_match <- age_dist_match %>% 
  pivot_longer(cols = 2:6, names_to = "eff_div", values_to = "prop") %>%
  group_by(age_class, eff_div)%>%
  summarise(prop = sum(prop))


# read in seroprevalence estimates
sp <- as.data.frame(readRDS(
  here(
  "outputs/sp_for_catalytic_model_fitting.rds")
  )%>% ungroup())
sp$mean <- rowMeans(sp[,c(3:1002)])
mch <- unique(sp$MCH)

# summarise seroprevalence estimates by age and MCH
sp <- sp %>% dplyr::select(age_class, age1, age2, mean, MCH) %>%
  mutate(
    age_class = case_when(
      (age_class == "(15,20]" | age_class == "(20,25]" | age_class == "(25,30]" | 
         age_class == "(30,35]" | age_class == "(35,40]" | age_class == "(40,45]") ~ "14-44",
      age_class == "(45,110]" ~ "45-59",
      TRUE ~ as.character(NA)))%>%
  filter(!is.na(age_class)) %>%
  group_by(age_class, MCH) %>%
  summarise(sp = mean(mean))

# stick on an end for 60+ using same seroprev as 45-110
sp60 <- sp %>% filter(age_class == "45-59")%>%
  mutate(age_class = "60+")
sp <- bind_rows(sp, sp60) %>%
  mutate(preg_status = "not_pregnant")

spPreg <- sp %>% filter(age_class == "14-44")%>%
  mutate(preg_status = "pregnant")

sp <- bind_rows(sp, spPreg) %>%
  mutate(eff_div = case_when(MCH == "Barisal (SBMCH)" ~ "Barisal.Division",
                             MCH == "Bogra (SZRMCH)" ~ "Rajshahi.Division",
                             MCH == "Chittagong (CMCH)" ~ "Chittagong.Division",
                             MCH == "Sylhet (SOMCH)" ~ "Sylhet.Division",
                             MCH == "Kishoreganj (JIMCH)" ~ "Dhaka.Division",
                             MCH == "Mitford (SSMCH)"~ "Dhaka.Division"))

## merge age dist, population and sp --> 
## seronegative population in each age bin/pregnancy status and division

DAT <- merge(sp, age_dist_match, by = c("eff_div", "age_class"), all.x = T)

DAT2 <- merge(DAT, co, by = c("MCH"), all.x = T)

DAT3 <- DAT2 %>%
  mutate(susc_pop = if_else(
    preg_status == "not_pregnant", 
    nonpreg_pop14_14orOver * (1-sp) * prop, 
    preg_pop14 * (1-sp)))

# save for pop
saveRDS(DAT3, file = here(
  paste0(my_path, "/community_survey/pop_group.rds"))
)

# read in posterior estimates for the FoI of each division
FoI_samp <- matrix(ncol = 6, nrow = 1000)
colnames(FoI_samp) <- mch
div_lookup <- data.frame(
  Division = colnames(FoI_samp), idx = 1:6)

for(i in 1:6){
  FoI <- rstan::extract(readRDS(
    here(
      paste0("outputs/FoI/MCH", i, "age_pooling_chains_nosr.rds")))
    )$foi1
  FoI_samp[,i] <- sample(FoI, 1000)
}

FoI_samp <- as.data.frame(FoI_samp)
names(FoI_samp) <- mch
FoI_samp$iteration <- 1:dim(FoI_samp)[1]
FoI_samp <- FoI_samp %>% pivot_longer(
  cols = 1:6, values_to = "FoI", names_to = "MCH")
FoI_samp <- FoI_samp %>% pivot_wider(
  values_from = "FoI", names_from = "iteration")

# long_df <- FoI_samp %>%
#   pivot_longer(
#     cols = -MCH,
#     names_to = "iteration",
#     values_to = "foi"
#   ) %>%
#   mutate(iteration = as.integer(iteration))
# 
# ggplot(long_df)+
#   geom_histogram(aes(x = foi))+
#   facet_wrap(~MCH)

## merge FoI and add mean infections over 3 years
DAT4 <- merge(FoI_samp, DAT3, by = "MCH", all.x = T)

## make long by iteration
DAT4 <- DAT4 %>%
  pivot_longer(cols = 2:1001, names_to = "iteration", values_to = "foi")

DAT4 <- DAT4 %>%
  mutate(infections = as.integer(round(susc_pop*foi*3, 0)),
         MCH = recode(MCH, `Barisal (SBMCH)` = "Barisal",
                      `Bogra (SZRMCH)` = "Bogra",
                      `Chittagong (CMCH)` = "Chittagong",
                      `Kishoreganj (JIMCH)` = "Kishoreganj",
                      `Mitford (SSMCH)` = "Dhaka",
                      `Sylhet (SOMCH)` = "Sylhet"),
         iteration = as.integer(iteration))%>%
  dplyr::select(iteration, MCH, infections, age_class, preg_status)%>%
  arrange(iteration, preg_status, age_class, MCH)


# ## plot to check
# DAT4 <- DAT4%>%
#   mutate(group = paste0(age_class, preg_status))
# p4 <- ggplot(DAT4, aes(y = group, x = infections, fill = group))+
#   geom_density_ridges(stat="binline", binwidth = 1000, scale = 0.75,
#                       draw_baseline = F, alpha = 0.7, color = NA)+
#   scale_fill_viridis_d(name = " ")+
#   theme_minimal() +
#   theme(text = element_text(size = 20),
#         legend.position = "none")+
#   facet_wrap(~MCH, nrow = 6)+
#   xlab("HEV infections")+
#   ylab("")+
#   scale_y_discrete(
#         labels = c("14-44 yrs", "pregnant", "45-59 yrs", "60+ yrs"))

## save infections
saveRDS(DAT4, file = here(
  "outputs/infections_nosr.rds")
)
## summarise infections for text
foo <- readRDS("outputs/infections_agepooling_nosr.rds")%>%
  group_by(iteration, MCH)%>%
  summarise(infections = sum(infections))%>% 
  group_by(MCH) %>% summarise(mean = mean(infections), 
                              q2.5 = quantile(infections, 0.025), 
                              q97.5 = quantile(infections, 0.975))
saveRDS(foo, 
        file = here(
          "outputs/infections_summary_for_text_nosr.rds")
)
