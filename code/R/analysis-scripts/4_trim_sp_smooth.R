# This script takes the smoothed age stratified seroprevalence data and cuts out the 
# catchment areas pertaining to the mortality data

# It depends on:
# - posterior draws from 3_fit_sp_smooth: 
#   at "outputs/seroprevalence/sp_post_agepooling.rds"
# - shapefiles of the catchment areas at "outputs/catchment.rds"

# It produces:
# - Figure for S2B "figures/FigS2B.png"
# - Figure for S2C "figures/FigS2C.png"
# - seroprevalence for plotting at "outputs/seroprevalence/sp_for_plotting.rds"
# - seroprevalence reformatted for FoI estimation
#   at "outputs/sp_for_catalytic_model_fitting.rds")

source("code/dependencies.R")

#-------------------------------------------------------------------------------

# read in smoothed seroprevalence data - 
# posterior draws for each grid cell for each age class
pdraws <- readRDS(
  "outputs/seroprevalence/sp_post_agepooling.rds")

pdraws <- pdraws%>%
  mutate(age_class = case_when(age_id == 1 ~ "(0,5]",
                               age_id == 2 ~ "(5,10]",
                               age_id == 3 ~ "(10,15]",
                               age_id == 4 ~ "(15,20]",
                               age_id == 5 ~ "(20,25]",
                               age_id == 6 ~ "(25,30]",
                               age_id == 7 ~ "(30,35]",
                               age_id == 8 ~ "(35,40]",
                               age_id == 9 ~ "(40,45]",
                               age_id == 10 ~ "(45,110]"))%>%
  mutate(age_class = factor(age_class, levels = c("(0,5]","(5,10]","(10,15]",
                                                  "(15,20]",
                                                  "(20,25]",
                                                  "(25,30]",
                                                  "(30,35]",
                                                  "(35,40]",
                                                  "(40,45]",
                                                  "(45,110]")))

# cut out age stratified sp per catchment

### read in the catchment areas
bd3_merge <- readRDS(file = "outputs/seroprevalence/catchment.rds")%>%
  transform_to_btm()

### merge the upazilas in each catchment areas together
catchment <- bd3_merge %>%
  group_by(MCH) %>% 
  summarise()

### trim grids to the catchment areas
trm <- st_intersection(pdraws, catchment)

### quick sense check
catchment <- ggplot()+
  geom_sf(data = bgd0, colour = NA)+
  geom_sf(data = trm %>% filter(age_class == "(40,45]"), aes(fill = pred_1), colour = NA)+
  scale_fill_viridis_c(direction=1,alpha=.9, limits = c(0,1), name = "Seroprevalence")+
  theme_void()
ggsave(catchment, file = "figures/FigS2B.png", dpi =300)

# calculate mean for each catchment area per iteration and age class
t <- trm %>% 
  st_drop_geometry()%>% 
  group_by(MCH, age_class)%>%
  summarise(across(names(trm)[9:1008], mean))

# plot distribution by age for each catchment
t_long <- pivot_longer(t, cols = 3:1002, names_to = "iter", values_to = "sp")%>%
  mutate(age_class = recode(age_class, `(0,5]` = "0-5",
                            `(5,10]` = "5-10",
                            `(10,15]` = "10-15",
                            `(15,20]` = "15-20",
                            `(20,25]` = "20-25",
                            `(25,30]` = "25-30",
                            `(30,35]` = "30-35",
                            `(35,40]` = "35-40",
                            `(40,45]` = "40-45",
                            `(45,110]` = "45-110"))%>%
  mutate(age_class = factor(age_class, levels = c("0-5", "5-10",
                                                  "10-15", "15-20",
                                                  "20-25", "25-30",
                                                  "30-35", "35-40",
                                                  "40-45", "45-110")),
         pos = round(1000000 * sp, 0))

mch.labs <- c("Barisal", "Bogra", "Chittagong", "Kishoreganj", "Dhaka", "Sylhet")
names(mch.labs)<- c("Barisal (SBMCH)", "Bogra (SZRMCH)", "Chittagong (CMCH)",
                    "Kishoreganj (JIMCH)", "Mitford (SSMCH)","Sylhet (SOMCH)")
p <- ggplot(t_long)+
  geom_point(aes(x = age_class, y = sp), position = position_jitter(width = 0.25), alpha = 0.1)+
  facet_wrap(~MCH, nrow = 2, labeller = labeller(MCH = mch.labs))+
  theme_minimal()+
  ylab("Seroprevalence")+
  xlab("Age class")+
  theme(text = element_text(size = 24),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(p, file = "figures/FigS2C.png", height = 6, width = 9, unit = "in")

# reshape for stan
## add age class boundaries
t$age1 <- as.numeric(sapply( str_match_all(t$age_class, pattern = "\\d+"),"[[",1)) + 1
t$age2 <- as.numeric(sapply( str_match_all(t$age_class, pattern = "\\d+"),"[[",2)) + 0.99

# save
saveRDS(t_long, file = "outputs/seroprevalence/sp_for_plotting.rds")
saveRDS(t, file = "outputs/seroprevalence/sp_for_catalytic_model_fitting.rds")

