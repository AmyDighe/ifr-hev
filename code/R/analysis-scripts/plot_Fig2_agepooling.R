# This script plots Figure 2 in the manuscript draft

# It depends on:

# It produces:

# load libraries
source(here("code/dependencies.R"))

#-------------------------------------------------------------------------------


# read in pop sizes for plotting stats /10,000 people --------------------------

pop_group <- readRDS(file = 
                       here(
                         paste0(my_path,"/community_survey/pop_group.rds"))
                     )%>% mutate(
                       MCH = recode(MCH, `Chittagong (CMCH)` = "Chittagong",
                      `Kishoreganj (JIMCH)` = "Kishoreganj",
                      `Sylhet (SOMCH)` = "Sylhet",
                      `Barisal (SBMCH)` = "Barisal",
                      `Mitford (SSMCH)` = "Dhaka",
                      `Bogra (SZRMCH)` = "Bogra"),
         group = if_else(
           preg_status == "pregnant", "pregnant", age_class
           ),
         pop = if_else(
           preg_status == "pregnant", round(live_births,0), 
                       round(nonpreg_pop14_14orOver*prop, 0))
         )

# add in pop overall to match IFR panel
pop_overall <- pop_group %>%
  group_by(age_class, preg_status, group)%>%
  summarise(pop = sum(pop))%>%
  mutate(MCH = "Overall")
pop_group <- bind_rows(pop_group, pop_overall)


# read in deaths----------------------------------------------------------------
deaths <- as.data.frame(
  rstan::extract(
    readRDS(
      here(
        "outputs/joint_fit_age_pooling99.rds")
        ), pars = "deaths_sim_phevd"))
groups <- readRDS(
  here(
  "outputs/groups_ifr.rds")
)

deaths$iteration <- 1:dim(deaths)[1]
groups$groupNo <- 1:24
deaths <- deaths %>%
  pivot_longer(cols = 1:24, values_to = "deaths", names_to = "groupNo")%>%
  mutate(groupNo = as.numeric(gsub("[^0-9]", "", groupNo)))
tmp <- left_join(deaths, groups, by = "groupNo")%>%
  mutate(group = if_else(preg_status == "pregnant", "pregnant", age_class))

# add in deaths overall to match overall IFR panel
deaths_overall <- tmp %>%
  group_by(iteration, age_class, preg_status, group)%>%
  summarise(deaths = sum(deaths))%>%
  mutate(MCH = "Overall")
deaths <- bind_rows(tmp, deaths_overall)
plotDeaths <- left_join(deaths, pop_group, by = c("MCH", "group"))
plotDeaths <- plotDeaths %>%
  mutate(death_rate = 10000*deaths/(3*pop),
         MCH = factor(MCH, levels = c("Overall", "Barisal", "Bogra", "Chittagong",
                                      "Dhaka", "Kishoreganj", "Sylhet")))%>%
  group_by(MCH, group, pop)%>%
  summarise(mean_deaths = mean(death_rate),
            low_deaths = quantile(prob = 0.025, death_rate),
            high_deaths = quantile(prob = 0.975, death_rate),
            total_mean_deaths = mean(deaths),
            total_low_deaths = quantile(prob = 0.025, deaths),
            total_upp_deaths = quantile(prob = 0.975, deaths))



# read in infections -----------------------------------------------------------

infections <- readRDS(file = "outputs/infections.rds")%>%
  mutate(group = if_else(preg_status == "pregnant", "pregnant", age_class))

# merge deaths and infections with population size of sampled areas
dat <- merge(
  deaths, infections, by = c("iteration", "group", "MCH",
                             "preg_status", "age_class"),
             all.x = T)
dat <- merge(dat, pop_group %>% 
               dplyr::select(MCH, group, pop, susc_pop), 
             by = c("MCH", "group"), all.x = T) %>%
  mutate(death_rate = 10000*deaths/(3*pop),
         infection_rate = 10000*infections/pop,
         MCH = factor(
           MCH, levels = c("Overall", "Barisal", "Bogra", "Chittagong",
                           "Dhaka", "Kishoreganj", "Sylhet")))%>%
  group_by(MCH, group, pop)%>%
  summarise(mean_deaths = mean(death_rate),
            low_deaths = quantile(
              prob = 0.025, death_rate),
            high_deaths = quantile(
              prob = 0.975, death_rate),
            mean_infections = mean(infection_rate),
            low_infections = quantile(
              prob = 0.025, infection_rate, na.rm = T),
            high_infections = quantile(
              prob = 0.975, infection_rate, na.rm = T))




# read in FoI estimates---------------------------------------------------------
FoI <- matrix(nrow = 6, ncol = 3)
for(i in 1:6){
  FoI[i,] <- readRDS(
    here(
      paste0("outputs/FoI/MCH", i,
             "age_pooling_summary.rds"))
    )$summary["foi1", c("mean", "2.5%", "97.5%")]
}
FoI <- as.data.frame(FoI)
names(FoI) <- c("mean", "lower", "upper")
FoI$MCH <- c("Barisal", "Bogra", "Chittagong", "Kishoreganj", "Dhaka", "Sylhet")
FoI$group <- "60+"
FoI$y <- 1.0

saveRDS(FoI, file = here(
  "generated_data/FoI_summary_age_pooling.rds")
  )

# 2B. HEV deaths /10000 by catchment by group-----------------------------------
p1 <- ggplot(plotDeaths %>% filter(MCH != "Overall"))+
  geom_bar(aes(y = mean_deaths, x = group, fill = group), stat= "identity")+
  coord_flip()+
  geom_errorbar(aes(ymin = low_deaths, ymax = high_deaths, x = group), width = 0, 
                col = "grey70", linetype = 1)+
  facet_wrap(~MCH, nrow = 7)+
  ylab("Annual HEV deaths/10,000 people")+
  xlab(" ")+
  ylim(0,4)+
  scale_fill_viridis_d()+
  theme_minimal()+
  theme(legend.position = "none", strip.text = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.line.x = element_line(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

p1b <- ggplot(plotDeaths %>% filter(MCH == "Overall"))+
  geom_bar(aes(y = mean_deaths, x = group, fill = group), stat= "identity")+
  coord_flip()+
  geom_errorbar(aes(ymin = low_deaths, ymax = high_deaths, x = group), width = 0, 
                col = "grey70", linetype = 1)+
  ylab("Annual HEV deaths/10,000 people")+
  xlab(" ")+
  ylim(0,4)+
  scale_fill_viridis_d()+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.line.x = element_line(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
ptot1 <- plot_grid(p1b,p1, nrow = 2, rel_heights = c(1.2,6), align = "v")


# 2A. FoI by catchment----------------------------------------------------------
p2 <- ggplot(FoI)+
  geom_point(aes(x = mean, y = group), col = "red", size = 2)+
  geom_errorbar(aes(
    xmin = lower, xmax = upper, y = group), 
    width = 0, col = "red", linewidth = 1)+
  facet_wrap(~MCH, nrow = 6, strip.position = "left")+
  xlab("Annual risk of infection")+
  scale_x_continuous(breaks = c(0.05, 0.15, 0.25), 
                     labels = c("0.05", "0.15", "0.25"), 
                     expand = c(0.01, 0.01))+
  scale_y_discrete(limits=rev)+
  theme_minimal()+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_text(size = 11),
        axis.title.x = element_text(size = 10),
        axis.line.x = element_line())

# 2C. IFR overall and by catchment by group-------------------------------------

ifr_fit <- readRDS("outputs/joint_fit_age_pooling99.rds")
ifr_summary <- as.data.frame(summary(ifr_fit, pars = "p")$summary)
ifr_summary$parameter <- rownames(ifr_summary)

plotDat <- bind_cols(ifr_summary, groups) %>%
  mutate(group = if_else(preg_status == "pregnant", "pregnant", age_class))


p3 <- ggplot(data = plotDat, aes(y = group, color = group)) +
  geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0, size = 1) +
  scale_colour_viridis_d()+
  geom_point(aes(x = mean), size = 2) +
  theme_minimal()+
  coord_cartesian(xlim= c(0.00, 0.0031))+
  xlab("Infection fatality risk (IFR)")+
  facet_wrap(~MCH, nrow = 6)+
  theme(axis.text.y = element_blank(),
        strip.text = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10),
        legend.position = "none",
        axis.line.x = element_line(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

# Look at the betas-------------------------------------------------------------

as.data.frame(exp(summary(ifr_fit, pars = "beta")$summary))


# Extract mean IFR for each group across all places-----------------------------

chains <- rstan::extract(ifr_fit)
alpha_overall <- matrix(nrow = 12000, ncol = 3)

for(i in 1:3){
  alpha_overall[,i] <- chains$alpha_mu + chains$beta[,i]*1
}

# Add mean for reference group (14-45 non-pregnant)
p_overall.df <- as.data.frame(plogis(cbind(alpha_overall, chains$alpha_mu)))
names(p_overall.df) <- c("45-59", "60+", "pregnant", "14-44")

# Extract mean and 95% CrI
p_overall <- p_overall.df %>% 
  pivot_longer(cols = 1:4, names_to = "group", values_to = "mean_ifr") %>%
  group_by(group)%>%
  summarise(`2.5%` = quantile(mean_ifr, prob = 0.025),
            mean = mean(mean_ifr),
            `97.5%` = quantile(mean_ifr, prob = 0.975),
            MCH = "Overall")
saveRDS(p_overall, "outputs/p_overall_age_pooling.rds")

# Plot overall------------------------------------------------------------------
p4 <- ggplot(p_overall)+
  geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`, y = group), height = 0, size = 1) +
  geom_point(aes(x = mean, y = group), size = 5, pch = 18) +
  theme_minimal()+
  xlim(0.00, 0.0031)+
  xlab("Infection fatality risk (IFR)")+
  theme(strip.text = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none",
        axis.line.x = element_line(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

ptot2 <- plot_grid(p4,p3, nrow = 2, rel_heights = c(1.2,6), align = "v")

ptot3 <- plot_grid(ptot1, ptot2, nrow = 1, rel_widths = c(10/17,1))

ggsave(ptot3, file = "figures/Fig2BCagepooling.png", height = 7.2, width = 6.4, dpi = 300)

ggsave(p2, file = "figures/Fig2Aagepooling.png", height = 6, width = 2, dpi = 300)
