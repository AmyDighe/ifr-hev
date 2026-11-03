# This script takes age stratified seroprevalence posterior estimates and fits
# catalytic models of seroconversion to estimate the annual risk of infection
# among seronegative individuals, propagating uncertainty

# It depends on:
# - the age-stratified seroprevalence estimates for each catchment area 
#   at "outputs/seroprevalence/sp_for_catalytic_model_fitting.rds"
# - seroreversion rate estimates from Sitakunda 
#   at "public_data/serochit/sr_posterior_c10.rds"
#   and "public_data/serochit/sr_posterior_a10.rds"
# - serocatalytic model in rstan at "code/R/rstan_models/FoImod.stan"


# It produces:
# - chains and fit summaries in outputs/FoI/
#   note that chains plotted in FigS3
# - figures for FigS4 and FigS5 at 
#   "figures/FigS4.png" and
#   "figures/FigS5.png"

# load libraries
source(here::here(
  "code/dependencies.R")
)

#-------------------------------------------------------------------------------

# read in seroprevalence estimates
dat <- as.data.frame(
  readRDS(
    here::here(
    "outputs/seroprevalence/sp_for_catalytic_model_fitting.rds"))%>% 
    ungroup()
  )
# read in posterior estimates of the seroreversion rate
# from https://pubmed.ncbi.nlm.nih.gov/38497497/
sr_u10 <- rstan::extract(
  readRDS(
    here::here(
    "data_public/serochit/sr_posterior_c10.rds"))
    )$sigma
sr_o10 <- rstan::extract(
  readRDS(
    here::here(
      "data_public/serochit/sr_posterior_a10.rds"))
  )$sigma

## take mean of the posterior
### for low sr sensitivity analyses replace these with 0.05
### for no sr sensitivity analyses replace these with a very small number 
### e.g. 0.00000001
samp_u10 <- mean(sr_u10)
samp_o10 <- mean(sr_o10)

########################
# Fit catalytic model: #
########################

# Assumptions:
### constant FoI
### seroreversion based on empirical measures in Sitakunda in <10 and >10 yr-olds
### 1000 samples from the posterior estimate of seroprevalence for each 
### age class are used

mch <- unique(dat$MCH)

##dat <- dat[,c(1:12, 1003:1004)]

for(i in 1:(length(mch))){
  fit <- stan(
    file = here::here("code", "R", "rstan-models", "FoImod.stan"),
    data = list(
      A = length(unique(dat%>% filter(MCH == mch[i])%>% pull(age_class))),
      Z = 1000,
      sp = dat %>% 
        filter(MCH == mch[i]) %>% dplyr::select(-MCH, -age_class, -age1, -age2),
      age_foi = 200, # set age of change above oldest age class --> constant FoI
      age_sigma = 11, # age at which rate of seroreversion changes
      age1 = dat %>% filter(MCH == mch[i]) %>% pull(age1),
      age2 = dat %>% filter(MCH == mch[i]) %>% pull(age2),
      sigma1 = 0.05, # use mean of posterior sr_u10/o10 if using
      sigma2 = 0.05
    ),
    chains = 4,
    iter = 1000,
    warmup = 250,
    verbose = TRUE,
    cores = 4,
    control = list(adapt_delta = 0.99)
  )
  # save the model output
  sumfit <- rstan::summary(fit)
  saveRDS(
    sumfit, file = here(
      "outputs", "FoI", paste0("MCH", i, "age_pooling_summary_foi_lowsr.rds")
      ))
  saveRDS(
    fit, file = here(
      "outputs", "FoI", paste0("MCH", i, "age_pooling_chains_foi_lowsr.rds")
      ))
}

# diagnostics
p <- list(length = length(mch))

t_long <- readRDS(
  here(
  "outputs/seroprevalence/sp_for_plotting.rds")
)

for(i in 1:(length(mch))){
  ft <- as.data.frame(readRDS(
    here(
      "outputs", "FoI", paste0("MCH", i, "age_pooling_summary_lowsr.rds")
      ))$summary)
  ft <- ft %>% filter(
    grepl(rownames(ft), pattern = "sero")
    )%>% mutate(
      age_class = unique(dat$age_class)
      ) %>% mutate(
        age_class = recode(age_class, `(0,5]` = "0-5",
                              `(5,10]` = "5-10",
                              `(10,15]` = "10-15",
                              `(15,20]` = "15-20",
                              `(20,25]` = "20-25",
                              `(25,30]` = "25-30",
                              `(30,35]` = "30-35",
                              `(35,40]` = "35-40",
                              `(40,45]` = "40-45",
                              `(45,110]` = "45-110")
        )%>% mutate(
          age_class = factor(age_class, 
                             levels = c("0-5", "5-10",
                                        "10-15", "15-20",
                                        "20-25", "25-30",
                                        "30-35", "35-40",
                                        "40-45", "45-110")))
  
  
  p[[i]] <- ggplot()+
    geom_point(data = t_long 
               %>% filter(
                 MCH == mch[i]), aes(x = age_class, y = sp),
               position = position_jitter(width = 0.25), alpha = 0.1)+
    geom_point(data = ft, aes(
      x = age_class, y = mean), col = "red", pch = 17, fill = "red", size = 3)+
    geom_linerange(data = ft, aes(
      x = age_class, ymin = `2.5%`, ymax = `97.5%`), col = "red")+
    theme_classic()+
    xlab("age class")+
    ylab("seroprevalence")+
    ggtitle(mch[i])+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ylim(0,1)
  
}

ptot <- cowplot::plot_grid(
  p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], p[[6]], nrow = 2
  )

ggsave(here(
  "figures/FigS5.png"), 
  ptot, width = 12, height = 10, unit = "in")
