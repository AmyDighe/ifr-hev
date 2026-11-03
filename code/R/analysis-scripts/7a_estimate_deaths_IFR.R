# This script fits the joint model estimating the probability of HEV infection
# given death based on the hospital data [P(HEV|DAJ)]], the resultant HEV deaths
# among the acute jaundice deaths recalled by the community, and the IFR.

# It depends on:
## - mortality data at my_path/model_input/deaths.rds
## - estimated infections at outputs/infections.rds
## - hierarchical model at "code/R/rstan-models/joint_mod_deaths_ifr.stan"

# It produces:
## model fits at "outputs/joint_fit_age_pooling.rds"


# load libraries
source(here::here("code/dependencies.R"))

#-------------------------------------------------------------------------------

# Load data --------------------------------------------------------------------

deaths <- readRDS(file = here(
  paste0(
    my_path, "/model_input/deaths.rds")))%>% 
  rename(
    "preg_status" = maternal)%>%
  mutate(
    age_class = recode(
      age_class, `(13,44]` = "14-44",
      `(44,59]` = "45-59",
      `(59,101]` = "60+"),
    preg_status = recode(
      preg_status, NO = "not_pregnant",
      YES = "pregnant"),
    group = if_else(
      preg_status == "not_pregnant", age_class, "pregnant")
    )

infections <- readRDS(here(
  file = here(
    "outputs/infections.rds")
))%>% mutate(
    group = if_else(
      preg_status == "not_pregnant", age_class, "pregnant")
    )
infections_summary <- infections %>%
  group_by(
    MCH, age_class, preg_status, group)%>%
  summarise(
    mean_infections = mean(infections),
    median_infections = median(infections),
    sd_infections = sd(infections))%>%
  mutate(
    mean_infections = if_else(
    group == "pregnant", mean_infections, mean_infections)
    )

# Reshape ----------------------------------------------------------------------
## Merge (safety-net ensuring order is consistent)
data <- full_join(infections_summary, deaths,
                  by = c("MCH", "group", "age_class", "preg_status"))

## Unique groups (n MCH x n age classes + n places)
groups <- data %>% 
  dplyr::select(-DAJ, -DHEV, -DcAJ, -mean_infections, - median_infections, -sd_infections)
places <- unique(groups$MCH) # for indexing MCH correctly

# Make model matrix -------------------------------------------------------
X_full <- model.matrix(~ MCH + age_class + preg_status - 1, data = groups)

## Remove MCH as this is modeled as random effect
X <- X_full[, -which(str_detect(colnames(X_full), "MCH"))]


# Fit ---------------------------------------------------------------------

## Specify imperfect test parameters
sens <- 0.99 # Update with Aybuke's work ~0.95 - divergent transitions issue
spec <- 0.99 # Update with Aybuke's work ~0.95 - divergent transitions issue

ifr_fit <- stan(
  file = here::here("code/R/rstan-models/joint_mod_deaths_ifr.stan"),
  data = list(
    J = nrow(groups),
    L = length(unique(data$MCH)),
    M = ncol(X),
    DHEV = data$DHEV,
    DAJ = data$DAJ,
    DcAJ = data$DcAJ,
    I_post_mean = data$mean_infections,
    I_post_sd = data$sd_infections,
    X = X,
    sens = sens,
    spec = spec,
    map_to_place = map_dbl(groups$MCH, ~ (which(places == .)))),
  chains = 4,
  warmup = 1000,
  iter = 4000,
  save_warmup = FALSE,
  verbose = TRUE,
  cores = 4,
  control = list(adapt_delta = 0.98, 
                 max_treedepth = 15))

# Save --------------------------------------------------------------------
saveRDS(ifr_fit, file = here(
  "outputs/joint_fit_age_pooling99.rds")
)

