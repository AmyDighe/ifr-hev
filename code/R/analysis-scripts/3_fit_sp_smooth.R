# This script fits a geospatial model of age-stratified seroprevalence to 
# observations from villages sampled in the 2015 national serosurvey, 
# with partial pooling across age groups using INLA and SPDE

# It depends on:
## - code/dependencies.R
## - grid from serosurvey paper at
##   my_path, "/serosurvey/grid_with_covs.rds"
## - data from serosurvey paper at
##   my_path, "/serosurvey/dat_with_raster_covs.rds"


# It produces:
## - posterior draws "outputs/seroprevalence/sp_post_agepooling.rds"
## - cross validation metrics "outputs/seroprevalence/age_pool_cv_all_covs.rds" 
##   and "outputs/seroprevalence/age_pool_cv.rds"
## - the best fit "outputs/seroprevalence/sp_fit_agepooling.rds"
## - a sanity check figure "figures/FigS2A.tif"

source(here::here(
  "code/dependencies.R")
)

## -----------------------------------------------------------
## 1. LOAD SPATIAL GRID
## -----------------------------------------------------------

my_cov_grid <- readRDS(
  here(
  paste0(my_path, "/serosurvey/grid_with_covs.rds"))
)

stopifnot(inherits(my_cov_grid, "sf"))

## centroids used everywhere for INLA locations
grid_centroids <- st_centroid(my_cov_grid)
grid_coords <- st_coordinates(grid_centroids)

## -----------------------------------------------------------
## 2. LOAD RAW DATA
## -----------------------------------------------------------

dat_raw <- readRDS(
  here(
  paste0(my_path, "/serosurvey/dat_with_raster_covs.rds"))
)

## -----------------------------------------------------------
## 3. ASSIGN EACH VILLAGE TO DOMINANT GRID CELL (16 are split)
## -----------------------------------------------------------

village_grid <- dat_raw %>%
  count(villid, grid_id, name = "n") %>%
  group_by(villid) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(villid, grid_id_new = grid_id)

dat_fixed_raw <- dat_raw %>%
  left_join(village_grid, by = "villid") %>%
  mutate(grid_id = grid_id_new) %>%
  dplyr::select(-grid_id_new)

## sanity check: villages no longer split
check <- dat_fixed_raw %>%
  group_by(villid, CommLong, CommLat) %>%
  summarise(n_grid = n_distinct(grid_id), .groups = "drop") %>%
  count(n_grid)

print(check)   # MUST be only n_grid = 1

## -----------------------------------------------------------
## 4. GROUP BY VILLAGE/GRID CELL 
## -----------------------------------------------------------

dat <- dat_fixed_raw %>%
  group_by(grid_id, villid, age_cat_1) %>%
  summarise(
    n_pos = sum(pos),
    n_tested = n(),
    .groups = "drop"
  )

## -----------------------------------------------------------
## 5. ATTACH GRID COVARIATES (since we've updated grid cells)
## -----------------------------------------------------------

dat <- dat %>%
  left_join(
    st_drop_geometry(my_cov_grid),
    by = "grid_id"
  )


##------------------------------------------------------------
## 6. INDEX AGES FOR INLA (it can't use characters)
##------------------------------------------------------------

age_levels <- sort(unique(dat$age_cat_1))

dat$age_cat_1 <- factor(dat$age_cat_1, levels = age_levels, ordered = TRUE)
dat$age_id <- as.integer(dat$age_cat_1)
n_age <- length(age_levels)

## -----------------------------------------------------------
## 7. OBSERVATION COORDINATES (used to inform mesh)
## -----------------------------------------------------------

## attach coordinates to each observation row via grid_id
dat$X <- grid_coords[match(dat$grid_id, my_cov_grid$grid_id), 1]
dat$Y <- grid_coords[match(dat$grid_id, my_cov_grid$grid_id), 2]

stopifnot(!any(is.na(dat$X)))  # sanity check
stopifnot(!any(is.na(dat$Y)))

## -----------------------------------------------------------
## 8. BUILD MESH
## -----------------------------------------------------------
locs <- as.matrix(dat[, c("X", "Y")])

bnd <- inla.nonconvex.hull(
  locs,
  convex = -0.1
)

mesh <- inla.mesh.2d(
  boundary = bnd,
  loc = locs,
  cutoff = 500,
  max.edge = c(20000, 60000),
  offset=c(-0.05, -0.05)
)

plot(mesh)
points(locs, col = "red", pch = 16, cex = 0.5)

## -----------------------------------------------------------
## 9. DEFINE SPDE MODEL
## -----------------------------------------------------------

spde <- inla.spde2.matern(
  mesh = mesh,
  alpha = 2
)

##------------------------------------------------------------
## BUILD PREDICTION MATRIX
##------------------------------------------------------------

pred_grid <- expand.grid(
  grid_id = my_cov_grid$grid_id,
  age_id = 1:n_age
)

## attach coordinates
pred_grid$X <- grid_coords[
  match(pred_grid$grid_id, my_cov_grid$grid_id), 1
]

pred_grid$Y <- grid_coords[
  match(pred_grid$grid_id, my_cov_grid$grid_id), 2
]

## -----------------------------------------------------------
## 10. CREATE SPATIAL INDEX
## -----------------------------------------------------------

spatial_index <- inla.spde.make.index(
  name = "spatial",
  n.spde = spde$n.spde
)

## -----------------------------------------------------------
## 11. PROJECTOR MATRICES
## -----------------------------------------------------------

A_est <- inla.spde.make.A(
  mesh = mesh,
  loc = as.matrix(dat[, c("X","Y")])
)

A_pred <- inla.spde.make.A(
  mesh = mesh,
  loc = as.matrix(pred_grid[, c("X", "Y")])
)

## -----------------------------------------------------------
## BUILD STACKS
## -----------------------------------------------------------

## estimation stack
stack_est <- inla.stack(
  data = list(
    y = dat$n_pos,
    n = dat$n_tested
  ),
  A = list(A_est, 1),
  effects = list(
    list(
      spatial = spatial_index$spatial
    ),
    data.frame(
      Intercept = 1,
      age_id = dat$age_id
    )
  ),
  tag = "est"
)

## prediction stack
stack_pred <- inla.stack(
  data = list(
    y = NA,
    n = 1
  ),
  A = list(A_pred, 1),
  effects = list(
    list(
      spatial = spatial_index$spatial
    ),
    data.frame(
      Intercept = 1,
      age_id = pred_grid$age_id
    )
  ),
  tag = "pred"
)

## combine stacks

stk <- inla.stack(stack_est, stack_pred)

id.prd <- inla.stack.index(stk, "pred")$data

## -----------------------------------------------------------
## DEFINE MODEL
## -----------------------------------------------------------

formula <- y ~ -1 + Intercept +
  
  ## smooth age intercepts
  f(age_id, model = "rw1") +
  
  ## shared Bangladesh-wide spatial field
  f(spatial,
    model = spde
  )

## -----------------------------------------------------------
## FIT MODEL
## -----------------------------------------------------------

fit <- inla(
  
  formula,
  
  family = "binomial",
  
  data = inla.stack.data(stk),
  
  Ntrials = n,
  
  control.predictor = list(
    A = inla.stack.A(stk),
    compute = TRUE,
    link = 1
  ),
  
  control.compute = list(
    dic = TRUE,
    waic = TRUE,
    cpo = TRUE,
    config = TRUE
  ),
  
  verbose = TRUE
)

## -----------------------------------------------------------
## POSTERIOR DRAWS
## -----------------------------------------------------------

sim <- inla.posterior.sample(
  n = 1000,
  result = fit
)

post_draws <- sapply(sim, function(x) {
  plogis(x$latent[id.prd])
})

##------------------------------------------------------------
## FORMAT OUTPUT
##------------------------------------------------------------

post_draws <- as.data.frame(post_draws)

# rename for niceness
colnames(post_draws) <- paste0("pred_", 1:ncol(post_draws))

# attach prediction metadata
post_draws$grid_id <- pred_grid$grid_id
post_draws$age_id  <- pred_grid$age_id

## attach coordinates
post_draws <- merge(my_cov_grid, pdraws, by = "grid_id")%>%
  st_as_sf()

saveRDS(post_draws, here(
  "outputs/sp_post_agepooling.rds")
)
##------------------------------------------------------------
## MAPPING AGE SPECIFIC OUTPUTS
##------------------------------------------------------------

spatial_means <- fit$summary.fitted.values$mean[id.prd]
spatial_sds   <- fit$summary.fitted.values$sd[id.prd]

my_grid_sf <- pred_grid %>%
  mutate(
    mean = spatial_means,
    sd   = spatial_sds
  ) %>%
  left_join(my_cov_grid, by = "grid_id") %>%
  st_as_sf()

my_grid_sf$age_class <- factor(
  my_grid_sf$age_id,
  levels = 1:n_age,
  labels = age_levels
)

## load Bangladesh shapefile for mapping
bgd0 <- load_bgd_shp(level=0) %>%
  transform_to_btm() 

## cut out 
map_grid <-  st_intersection(my_grid_sf, bgd0)

##-------------------------------------------------------------
## PLOT AS A SANITY CHECK
##-------------------------------------------------------------

# base map layer (Bangladesh outline)
base_map <- bgd0

# plot
p <- ggplot() +
  
  geom_sf(data = map_grid, aes(fill = mean), color = NA) +
  
  geom_sf(data = base_map, fill = NA, color = "black", linewidth = 0.3) +
  
  scale_fill_viridis_c(
    option = "viridis",
    limits = c(0, 1),
    name = "Seroprevalence"
  ) +
  
  facet_wrap(~ age_class, ncol = 5) +
  
  coord_sf(datum = NA) +
  
  theme_void() +
  
  theme(
    strip.text = element_text(size = 10),
    legend.position = "bottom"
  )

##------------------------------------------------------------
## ASSESS MODEL FIT USING GRID CELL LOO-CV
##------------------------------------------------------------

# drop age classes with n_pos = NA (n=9)
# dat <- dat %>%
#   filter(
#     !is.na(n_pos)
#   )

############################################################
## GRID-LEVEL CROSS VALIDATION
############################################################

obs_grids <- sort(unique(dat$grid_id))

cv_results <- vector("list", length(obs_grids))

for(k in seq_along(obs_grids)) {
  
  cat("Fold", k, "of", length(obs_grids), "\n")
  
  holdout_grid <- obs_grids[k]
  
  ##--------------------------------------------------------
  ## MASK OUT ONE GRID CELL
  ##--------------------------------------------------------
  
  dat_cv <- dat
  
  dat_cv$n_pos[dat_cv$grid_id == holdout_grid] <- NA
  
  ##--------------------------------------------------------
  ## REBUILD STACK (same structure as main model)
  ##--------------------------------------------------------
  
  stack_cv <- inla.stack(
    
    data = list(
      y = dat_cv$n_pos,
      n = dat_cv$n_tested
    ),
    
    A = list(A_est, 1),
    
    effects = list(
      
      list(
        spatial = spatial_index$spatial
      ),
      
      data.frame(
        Intercept = 1,
        age_id = dat_cv$age_id
      )
      
    ),
    
    tag = "est"
    
  )
  
  ##--------------------------------------------------------
  ## FIT MODEL (same as main model)
  ##--------------------------------------------------------
  
  fit_cv <- inla(
    
    formula,
    
    family = "binomial",
    
    data = inla.stack.data(stack_cv),
    
    Ntrials = n,
    
    control.predictor = list(
      A = inla.stack.A(stack_cv),
      compute = TRUE,
      link = 1
    ),
    
    control.compute = list(
      config = TRUE
    )
    
  )
  
  ##--------------------------------------------------------
  ## EXTRACT HELD-OUT ROWS
  ##--------------------------------------------------------
  
  holdout_rows <- which(dat$grid_id == holdout_grid)
  
  pred <- fit_cv$summary.fitted.values$mean[holdout_rows]
  
  truth <- dat$n_pos[holdout_rows] /
    dat$n_tested[holdout_rows]
  
  cv_results[[k]] <- data.frame(
    grid_id = holdout_grid,
    age_id = dat$age_id[holdout_rows],
    pred = pred,
    truth = truth
  )
  
}

##------------------------------------------------------------
## COMBINE RESULTS
##------------------------------------------------------------

cv_results <- bind_rows(cv_results)

##------------------------------------------------------------
## METRICS
##------------------------------------------------------------

mae <- mean(
  abs(cv_results$pred - cv_results$truth),
  na.rm = TRUE
)

bias <- mean(
  cv_results$pred - cv_results$truth,
  na.rm = TRUE
)

pearson_cor <- cor(
  cv_results$pred,
  cv_results$truth,
  use = "complete.obs"
)

fit_metrics <- c("mae"= mae, "bias"= bias,"pearson_cor" = pearson_cor)

saveRDS(fit_metrics, file = here(
  "outputs/age_pool_cv.rds")
)

## -----------------------------------------------------------
## SAVE OUTPUTS
## -----------------------------------------------------------

saveRDS(
  fit,
  file = here(
    "outputs/sp_fit_agepooling.rds")
)

saveRDS(
  post_draws, 
  file = here(
    "outputs/sp_post_agepooling.rds")
  )

ggsave(
  p,
  filename = here(
    "figures/FigS2A.tif"),
  dpi = 300,
  width = 9, height = 6, unit = "in"
)

############################################################
## GRID-LEVEL CROSS VALIDATION
############################################################

obs_grids <- sort(unique(dat$grid_id))

cv_results_all_covs <- vector("list", length(obs_grids))

formula_all_covs <- y ~ -1 + Intercept +
  
  ## full spatial covariates
  logpop + poverty + travel_time + distance_to_water +
  
  ## smooth age intercepts
  f(age_id, model = "rw1") +
  
  ## shared Bangladesh-wide spatial field
  f(spatial,
    model = spde
  )

for(k in seq_along(obs_grids)) {
  
  cat("Fold", k, "of", length(obs_grids), "\n")
  
  holdout_grid <- obs_grids[k]
  
  ##--------------------------------------------------------
  ## MASK OUT ONE GRID CELL
  ##--------------------------------------------------------
  
  dat_cv <- dat
  
  dat_cv$n_pos[dat_cv$grid_id == holdout_grid] <- NA
  
  ##--------------------------------------------------------
  ## REBUILD STACK (same structure as main model)
  ##--------------------------------------------------------
  
  stack_cv <- inla.stack(
    
    data = list(
      y = dat_cv$n_pos,
      n = dat_cv$n_tested
    ),
    
    A = list(A_est, 1),
    
    effects = list(
      
      list(
        spatial = spatial_index$spatial
      ),
      
      data.frame(
        Intercept = 1,
        age_id = dat_cv$age_id,
        logpop = dat_cv$logpop,
        poverty = dat_cv$poverty,
        travel_time = dat_cv$travel_time,
        distance_to_water = dat_cv$distance_to_water
      )
      
    ),
    
    tag = "est"
    
  )
  
  ##--------------------------------------------------------
  ## FIT MODEL (same as main model)
  ##--------------------------------------------------------
  
  fit_cv <- inla(
    
    formula_all_covs,
    
    family = "binomial",
    
    data = inla.stack.data(stack_cv),
    
    Ntrials = n,
    
    control.predictor = list(
      A = inla.stack.A(stack_cv),
      compute = TRUE,
      link = 1
    ),
    
    control.compute = list(
      config = TRUE
    )
    
  )
  
  id_est <- inla.stack.index(stack_cv, "est")$data
  
  ##--------------------------------------------------------
  ## EXTRACT HELD-OUT ROWS
  ##--------------------------------------------------------
  
  pred_all <- fit_cv$summary.fitted.values$mean[id_est]
  
  holdout_rows <- which(dat$grid_id == holdout_grid)
  
  pred <- pred_all[holdout_rows]
  
  truth <- dat$n_pos[holdout_rows] /
    dat$n_tested[holdout_rows]
  
  cv_results_all_covs[[k]] <- data.frame(
    grid_id = holdout_grid,
    age_id = dat$age_id[holdout_rows],
    pred = pred,
    truth = truth
  )
  
}

##------------------------------------------------------------
## COMBINE RESULTS
##------------------------------------------------------------

cv_results_all_covs <- bind_rows(cv_results_all_covs)

##------------------------------------------------------------
## METRICS
##------------------------------------------------------------

mae <- mean(
  abs(cv_results_all_covs$pred - cv_results_all_covs$truth),
  na.rm = TRUE
)

bias <- mean(
  cv_results_all_covs$pred - cv_results_all_covs$truth,
  na.rm = TRUE
)

pearson_cor <- cor(
  cv_results_all_covs$pred,
  cv_results_all_covs$truth,
  use = "complete.obs"
)

fit_metrics <- c("mae"= mae, "bias"= bias,"pearson_cor" = pearson_cor)

saveRDS(fit_metrics, file = here(
  "outputs/age_pool_cv_all_covs.rds")
)
