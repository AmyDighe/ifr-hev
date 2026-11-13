# Semi-synthetic test of the sens/spec handling in the joint deaths-IFR model.
#
# The private mortality data (my_path/model_input/deaths.rds) are NOT used.
# Instead:
#   - infections: real, from outputs/infections.rds
#   - DcAJ: real, recovered exactly from the saved 0.99 fit
#           (deaths_exp_phevd / p_HEV_d)
#   - DAJ, DHEV: SIMULATED. DAJ is made up; DHEV ~ Binomial(DAJ, posterior mean
#           p_HEV_d from the saved 0.99 fit).
#
# Usage: Rscript sens_spec_test/run_test.R   (from the repo root)

suppressMessages({
  library(cmdstanr); library(rstan); library(dplyr); library(stringr); library(purrr)
})

groups <- readRDS("outputs/groups_ifr.rds")
infections_summary <- readRDS("outputs/infections.rds") %>%
  group_by(MCH, age_class, preg_status) %>%
  summarise(mean_infections = mean(infections),
            sd_infections = sd(infections), .groups = "drop")
data <- left_join(groups, infections_summary,
                  by = c("MCH", "age_class", "preg_status"))

old_fit <- readRDS("outputs/joint_fit_age_pooling99.rds")
e <- rstan::extract(old_fit, c("p_HEV_d", "deaths_exp_phevd"))
DcAJ <- as.integer(round(e$deaths_exp_phevd[1, ] / e$p_HEV_d[1, ]))
p_hat <- colMeans(e$p_HEV_d)

set.seed(20261112)
DAJ <- sample(10:60, nrow(groups), replace = TRUE)
DHEV <- rbinom(nrow(groups), DAJ, p_hat)

X_full <- model.matrix(~ MCH + age_class + preg_status - 1, data = groups)
X <- X_full[, -which(str_detect(colnames(X_full), "MCH"))]
places <- unique(groups$MCH)

base_data <- list(
  J = nrow(groups), L = length(places), M = ncol(X),
  DHEV = DHEV, DAJ = DAJ, DcAJ = DcAJ,
  I_post_mean = data$mean_infections, I_post_sd = data$sd_infections,
  X = X, map_to_place = map_dbl(groups$MCH, ~ which(places == .)))

mod_old <- cmdstan_model("code/R/rstan-models/joint_mod_deaths_ifr.stan")
mod_new <- cmdstan_model("code/R/rstan-models/joint_mod_deaths_ifr_fwd.stan")

pars <- c("p", "p_HEV_d", "deaths_exp_phevd", "alpha_mu", "alpha_sigma",
          "alpha_mu_d", "alpha_sigma_d", "beta", "beta_d")

run <- function(label, mod, sens, spec, d_min = NULL) {
  d <- c(base_data, list(sens = sens, spec = spec))
  if (!is.null(d_min)) d$d_min <- d_min
  t0 <- Sys.time()
  fit <- mod$sample(data = d, chains = 4, parallel_chains = 4,
                    iter_warmup = 1000, iter_sampling = 3000, seed = 1,
                    refresh = 0, adapt_delta = 0.98, max_treedepth = 15,
                    show_messages = FALSE, show_exceptions = FALSE)
  s <- fit$summary(pars, mean, median,
                   ~posterior::quantile2(.x, c(0.025, 0.975)),
                   rhat = posterior::rhat, ess_bulk = posterior::ess_bulk)
  dg <- fit$diagnostic_summary(quiet = TRUE)
  out <- list(
    label = label, sens = sens, spec = spec,
    divergent = sum(dg$num_divergent),
    max_treedepth_hits = sum(dg$num_max_treedepth),
    min_bfmi = min(dg$ebfmi),
    max_rhat = max(s$rhat, na.rm = TRUE),
    min_ess = min(s$ess_bulk, na.rm = TRUE),
    mean_stepsize = mean(fit$metadata()$step_size_adaptation),
    minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")),
    summary = s)
  cat(sprintf("%-20s se=%.2f sp=%.2f | div=%4d tree=%d Rhat=%.3f ess=%6.0f step=%.4f | %.1f min\n",
              label, sens, spec, out$divergent, out$max_treedepth_hits,
              out$max_rhat, out$min_ess, out$mean_stepsize, out$minutes))
  out
}

res <- list(
  run("original",            mod_old, 0.99, 0.99),
  run("original",            mod_old, 0.95, 0.95),
  run("forward (d from 1)",  mod_new, 0.99, 0.99, d_min = 1),
  run("forward (d from 1)",  mod_new, 0.95, 0.95, d_min = 1),
  run("forward (d from 0)",  mod_new, 0.95, 0.95, d_min = 0)
)
saveRDS(list(res = res, DAJ = DAJ, DHEV = DHEV, DcAJ = DcAJ, groups = groups),
        "sens_spec_test/results.rds")
