# Sanity check of the scaled simulation: on the SAME scaled datasets, does the Poisson
# model recover the truth (so the data generation is fine) while the current term does not?
# Also tracks the HEV fraction p_HEV_d, to see where the current term's bias enters.
suppressMessages({library(cmdstanr); library(dplyr); library(stringr); library(purrr); library(parallel)})
real <- readRDS("sens_spec_test/results_real.rds")
groups <- real$groups; dat <- real$data; places <- unique(groups$MCH)
X_full <- model.matrix(~ MCH + age_class + preg_status - 1, data = groups)
X <- X_full[, -which(str_detect(colnames(X_full), "MCH"))]
med <- function(fit, par) { s <- real$res[[fit]]$summary; s$median[startsWith(s$variable, paste0(par, "["))] }
mods <- list(current = cmdstan_model("sens_spec_test/joint_mod_fwd_priorsd.stan"),
             poisson = cmdstan_model("sens_spec_test/joint_mod_poisson_priorsd.stan"))
scale <- 5; ifr <- med("pois95", "p"); ph <- med("pois95", "p_HEV_d"); J <- length(ifr)
one <- function(i, ss) {
  set.seed(1000 + i)
  I_mean <- dat$mean_infections * scale; I_sd <- dat$sd_infections * scale
  I_true <- I_mean * exp((I_sd / I_mean) * rnorm(J))
  DcAJ <- rpois(J, I_true * ifr / ph)
  DHEV <- rbinom(J, dat$DAJ, ss * ph + (1 - ss) * (1 - ph))
  map_dfr(names(mods), function(m) {
    d <- list(J = J, L = length(places), M = ncol(X), DHEV = DHEV, DAJ = dat$DAJ, DcAJ = DcAJ,
              I_post_mean = I_mean, I_post_sd = I_sd, X = X, map_to_place = map_dbl(groups$MCH, ~ which(places == .)),
              sens = ss, spec = ss, d_min = 1, beta_sd = 0.5, beta_d_sd = 0.5)
    fit <- mods[[m]]$sample(data = d, chains = 2, parallel_chains = 1, iter_warmup = 750, iter_sampling = 1500,
              seed = i, refresh = 0, adapt_delta = 0.98, max_treedepth = 15, show_messages = FALSE, show_exceptions = FALSE)
    mm <- function(v) apply(posterior::as_draws_matrix(fit$draws(v)), 2, median)
    tibble(sim = i, sens_spec = ss, model = m, age = groups$group, ifr_ratio = mm("p") / ifr, phev_ratio = mm("p_HEV_d") / ph,
           hev_deaths_ratio = mm("deaths_exp_phevd") / (DcAJ * ph))
  })
}
jobs <- expand.grid(i = 1:8, ss = c(0.95, 0.99))
out <- bind_rows(mclapply(seq_len(nrow(jobs)), function(k) one(jobs$i[k], jobs$ss[k]), mc.cores = 10))
saveRDS(out, "sens_spec_test/dmin_scaled_check.rds")
print(as.data.frame(out %>% group_by(sens_spec, model) %>% summarise(across(ends_with("ratio"), ~round(median(.), 2)), .groups = "drop")))
print(as.data.frame(out %>% group_by(sens_spec, model, age) %>% summarise(across(ends_with("ratio"), ~round(median(.), 2)), .groups = "drop")))
