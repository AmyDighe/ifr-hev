# Why does the Poisson community term underestimate IFR in pregnant women?
# Poisson DGP and Poisson model at the real community sample size; the number of
# hospital jaundice deaths tested (DAJ) is multiplied by `hosp_scale`. Tracks both
# IFR and the HEV fraction p_HEV_d, with posterior medians and posterior means.
# `inf_scale` multiplies infections (and so community deaths) as well; default 1.
# Usage: Rscript sens_spec_test/sim_poisson_hospital.R [n_sims] [inf_scale]   (from repo root)
suppressMessages({library(cmdstanr); library(dplyr); library(stringr); library(purrr); library(parallel)})
n_sims <- as.integer(commandArgs(trailingOnly = TRUE)[1]); if (is.na(n_sims)) n_sims <- 40
inf_scale <- as.numeric(commandArgs(trailingOnly = TRUE)[2]); if (is.na(inf_scale)) inf_scale <- 1
real <- readRDS("sens_spec_test/results_real.rds")
groups <- real$groups; dat <- real$data
dat$mean_infections <- dat$mean_infections * inf_scale; dat$sd_infections <- dat$sd_infections * inf_scale
 places <- unique(groups$MCH)
X_full <- model.matrix(~ MCH + age_class + preg_status - 1, data = groups)
X <- X_full[, -which(str_detect(colnames(X_full), "MCH"))]
med <- function(fit, par) { s <- real$res[[fit]]$summary; s$median[startsWith(s$variable, paste0(par, "["))] }
mod <- cmdstan_model("sens_spec_test/joint_mod_poisson_priorsd.stan")
ifr <- med("pois95", "p"); ph <- med("pois95", "p_HEV_d"); J <- length(ifr); ss <- 0.95
one <- function(i, hs, psd) {
  set.seed(1000 + i)
  I_true <- dat$mean_infections * exp((dat$sd_infections / dat$mean_infections) * rnorm(J))
  DcAJ <- rpois(J, I_true * ifr / ph)
  DAJ <- dat$DAJ * hs
  DHEV <- rbinom(J, DAJ, ss * ph + (1 - ss) * (1 - ph))
  d <- list(J = J, L = length(places), M = ncol(X), DHEV = DHEV, DAJ = DAJ, DcAJ = DcAJ,
            I_post_mean = dat$mean_infections, I_post_sd = dat$sd_infections, X = X,
            map_to_place = map_dbl(groups$MCH, ~ which(places == .)),
            sens = ss, spec = ss, d_min = 1, beta_sd = psd, beta_d_sd = psd)
  fit <- tryCatch(mod$sample(data = d, chains = 2, parallel_chains = 1, iter_warmup = 750, iter_sampling = 1500,
            seed = i, refresh = 0, adapt_delta = 0.98, max_treedepth = 15,
            show_messages = FALSE, show_exceptions = FALSE), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  f <- function(v, fn) apply(posterior::as_draws_matrix(fit$draws(v)), 2, fn)
  q <- apply(posterior::as_draws_matrix(fit$draws("p")), 2, quantile, c(0.025, 0.975))
  tibble(sim = i, hosp_scale = hs, prior_sd = psd, age = groups$group,
         ifr_med = f("p", median) / ifr, ifr_mean = f("p", mean) / ifr,
         ph_med = f("p_HEV_d", median) / ph, ph_mean = f("p_HEV_d", mean) / ph,
         cover = q[1, ] <= ifr & ifr <= q[2, ], fold = q[2, ] / q[1, ],
         hosp_pos_true = rep(sum((DAJ * ph)[groups$group == "pregnant"]), J))
}
jobs <- expand.grid(i = seq_len(n_sims), hs = if (inf_scale == 1) c(1, 5, 25) else 25, psd = c(0.5, 1.5))
out <- bind_rows(mclapply(seq_len(nrow(jobs)), function(k) one(jobs$i[k], jobs$hs[k], jobs$psd[k]), mc.cores = 10))
saveRDS(out, sprintf("sens_spec_test/sim_poisson_hospital_inf%g.rds", inf_scale))
print(as.data.frame(out %>% group_by(prior_sd, hosp_scale, age) %>% summarise(
  n = n() / 6, ifr_median_of_medians = median(ifr_med), ifr_mean_of_means = mean(ifr_mean),
  ph_median_of_medians = median(ph_med), ph_mean_of_means = mean(ph_mean),
  coverage = mean(cover), fold = exp(median(log(fold))), .groups = "drop") %>%
  mutate(across(where(is.numeric), ~round(., 2)))))
