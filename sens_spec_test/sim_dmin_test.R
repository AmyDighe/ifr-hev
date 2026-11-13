# Does the d >= 1 lower limit of the community term cause the pregnancy overshoot?
# Intervention test: multiply infections (and hence expected deaths) by `scale`
# with the SAME true IFRs. If forcing >= 1 HEV death per group is the cause, the
# overshoot under wide priors should vanish once every group expects several deaths.
# Baseline (scale = 1) results are in sim_recovery.rds / sim_recovery_betasd1.5.rds.
# Usage: Rscript sens_spec_test/sim_dmin_test.R [n_sims] [scale]   (from repo root)
suppressMessages({library(cmdstanr); library(dplyr); library(stringr); library(purrr); library(parallel)})
a <- commandArgs(trailingOnly = TRUE)
n_sims <- as.integer(a[1]); if (is.na(n_sims)) n_sims <- 50
scale <- as.numeric(a[2]); if (is.na(scale)) scale <- 10
real <- readRDS("sens_spec_test/results_real.rds")
groups <- real$groups; dat <- real$data; places <- unique(groups$MCH)
X_full <- model.matrix(~ MCH + age_class + preg_status - 1, data = groups)
X <- X_full[, -which(str_detect(colnames(X_full), "MCH"))]
med <- function(fit, par) { s <- real$res[[fit]]$summary; s$median[startsWith(s$variable, paste0(par, "["))] }
sens <- 0.95; spec <- 0.95
mod <- cmdstan_model("sens_spec_test/joint_mod_fwd_priorsd.stan")
scenarios <- list(A_truth_agree_fit = "fwd95", B_truth_poisson_fit = "pois95")

one_sim <- function(i, truth) {
  set.seed(1000 + i)
  ifr <- med(truth, "p"); ph <- med(truth, "p_HEV_d"); J <- length(ifr)
  I_mean <- dat$mean_infections * scale; I_sd <- dat$sd_infections * scale
  I_true <- I_mean * exp((I_sd / I_mean) * rnorm(J))
  DcAJ <- rpois(J, I_true * ifr / ph)
  DHEV <- rbinom(J, dat$DAJ, sens * ph + (1 - spec) * (1 - ph))
  map_dfr(c(0.5, 1.5), function(bs) {
    d <- list(J = J, L = length(places), M = ncol(X), DHEV = DHEV, DAJ = dat$DAJ, DcAJ = DcAJ,
              I_post_mean = I_mean, I_post_sd = I_sd, X = X,
              map_to_place = map_dbl(groups$MCH, ~ which(places == .)),
              sens = sens, spec = spec, d_min = 1, beta_sd = bs, beta_d_sd = bs)
    fit <- tryCatch(mod$sample(data = d, chains = 2, parallel_chains = 1, iter_warmup = 750,
              iter_sampling = 1500, seed = i, refresh = 0, adapt_delta = 0.98, max_treedepth = 15,
              show_messages = FALSE, show_exceptions = FALSE), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    q <- apply(posterior::as_draws_matrix(fit$draws("p")), 2, quantile, c(0.025, 0.5, 0.975))
    tibble(sim = i, prior_sd = bs, group = seq_len(J), age = groups$group, truth = ifr,
           exp_hev_deaths = dat$mean_infections * scale * ifr, lo = q[1, ], med = q[2, ], hi = q[3, ])
  })
}
out <- imap_dfr(scenarios, function(truth, nm) {
  r <- bind_rows(mclapply(seq_len(n_sims), one_sim, truth = truth, mc.cores = 10))
  cat("done", nm, nrow(r), "rows\n"); mutate(r, scenario = nm)
})
saveRDS(out, sprintf("sens_spec_test/sim_dmin_scale%s.rds", scale))
print(as.data.frame(out %>% group_by(scenario, prior_sd, age) %>% summarise(
  exp_hev_deaths = round(median(exp_hev_deaths), 1), median_ratio = round(median(med / truth), 2),
  coverage95 = round(mean(lo <= truth & truth <= hi), 2), .groups = "drop")))
