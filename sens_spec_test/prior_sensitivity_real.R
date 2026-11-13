# Sensitivity of the real-data results to the prior sd on age/pregnancy effects.
# Usage: Rscript sens_spec_test/prior_sensitivity_real.R   (from repo root; needs results_real.rds)
suppressMessages({library(cmdstanr); library(dplyr); library(stringr); library(purrr)})
real <- readRDS("sens_spec_test/results_real.rds")
groups <- real$groups; dat <- real$data; places <- unique(groups$MCH)
X_full <- model.matrix(~ MCH + age_class + preg_status - 1, data = groups)
X <- X_full[, -which(str_detect(colnames(X_full), "MCH"))]
mods <- list(current = cmdstan_model("sens_spec_test/joint_mod_fwd_priorsd.stan"),
             poisson = cmdstan_model("sens_spec_test/joint_mod_poisson_priorsd.stan"))
f3 <- function(x, k = 1) sprintf("%.*f (%.*f-%.*f)", k, median(x), k, quantile(x, .025), k, quantile(x, .975))
out <- list()
for (s in c(0.99, 0.95)) for (m in names(mods)) for (bs in c(0.5, 1, 1.5, 2.5)) {
  d <- list(J = nrow(groups), L = length(places), M = ncol(X), DHEV = dat$DHEV, DAJ = dat$DAJ, DcAJ = dat$DcAJ,
            I_post_mean = dat$mean_infections, I_post_sd = dat$sd_infections, X = X,
            map_to_place = map_dbl(groups$MCH, ~ which(places == .)), sens = s, spec = s, d_min = 1,
            beta_sd = bs, beta_d_sd = bs)
  fit <- mods[[m]]$sample(data = d, chains = 4, parallel_chains = 4, iter_warmup = 1000, iter_sampling = 3000,
            seed = 1, refresh = 0, adapt_delta = 0.98, max_treedepth = 15, show_messages = FALSE, show_exceptions = FALSE)
  dr <- posterior::as_draws_df(fit$draws(c("alpha_mu", "beta", "beta_d", "alpha_mu_d", "deaths_exp_phevd")))
  ifr <- function(b) plogis(dr$alpha_mu + b) * 1e4
  tot <- rowSums(as.matrix(dr[, grep("^deaths_exp_phevd", names(dr))]))
  out[[length(out) + 1]] <- tibble(sens_spec = s, model = m, prior_sd = bs,
    div = sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent),
    IFR_14_44 = f3(ifr(0)), IFR_45_59 = f3(ifr(dr$`beta[1]`)), IFR_60 = f3(ifr(dr$`beta[2]`)),
    IFR_preg = f3(ifr(dr$`beta[3]`)), OR_preg = f3(exp(dr$`beta[3]`)), OR_60 = f3(exp(dr$`beta[2]`)),
    pHEV_14_44 = f3(plogis(dr$alpha_mu_d), 2), pHEV_60 = f3(plogis(dr$alpha_mu_d + dr$`beta_d[2]`), 2),
    pHEV_preg = f3(plogis(dr$alpha_mu_d + dr$`beta_d[3]`), 2), HEV_deaths = f3(tot, 0))
}
out <- bind_rows(out); saveRDS(out, "sens_spec_test/prior_sensitivity_real.rds")
options(width = 250); print(as.data.frame(out))
