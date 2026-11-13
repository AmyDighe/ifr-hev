# Poisson community term on the real data with effect prior sd 0.5 and 1.5:
# how do the pregnancy effects on IFR (beta[3]) and on the HEV fraction (beta_d[3]) move?
suppressMessages({library(cmdstanr); library(dplyr); library(stringr); library(purrr)})
real <- readRDS("sens_spec_test/results_real.rds")
groups <- real$groups; dat <- real$data; places <- unique(groups$MCH)
X_full <- model.matrix(~ MCH + age_class + preg_status - 1, data = groups)
X <- X_full[, -which(str_detect(colnames(X_full), "MCH"))]
mod <- cmdstan_model("sens_spec_test/joint_mod_poisson_priorsd.stan")
for (psd in c(0.5, 1.5)) {
  d <- list(J = nrow(groups), L = length(places), M = ncol(X), DHEV = dat$DHEV, DAJ = dat$DAJ, DcAJ = dat$DcAJ,
            I_post_mean = dat$mean_infections, I_post_sd = dat$sd_infections, X = X,
            map_to_place = map_dbl(groups$MCH, ~ which(places == .)),
            sens = 0.95, spec = 0.95, d_min = 1, beta_sd = psd, beta_d_sd = psd)
  fit <- mod$sample(data = d, chains = 4, parallel_chains = 4, iter_warmup = 1000, iter_sampling = 3000, seed = 1,
                    refresh = 0, adapt_delta = 0.98, max_treedepth = 15, show_messages = FALSE, show_exceptions = FALSE)
  cat("\nprior sd", psd, " divergences", sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent), "\n")
  print(as.data.frame(fit$summary(c("beta", "beta_d"), median, ~quantile(.x, c(0.025, 0.975)))), digits = 2)
  dr <- posterior::as_draws_matrix(fit$draws(c("p", "p_HEV_d")))
  pg <- which(groups$group == "pregnant")
  s <- function(v) round(quantile(rowMeans(dr[, paste0(v, "[", pg, "]")]), c(0.5, 0.025, 0.975)), 5)
  cat("pregnant mean HEV fraction:", s("p_HEV_d"), "\npregnant mean IFR:", s("p"), "\n")
  cat("pregnancy OR on IFR:", round(exp(quantile(posterior::as_draws_matrix(fit$draws("beta[3]")), c(0.5, 0.025, 0.975))), 2), "\n")
}
