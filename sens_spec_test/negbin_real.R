# Is there enough information in the real data for a negative binomial on DcAJ?
# Fits Poisson and NB versions, compares the overdispersion posterior with its
# prior, and runs a posterior predictive check on the Poisson fit.
suppressMessages({library(cmdstanr); library(dplyr); library(stringr); library(purrr)})
real <- readRDS("sens_spec_test/results_real.rds")
groups <- real$groups; dat <- real$data; places <- unique(groups$MCH)
X_full <- model.matrix(~ MCH + age_class + preg_status - 1, data = groups)
X <- X_full[, -which(str_detect(colnames(X_full), "MCH"))]
mk <- function(s) list(J = nrow(groups), L = length(places), M = ncol(X), DHEV = dat$DHEV, DAJ = dat$DAJ,
  DcAJ = dat$DcAJ, I_post_mean = dat$mean_infections, I_post_sd = dat$sd_infections, X = X,
  map_to_place = map_dbl(groups$MCH, ~ which(places == .)), sens = s, spec = s, d_min = 1)
mp <- cmdstan_model("sens_spec_test/joint_mod_poisson_dcaj.stan")
mn <- cmdstan_model("sens_spec_test/joint_mod_negbin_dcaj.stan")
g <- paste(substr(groups$MCH, 1, 4), groups$group)
for (s in c(0.99, 0.95)) {
  f <- function(m) m$sample(data = mk(s), chains = 4, parallel_chains = 4, iter_warmup = 1000,
    iter_sampling = 3000, seed = 1, refresh = 0, adapt_delta = 0.98, max_treedepth = 15,
    show_messages = FALSE, show_exceptions = FALSE)
  fp <- f(mp); fn <- f(mn)
  cat(sprintf("\n===== sens = spec = %.2f =====\ndivergences: poisson %d, negbin %d\n", s,
      sum(fp$diagnostic_summary(quiet = TRUE)$num_divergent), sum(fn$diagnostic_summary(quiet = TRUE)$num_divergent)))
  isp <- as.numeric(posterior::as_draws_matrix(fn$draws("inv_sqrt_phi")))
  prior <- abs(rnorm(1e5))
  cat("inv_sqrt_phi (0 = Poisson): posterior median (95% CrI)", round(quantile(isp, c(.5, .025, .975)), 2),
      "| prior", round(quantile(prior, c(.5, .025, .975)), 2), "\n")
  cat("phi: posterior median (95% CrI)", round(quantile(1 / isp^2, c(.5, .025, .975)), 1), "\n")
  cat("P(inv_sqrt_phi < 0.2), i.e. close to Poisson: posterior", round(mean(isp < 0.2), 3), "prior", round(mean(prior < 0.2), 3), "\n")
  # posterior predictive check of the Poisson fit: variance/mean discrepancy
  mu <- posterior::as_draws_matrix(fp$draws("I_latent")) * posterior::as_draws_matrix(fp$draws("p")) /
        posterior::as_draws_matrix(fp$draws("p_HEV_d"))
  mu <- matrix(as.numeric(mu), nrow = nrow(mu))
  chi <- function(y, m) rowSums((y - m)^2 / m)
  yrep <- matrix(rpois(length(mu), mu), nrow = nrow(mu))
  obs <- chi(matrix(dat$DcAJ, nrow(mu), ncol(mu), byrow = TRUE), mu)
  cat("Poisson fit PPC: chi-square discrepancy, observed median", round(median(obs), 1), "replicated median",
      round(median(chi(yrep, mu)), 1), "| posterior predictive p =", round(mean(chi(yrep, mu) >= obs), 3), "\n")
  q <- function(fit) { d <- posterior::as_draws_matrix(fit$draws("p")); apply(d, 2, quantile, c(.5, .025, .975)) * 1e4 }
  a <- q(fp); b <- q(fn)
  print(data.frame(group = g, DcAJ = dat$DcAJ, fitted_pois = round(apply(mu, 2, median), 1),
    IFR_pois = sprintf("%.2f (%.2f-%.2f)", a[1, ], a[2, ], a[3, ]),
    IFR_negbin = sprintf("%.2f (%.2f-%.2f)", b[1, ], b[2, ], b[3, ])))
  cat("median ratio of IFR CrI width, negbin/poisson:", round(median((b[3, ] / b[2, ]) / (a[3, ] / a[2, ])), 2),
      "| median ratio of IFR medians:", round(median(b[1, ] / a[1, ]), 2), "\n")
}
