# Simulation-recovery test with WIDER priors on age/pregnancy effects (arg 2 = prior sd).
# Simulate data from known IFRs, fit three formulations, compare recovery.
#   agree   : current term, sum_d Binomial(d | DcAJ, p) * Poisson(d | I * IFR), d from 1
#   poisson : DcAJ ~ Poisson(I * IFR / p)
#   negbin  : DcAJ ~ NegBinomial(I * IFR / p, phi)
# Truth (IFR and true P(HEV | fatal jaundice) per group) is taken from fits to
# the real data; hospital denominators (DAJ) and infections are the real ones.
# Usage: Rscript sens_spec_test/sim_recovery.R [n_sims]   (from repo root)

suppressMessages({library(cmdstanr); library(dplyr); library(stringr); library(purrr); library(parallel)})
n_sims <- as.integer(commandArgs(trailingOnly = TRUE)[1]); if (is.na(n_sims)) n_sims <- 50
real <- readRDS("sens_spec_test/results_real.rds")
groups <- real$groups; dat <- real$data; places <- unique(groups$MCH)
X_full <- model.matrix(~ MCH + age_class + preg_status - 1, data = groups)
X <- X_full[, -which(str_detect(colnames(X_full), "MCH"))]
med <- function(fit, par) { s <- real$res[[fit]]$summary; s$median[startsWith(s$variable, paste0(par, "["))] }
sens <- 0.95; spec <- 0.95

beta_sd <- as.numeric(commandArgs(trailingOnly = TRUE)[2]); if (is.na(beta_sd)) beta_sd <- 1.5
mods <- list(
  agree   = cmdstan_model("sens_spec_test/joint_mod_fwd_priorsd.stan"),
  poisson = cmdstan_model("sens_spec_test/joint_mod_poisson_priorsd.stan"))

scenarios <- list(
  A_truth_agree_fit_poisson_dgp   = list(truth = "fwd95",  phi = Inf),
  B_truth_poisson_fit_poisson_dgp = list(truth = "pois95", phi = Inf),
  C_truth_poisson_fit_negbin_dgp  = list(truth = "pois95", phi = 3))

one_sim <- function(i, sc) {
  set.seed(1000 + i)
  ifr <- med(sc$truth, "p"); ph <- med(sc$truth, "p_HEV_d"); J <- length(ifr)
  I_true <- dat$mean_infections * exp((dat$sd_infections / dat$mean_infections) * rnorm(J))
  mu_all <- I_true * ifr / ph                       # expected community jaundice deaths
  DcAJ <- if (is.finite(sc$phi)) rnbinom(J, mu = mu_all, size = sc$phi) else rpois(J, mu_all)
  DHEV <- rbinom(J, dat$DAJ, sens * ph + (1 - spec) * (1 - ph))
  d <- list(J = J, L = length(places), M = ncol(X), DHEV = DHEV, DAJ = dat$DAJ, DcAJ = DcAJ,
            I_post_mean = dat$mean_infections, I_post_sd = dat$sd_infections, X = X,
            map_to_place = map_dbl(groups$MCH, ~ which(places == .)),
            sens = sens, spec = spec, d_min = 1, beta_sd = beta_sd, beta_d_sd = beta_sd)
  map_dfr(names(mods), function(m) {
    fit <- tryCatch(mods[[m]]$sample(data = d, chains = 2, parallel_chains = 1,
                      iter_warmup = 750, iter_sampling = 1500, seed = i, refresh = 0,
                      adapt_delta = 0.98, max_treedepth = 15,
                      show_messages = FALSE, show_exceptions = FALSE), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    dr <- posterior::as_draws_matrix(fit$draws("p"))
    q <- apply(dr, 2, quantile, c(0.025, 0.5, 0.975))
    tibble(sim = i, model = m, group = seq_len(J), age = groups$group, truth = ifr,
           lo = q[1, ], med = q[2, ], hi = q[3, ],
           div = sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent))
  })
}

out <- imap_dfr(scenarios, function(sc, nm) {
  r <- bind_rows(mclapply(seq_len(n_sims), one_sim, sc = sc, mc.cores = 6))
  cat("done", nm, nrow(r), "rows\n"); mutate(r, scenario = nm)
})
saveRDS(out, sprintf("sens_spec_test/sim_recovery_betasd%s.rds", beta_sd))

summ <- function(df, ...) df %>% group_by(scenario, model, ...) %>% summarise(
  median_ratio = median(med / truth),              # 1 = unbiased
  rmse_log = sqrt(mean(log(med / truth)^2)),
  coverage95 = mean(lo <= truth & truth <= hi),
  ci_fold_width = exp(median(log(hi / lo))),
  .groups = "drop")
cat("\nOverall\n"); print(as.data.frame(summ(out) %>% mutate(across(where(is.numeric), ~round(., 2)))))
cat("\nBy age/pregnancy group\n"); print(as.data.frame(summ(out, age) %>% mutate(across(where(is.numeric), ~round(., 2)))))
cat("\nFits with any divergence\n")
print(as.data.frame(out %>% distinct(scenario, model, sim, div) %>% group_by(scenario, model) %>%
  summarise(fits = n(), with_div = sum(div > 0), .groups = "drop")))
