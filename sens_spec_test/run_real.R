# Refit of script 7a on the REAL mortality data, via cmdstanr (rstan cannot
# load compiled models on this machine). Data prep copied from
# code/R/analysis-scripts/7a_estimate_deaths_IFR.R.
#
# Usage: Rscript sens_spec_test/run_real.R /path/to/data_private   (from repo root)

suppressMessages({library(cmdstanr); library(dplyr); library(stringr); library(purrr)})
my_path <- commandArgs(trailingOnly = TRUE)[1]

deaths <- readRDS(file.path(my_path, "model_input/deaths.rds")) %>%
  rename("preg_status" = maternal) %>%
  mutate(
    age_class = recode(age_class, `(13,44]` = "14-44", `(44,59]` = "45-59", `(59,101]` = "60+"),
    preg_status = recode(preg_status, NO = "not_pregnant", YES = "pregnant"),
    group = if_else(preg_status == "not_pregnant", age_class, "pregnant"))

infections <- readRDS("outputs/infections.rds") %>%
  mutate(group = if_else(preg_status == "not_pregnant", age_class, "pregnant"))
infections_summary <- infections %>%
  group_by(MCH, age_class, preg_status, group) %>%
  summarise(mean_infections = mean(infections), sd_infections = sd(infections), .groups = "drop")

data <- full_join(infections_summary, deaths, by = c("MCH", "group", "age_class", "preg_status"))
groups <- data %>% dplyr::select(-DAJ, -DHEV, -DcAJ, -mean_infections, -sd_infections)
places <- unique(groups$MCH)
stopifnot(identical(as.data.frame(groups[, c("MCH","age_class","preg_status","group")]),
                    as.data.frame(readRDS("outputs/groups_ifr.rds")[, c("MCH","age_class","preg_status","group")])))

X_full <- model.matrix(~ MCH + age_class + preg_status - 1, data = groups)
X <- X_full[, -which(str_detect(colnames(X_full), "MCH"))]

base_data <- list(
  J = nrow(groups), L = length(places), M = ncol(X),
  DHEV = data$DHEV, DAJ = data$DAJ, DcAJ = data$DcAJ,
  I_post_mean = data$mean_infections, I_post_sd = data$sd_infections,
  X = X, map_to_place = map_dbl(groups$MCH, ~ which(places == .)))

mod_old  <- cmdstan_model("code/R/rstan-models/joint_mod_deaths_ifr.stan")
mod_new  <- cmdstan_model("code/R/rstan-models/joint_mod_deaths_ifr_fwd.stan")
mod_pois <- cmdstan_model("sens_spec_test/joint_mod_poisson_dcaj.stan")

pars <- c("p", "p_HEV_d", "deaths_exp_phevd", "alpha_mu", "alpha_sigma",
          "alpha_mu_d", "alpha_sigma_d", "beta", "beta_d")

run <- function(label, mod, sens, spec, d_min = NULL) {
  d <- c(base_data, list(sens = sens, spec = spec))
  if (!is.null(d_min)) d$d_min <- d_min
  fit <- mod$sample(data = d, chains = 4, parallel_chains = 4,
                    iter_warmup = 1000, iter_sampling = 3000, seed = 1,
                    refresh = 0, adapt_delta = 0.98, max_treedepth = 15,
                    show_messages = FALSE, show_exceptions = FALSE)
  s <- fit$summary(pars, mean, median, ~posterior::quantile2(.x, c(0.025, 0.975)),
                   rhat = posterior::rhat, ess_bulk = posterior::ess_bulk)
  dg <- fit$diagnostic_summary(quiet = TRUE)
  out <- list(label = label, sens = sens, spec = spec,
              divergent = sum(dg$num_divergent), max_treedepth_hits = sum(dg$num_max_treedepth),
              min_bfmi = min(dg$ebfmi), max_rhat = max(s$rhat, na.rm = TRUE),
              min_ess = min(s$ess_bulk, na.rm = TRUE),
              mean_stepsize = mean(fit$metadata()$step_size_adaptation), summary = s)
  cat(sprintf("%-14s se=%.2f sp=%.2f | div=%4d tree=%d Rhat=%.3f ess=%6.0f step=%.4f\n",
              label, sens, spec, out$divergent, out$max_treedepth_hits,
              out$max_rhat, out$min_ess, out$mean_stepsize))
  out
}

res <- list(
  orig99 = run("original", mod_old, 0.99, 0.99),
  orig95 = run("original", mod_old, 0.95, 0.95),
  fwd99  = run("forward",  mod_new, 0.99, 0.99, d_min = 1),
  fwd95  = run("forward",  mod_new, 0.95, 0.95, d_min = 1),
  pois99 = run("poisson DcAJ", mod_pois, 0.99, 0.99, d_min = 1),
  pois95 = run("poisson DcAJ", mod_pois, 0.95, 0.95, d_min = 1))
saveRDS(list(res = res, groups = groups, data = data), "sens_spec_test/results_real.rds")

g <- paste(substr(groups$MCH, 1, 4), groups$group)
tab <- function(par, col = "median", f = 1) {
  m <- sapply(res, function(x) { s <- x$summary; s[[col]][startsWith(s$variable, paste0(par, "["))] * f })
  rownames(m) <- g; m }
cat("\nData (DHEV/DAJ hospital, DcAJ community)\n")
print(data.frame(g, DHEV = data$DHEV, DAJ = data$DAJ, raw = round(data$DHEV/data$DAJ, 3), DcAJ = data$DcAJ))
cat("\nIFR median per 10,000\n"); print(round(tab("p", f = 1e4), 2))
cat("\np_HEV_d median (orig = apparent, others = true)\n"); print(round(tab("p_HEV_d"), 3))
cat("\nExpected community HEV deaths, total\n"); print(round(colSums(tab("deaths_exp_phevd", "mean")), 1))
k <- c("alpha_mu", "alpha_sigma", "alpha_mu_d", "alpha_sigma_d", paste0("beta[", 1:3, "]"), paste0("beta_d[", 1:3, "]"))
m <- sapply(res, function(x) { s <- x$summary; i <- match(k, s$variable)
  sprintf("%.2f (%.2f, %.2f)", s$median[i], s$q2.5[i], s$q97.5[i]) }); rownames(m) <- k
cat("\n"); print(colnames(X)); print(m, quote = FALSE)
