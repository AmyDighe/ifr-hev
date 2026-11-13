# Simulation recovery for the candidate model variants against the current model.
# Truth = per-group posterior medians from a real-data fit of the candidate
# ("free": pregnancy effect large; "original_priors": pregnancy effect modest).
# Data: Poisson community jaundice deaths, binomial hospital results, real sample sizes.
# Usage: Rscript sens_spec_test/candidate_sim.R [n_sims]   (after candidate_real.R)
source("sens_spec_test/candidate_setup.R"); library(parallel)
n_sims <- as.integer(commandArgs(trailingOnly = TRUE)[1]); if (is.na(n_sims)) n_sims <- 40
truths <- readRDS("sens_spec_test/candidate_real.rds")$truth[c("free", "original_priors")]
cand <- cmdstan_model("sens_spec_test/joint_mod_candidate.stan")
cur <- cmdstan_model("code/R/rstan-models/joint_mod_deaths_ifr_fwd.stan")
ref <- match(paste(groups$MCH, "14-44"), paste(groups$MCH, groups$group))   # same-site 14-44 group
contrast <- function(p, g) { k <- which(groups$group == g); rowMeans(log(p[, k, drop = FALSE] / p[, ref[k], drop = FALSE])) }
one <- function(i, tn) {
  set.seed(1000 + i); tr <- truths[[tn]]
  I_true <- dat$mean_infections * exp((dat$sd_infections / dat$mean_infections) * rnorm(J))
  DcAJ <- rpois(J, I_true * tr$ifr / tr$ph)
  DHEV <- rbinom(J, dat$DAJ, 0.95 * tr$ph + 0.05 * (1 - tr$ph))
  b <- base_list(DHEV, dat$DAJ, DcAJ)
  map_dfr(c("current", variants), function(m) {
    fit <- tryCatch((if (m == "current") cur else cand)$sample(
      data = if (m == "current") c(b, list(d_min = 1)) else c(b, variant(m)),
      chains = 2, parallel_chains = 1, iter_warmup = 750, iter_sampling = 1500, seed = i, refresh = 0,
      adapt_delta = 0.98, max_treedepth = 15, show_messages = FALSE, show_exceptions = FALSE), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    p <- unclass(posterior::as_draws_matrix(fit$draws("p")))
    q <- apply(p, 2, quantile, c(0.025, 0.5, 0.975))
    grp <- tibble(kind = "ifr", what = groups$group, truth = tr$ifr, lo = q[1, ], med = q[2, ], hi = q[3, ], mean = colMeans(p))
    con <- map_dfr(c("45-59", "60+", "pregnant"), function(g) { x <- exp(contrast(p, g)); k <- which(groups$group == g)
      tibble(kind = "ratio_vs_1444", what = g, truth = exp(mean(log(tr$ifr[k] / tr$ifr[ref[k]]))),
             lo = quantile(x, 0.025), med = median(x), hi = quantile(x, 0.975), mean = mean(x)) })
    bind_rows(grp, con) %>% mutate(sim = i, truth_set = tn, model = m,
                                   div = sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent))
  })
}
jobs <- expand.grid(i = seq_len(n_sims), tn = names(truths), stringsAsFactors = FALSE)
out <- bind_rows(mclapply(seq_len(nrow(jobs)), function(k) one(jobs$i[k], jobs$tn[k]), mc.cores = 10))
saveRDS(out, "sens_spec_test/candidate_sim.rds")
res <- out %>% mutate(model = factor(model, c("current", variants))) %>% group_by(truth_set, kind, what, model) %>% summarise(
  fits = n_distinct(sim), true_value = signif(exp(mean(log(truth))), 2), median_ratio = median(med / truth), mean_ratio = mean(mean / truth),
  rmse_log = sqrt(mean(log(med / truth)^2)), coverage = mean(lo <= truth & truth <= hi),
  fold_width = exp(median(log(hi / lo))), .groups = "drop") %>% mutate(across(where(is.numeric), ~round(., 2)))
options(width = 200); print(as.data.frame(res))
print(as.data.frame(out %>% distinct(truth_set, model, sim, div) %>% group_by(truth_set, model) %>% summarise(fits = n(), with_div = sum(div > 0), .groups = "drop")))
