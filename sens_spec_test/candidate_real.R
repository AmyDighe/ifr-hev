# Candidate model on the real data, four variants, compared with the current model (fwd95).
# Usage: Rscript sens_spec_test/candidate_real.R   (from repo root, after run_real.R)
source("sens_spec_test/candidate_setup.R")
mod <- cmdstan_model("sens_spec_test/joint_mod_candidate.stan")
pg <- groups$group == "pregnant"
fits <- map(set_names(variants), function(v) {
  d <- c(base_list(dat$DHEV, dat$DAJ, dat$DcAJ), variant(v))
  fit <- mod$sample(data = d, chains = 4, parallel_chains = 4, iter_warmup = 1000, iter_sampling = 3000, seed = 1,
                    refresh = 0, adapt_delta = 0.98, max_treedepth = 15, show_messages = FALSE, show_exceptions = FALSE)
  dr <- posterior::as_draws_matrix(fit$draws(c("p", "p_HEV_d", "beta", "beta_d", "alpha_mu", "alpha_sigma",
                                               "alpha_mu_d", "alpha_sigma_d", "deaths_exp_ifr", "DcAJ_rep", "DHEV_rep")))
  list(div = sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent),
       rhat = max(fit$summary(c("p", "p_HEV_d"), posterior::rhat)[[2]]), draws = dr)
})
q <- function(x) quantile(x, c(0.5, 0.025, 0.975))
fmt <- function(x, d = 2) sprintf("%s (%s-%s)", signif(x[1], d), signif(x[2], d), signif(x[3], d))
col <- function(dr, v) dr[, startsWith(colnames(dr), paste0(v, "[")), drop = FALSE]
tab <- imap_dfr(fits, function(f, v) { dr <- f$draws; p <- col(dr, "p"); ph <- col(dr, "p_HEV_d"); b <- col(dr, "beta")
  grp <- function(m, g) rowMeans(m[, groups$group == g, drop = FALSE])
  rep_c <- col(dr, "DcAJ_rep"); rep_h <- col(dr, "DHEV_rep")
  tibble(variant = v, divergent = f$div, rhat = round(f$rhat, 3),
         ifr_1444_per100k = fmt(q(grp(p, "14-44")) * 1e5), ifr_4559 = fmt(q(grp(p, "45-59")) * 1e5),
         ifr_60 = fmt(q(grp(p, "60+")) * 1e5), ifr_preg = fmt(q(grp(p, "pregnant")) * 1e5),
         or_4559 = fmt(q(exp(b[, 1]))), or_60 = fmt(q(exp(b[, 2]))), or_preg = fmt(q(exp(b[, preg_col]))),
         hevfrac_1444 = fmt(q(grp(ph, "14-44"))), hevfrac_preg = fmt(q(grp(ph, "pregnant"))),
         hev_deaths_total = fmt(q(rowSums(col(dr, "deaths_exp_ifr"))), 3),
         alpha_sigma = fmt(q(dr[, "alpha_sigma"])),
         ppc_comm_preg = round(mean(rowSums(rep_c[, pg]) >= sum(dat$DcAJ[pg])), 2),
         ppc_hosp_preg = round(mean(rowSums(rep_h[, pg]) >= sum(dat$DHEV[pg])), 2),
         ppc_comm_var = round(mean(apply(rep_c, 1, var) >= var(dat$DcAJ)), 2)) })
print(as.data.frame(t(tab)), quote = FALSE)
# current model for reference (fwd95 from run_real.R)
s <- real$res[["fwd95"]]$summary
cur <- function(v, g) signif(mean(s$median[startsWith(s$variable, paste0(v, "["))][groups$group == g]), 2)
cat("\nCurrent model (fwd95), mean of group medians, IFR per 100k: ",
    paste(c("14-44", "45-59", "60+", "pregnant"), sapply(c("14-44", "45-59", "60+", "pregnant"), function(g) cur("p", g) * 1e5), collapse = "; "), "\n")
print(as.data.frame(s[grepl("^beta\\[", s$variable), c("variable", "median", "q2.5", "q97.5")]))
# per-group medians of the free variant = truth for the simulation
truth <- map(fits, ~ list(ifr = apply(col(.x$draws, "p"), 2, median), ph = apply(col(.x$draws, "p_HEV_d"), 2, median)))
saveRDS(list(table = tab, truth = truth), "sens_spec_test/candidate_real.rds")
