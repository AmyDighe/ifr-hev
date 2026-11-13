# Two no-fit diagnostics for the d >= 1 hypothesis.
suppressMessages(library(dplyr))
real <- readRDS("sens_spec_test/results_real.rds"); dat <- real$data
# 1. Bias vs true expected HEV deaths per group, in the existing simulations (current term)
s05 <- readRDS("sens_spec_test/sim_recovery.rds") %>% filter(model == "agree") %>% mutate(prior_sd = 0.5)
s15 <- readRDS("sens_spec_test/sim_recovery_betasd1.5.rds") %>% filter(model == "agree") %>% mutate(prior_sd = 1.5)
b <- bind_rows(s05, s15) %>% filter(!grepl("negbin", scenario)) %>%
  mutate(exp_d = dat$mean_infections[group] * truth,
         bin = cut(exp_d, c(0, 0.5, 1, 2, 5, Inf), labels = c("<0.5", "0.5-1", "1-2", "2-5", ">5")))
cat("1. Estimated/true IFR by true expected HEV deaths in the group (current term)\n")
print(as.data.frame(b %>% group_by(prior_sd, bin) %>% summarise(n_groups = n_distinct(paste(scenario, group)),
  pregnant_share = round(mean(age == "pregnant"), 2), median_ratio = round(median(med / truth), 2), .groups = "drop")))
cat("\n   Non-pregnant groups only\n")
print(as.data.frame(b %>% filter(age != "pregnant") %>% group_by(prior_sd, bin) %>% summarise(
  n_groups = n_distinct(paste(scenario, group)), median_ratio = round(median(med / truth), 2), .groups = "drop")))
# 2. What expected number of deaths (lambda) does the community term itself prefer?
term <- function(lambda, n, p, dmin = 1) sapply(lambda, function(l) { d <- dmin:n; sum(dbinom(d, n, p) * dpois(d, l)) })
cat("\n2. Lambda maximising sum_{d>=1} Binomial(d | n, p) * Poisson(d | lambda), against the binomial mean n*p\n")
g <- expand.grid(n = c(1, 3, 6, 10, 30, 60), p = c(0.05, 0.1, 0.2))
g$binomial_mean <- g$n * g$p
g$preferred_lambda <- mapply(function(n, p) optimize(function(l) -term(l, n, p), c(1e-4, 100))$minimum, g$n, g$p)
g$ratio <- g$preferred_lambda / g$binomial_mean
print(round(g, 2))
