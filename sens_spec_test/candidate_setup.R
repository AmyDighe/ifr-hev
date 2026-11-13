# Shared set-up for the candidate model: data, design matrices and the prior/pregnancy variants.
suppressMessages({library(cmdstanr); library(dplyr); library(stringr); library(purrr)})
real <- readRDS("sens_spec_test/results_real.rds")
groups <- real$groups; dat <- real$data; places <- unique(groups$MCH)
X_full <- model.matrix(~ MCH + age_class + preg_status - 1, data = groups)
X <- X_full[, -which(str_detect(colnames(X_full), "MCH"))]
preg_col <- which(str_detect(colnames(X), "preg")); M <- ncol(X); J <- nrow(groups)
stopifnot(length(preg_col) == 1)
base_list <- function(DHEV, DAJ, DcAJ, ss = 0.95) list(
  J = J, L = length(places), M = M, DHEV = DHEV, DAJ = DAJ, DcAJ = DcAJ,
  I_post_mean = dat$mean_infections, I_post_sd = dat$sd_infections, X = X,
  map_to_place = map_dbl(groups$MCH, ~ which(places == .)), sens = ss, spec = ss)
# Candidate priors: IFR centred on 1 in 3000 with a 95% range of about 1 in 150,000 to 1 in 60;
# site spread and covariate effects with sd 1 (the originals were 0.5).
variant <- function(name) {
  v <- list(X_d = X, alpha_mu_mean = -8, alpha_mu_sd = 2, alpha_sigma_sd = 1,
            alpha_mu_d_mean = -2, alpha_mu_d_sd = 1, alpha_sigma_d_sd = 1,
            beta_mean = rep(0, M), beta_sd = rep(1, M), beta_d_mean = rep(0, M), beta_d_sd = rep(1, M))
  if (name == "original_priors") {       # Poisson term, Amy's priors
    v$alpha_mu_mean <- -6; v$alpha_sigma_sd <- 0.5; v$alpha_sigma_d_sd <- 0.5
    v$beta_sd <- rep(0.5, M); v$beta_d_sd <- rep(0.5, M)
  } else if (name == "pooled") {         # pregnant HEV fraction = non-pregnant 14-44
    v$X_d[, preg_col] <- 0
  } else if (name == "informative") {    # ILLUSTRATIVE prior on the pregnancy odds ratio for IFR:
    v$beta_mean[preg_col] <- log(5); v$beta_sd[preg_col] <- 0.5   # median 5, 95% range 1.9 to 13
  } else stopifnot(name == "free")
  v
}
variants <- c("original_priors", "free", "pooled", "informative")
