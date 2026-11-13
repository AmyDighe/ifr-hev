# Sensitivity/specificity fix and model checks

Exploratory code from a review of `joint_mod_deaths_ifr.stan`. Nothing here
changes the existing analysis scripts or outputs.

All scripts run from the repository root and use **cmdstanr** (not rstan).
Result files (`*.rds`), logs and compiled models are not committed; the
real-data results contain the private death counts, so keep them out of the
public repository.

## The fix

`code/R/rstan-models/joint_mod_deaths_ifr_fwd.stan` parameterises the *true*
probability of HEV given fatal jaundice and applies test error in the forward
direction (`sens * p + (1 - spec) * (1 - p)`) for the hospital data. This
removes the clipping at 0 that caused divergent transitions at
sens = spec = 0.95.

To use it in `7a_estimate_deaths_IFR.R`, point `stan()` at the new file and
add `d_min = 1` to the data list (`d_min = 1` reproduces the original sum over
d = 1..DcAJ; do not use 0, it collapses the IFRs towards zero).

## Scripts

Run in this order; later scripts read `results_real.rds` written by `run_real.R`.

| Script | What it does | Needs private data |
| --- | --- | --- |
| `run_test.R` | Original vs forward model at 0.99 and 0.95 on semi-synthetic data (real infections and DcAJ, simulated hospital counts) | No |
| `compare.R` | Prints comparisons from `run_test.R` results | No |
| `run_real.R <path to data_private>` | Original, forward and Poisson models on the real data at 0.99 and 0.95 | Yes |
| `negbin_real.R` | Poisson vs negative binomial on DcAJ, overdispersion check | via `results_real.rds` |
| `prior_sensitivity_real.R` | Refits with prior sd 0.5, 1, 1.5, 2.5 on age/pregnancy effects | via `results_real.rds` |
| `sim_recovery.R [n_sims]` | Simulation-recovery test of the community term, original priors | via `results_real.rds` |
| `sim_recovery_priors.R [n_sims] [prior_sd]` | Same, with wider effect priors | via `results_real.rds` |
| `dmin_diagnostics.R` | No-fit checks of the d >= 1 lower limit: bias vs expected deaths per group (needs the two simulation result files) and the expected deaths the community term prefers | via simulation results |
| `sim_dmin_test.R [n_sims] [scale]` | Simulation with infections multiplied by `scale` and the same true IFRs, current term, prior sd 0.5 and 1.5 | via `results_real.rds` |
| `dmin_scaled_check.R` | Current vs Poisson term on the same scaled datasets, tracking IFR and HEV fraction | via `results_real.rds` |
| `sim_poisson_hospital.R [n_sims] [inf_scale]` | Poisson term with the hospital sample multiplied by 1, 5 and 25: does the pregnant underestimate shrink with more hospital data? | via `results_real.rds` |
| `poisson_prior_real.R` | Poisson term on the real data with effect prior sd 0.5 and 1.5 | via `results_real.rds` |
| `candidate_real.R` | Candidate model (`joint_mod_candidate.stan`) on the real data: original priors, pregnancy free, pregnant HEV fraction pooled, illustrative pregnancy prior. Sources `candidate_setup.R`, which defines the variants | via `results_real.rds` |
| `candidate_sim.R [n_sims]` | Simulation recovery of the candidate variants against the current model, IFRs and contrasts vs 14-44 | via `candidate_real.rds` |

## Exploratory models

- `joint_mod_poisson_dcaj.stan`: `DcAJ ~ Poisson(I * IFR / p_HEV_d)` in place of the current community term.
- `joint_mod_negbin_dcaj.stan`: negative binomial version of the above.
- `joint_mod_fwd_priorsd.stan`, `joint_mod_poisson_priorsd.stan`: forward and Poisson models with the effect prior sd passed as data (`beta_sd`, `beta_d_sd`).
- `joint_mod_candidate.stan`: candidate replacement. Forward test error, Poisson community term, all priors passed as data, separate covariate matrix for the HEV fraction.
