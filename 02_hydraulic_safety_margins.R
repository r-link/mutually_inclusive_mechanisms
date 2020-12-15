###############################################################################
#
#    Hydraulic safety margins
#
#    Script: Roman Link (roman.link@plant-ecology.de)	
#
###############################################################################

# This script uses the estimated vulnerability curve parameters from script 01
# and the measured mid-day water potentials from August 2018 to compute
# hydraulic safety margins on species level. Vulnerability curve parameters for
# the ring-porous Quercus species could not be measured in the laboratory in
# Würzburg. For this reason, literature estimates from Lobo et al. (2018) are
# used for these species.

# 1. Preparations -------------------------------------------------------------
# if pacman is not installed, install via install.packages("pacman").
# the package manager will take care of the installation of all subsequent 
# packages.
pacman::p_load(readxl, tidyverse)

# Define functions for pooled mean and standard error of the VC measurements
# (classical inverse variance weighting based on fixed effects meta-analytic
# model)

# averages of parameter estimates
pooled_mean <- function(mean, se) { 
  weighted.mean(mean, 1 / se ^ 2)
}

# standard errors of averages of parameter estimates
pooled_se <- function(se) {
  sqrt( 1 / sum( 1 / se ^ 2))
}

# 2. Load water potential data -----------------------------------------------
# load leaf water potentials from August 2020 and aggregate on species level
leaf_psi <- read_xlsx("data/IDENT_leaf_psi_2020-08.xlsx", sheet = 2) %>% 
  group_by(species) %>% 
  summarize(psimd = mean(-psi_md),
            se_psimd = sd(psi_md) / length(psi_md))

# 3. Load vulnerability curve data -------------------------------------------
# VC parameters for the long-vesseled Quercus species from Lobo et al. (2018) 
querc <- tibble(
  species  = c("QURO", "QURU"),
  slope    = c(67.07, 24.52),
  p50      = c(-4.74, -4.43),
  p12      = c(-3.81, -2.14),
  p88      = c(-5.66, -6.72),
  se_slope = c(14.09, 4.08),
  se_p50   = c(0.09, 0.25),
  se_p12   = c(0.20, 0.48),
  se_p88   = c(0.09, 0.29)
)

# load VC parameters from the curves measured for the IDENT dataset
VC <- read_csv("output/IDENT_VC_parameters.csv") %>% 
  select(-lower, -upper) %>%  # remove confidence intervals
  pivot_wider(names_from = term,  # convert to a wide table
              values_from = c(est, se)) %>%  
  group_by(species) %>%  # group by species
  # calculate aggregated estimates and standard errors
  summarize(slope = pooled_mean(est_slope, se_slope),
            p50   = pooled_mean(est_p50, se_p50),
            p12   = pooled_mean(est_p12, se_p12),
            p88   = pooled_mean(est_p88, se_p88),
            se_slope = pooled_se(se_slope),
            se_p50   = pooled_se(se_p50),
            se_p12   = pooled_se(se_p12),
            se_p88   = pooled_se(se_p88)) %>% 
  bind_rows(querc) # merge with Quercus data
VC

# combine leaf water potential and VC to calculate safety margins 
# (using P88 as a lethal threshold for angiosperms and p50 for conifers)
HSM <- full_join(leaf_psi, VC) %>% 
  mutate(
    conifer = species %in% c("LADE", "LALA", "PIAB", "PIPU", "PIST", "PISY"),
    HSM     = ifelse(conifer, psimd - p50, psimd - p88),
    se_HSM  = ifelse(conifer, sqrt(se_psimd^2 + se_p50^2), sqrt(se_psimd^2 + se_p88^2))
  ) %>% 
  select(-conifer)

# export hydraulic safety margins
write_csv(HSM, "output/IDENT_hydraulic_safety_margins.csv")
