###############################################################################
#
#    IDENT modeling with Stan
#
#    Script: Roman Link  (roman.link@plant-ecology.de)	
#
###############################################################################

# This script fits the hierarchical Bayesian model for tree mortality defined in
# "stan/neighbourhood_mortality_model.stan" and performs basic model checking.

# 1. Load packages ------------------------------------------------------------
# if pacman is not installed, install via install.packages("pacman").
# the package manager will take care of the installation of all subsequent 
# packages.
pacman::p_load(tidyverse, readxl, brms, shinystan, rstan, loo, tidybayes)

# rstan option settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# 2. Read and process raw data ------------------------------------------------
# mortality dataset
mort  <- read_xlsx("data/IDENT_mortality.xlsx", sheet = 2)

# neighbourhood matrix
neigh <- read_csv("output/IDENT_neighbourhood_matrix.csv") 

# neighbourhood-specific beetle counts
beet  <- read_csv("output/IDENT_neighbours_beetles.csv")

# hydraulic safety margins
hsm   <- read_csv("output/IDENT_hydraulic_safety_margins.csv") 

# relative sugar contents
relsug <- read_xlsx("data/IDENT_relative_sugar_contents.xlsx", sheet = 2) %>% 
  select(species, relsug = diffSugOct) %>% 
  na.omit()

# join and prepare datasets
data <- full_join(mort, neigh) %>% 
  left_join(beet) %>% 
  # convert mortality to integer (resprouted trees counted as functionally dead)
  mutate(mort = as.numeric(mortality_2019 != "alive")) %>% 
  # filter unnecessary rows
  filter(position == "center",   # only use central trees as a response
         is.na(mortality_2017),  # remove trees that were dead before the experiments
         !is.na(height_2018),    # remove trees without reliable height data
         height_2018 >= 30) %>%  # remove trees with a height below 30 cm (likely resprouts)
  # prepare predictor variables
  mutate(
    height = as.numeric(scale(log(height_2018))),  # height -> scaled, centered and log-transformed
    beetles_neigh = as.numeric(beetles_neigh > 0), # dichotomize beetles_neigh (almost no trees with more than 1-2 infested neighbours)
    plotspec = paste0(plot_cn, "_", species))      # define index variable for plot-species combinations
data


# 3. Construct dataset for model fitting in Stan using brms constructors ------

# use brms::make_standata to construct list of data for stan model fitting
stan_data <- make_standata(
  mort ~ height + beetles + beetles_neigh + fert + # main effects
    (height + beetles_neigh + fert | species) +     # species-wise effects
    (ACPL + ACSA + BEPA + BEPE + LADE + LALA + 
       PIAB + PIPU + PIST + PISY + QURO + QURU + 0 || species) + # neighbourhood effects
    (1 | plotspec),                                  # design effects
  data = data)


# position matrix for neighbour effects 
# (matrix of index variables used within Stan to handle within-pair correlations)
pos_mat <-  matrix(0, nrow = 12, ncol = 12)
k <- 1
for (i in 1:11){
  for (j in (i + 1):12){
    pos_mat[i, j] <- k;
    pos_mat[j, i] <- k;
    k <- k +1;
  }
}
pos_mat

# construct HSM observations with properly accounted uncertainty
hsm1 <- mutate(hsm,
               hsm_scaled = as.numeric(scale(HSM)),
               hsm_se_scaled = se_HSM/sd(HSM)) 

# append additional data for model fitting
data_fit <- c(stan_data,
              list(pos_mat = pos_mat,          # position matrix for neighbour effects
                   N_pairs = max(pos_mat),     # number of neighbour species pairs
                   N_spec  = 12,               # number of species
                   hsm_obs = hsm1$hsm_scaled,  # observed HSM values
                   sd_hsm  = hsm1$hsm_se_scaled, # standard error of HSM values
                   N_nsc   = nrow(relsug),       # number of NSC measurements
                   J_nsc   = as.numeric(factor(relsug$species)),  # species indicator for relative sugar data
                   nsc_obs = as.numeric(scale(relsug$relsug)))    # observed relative sugar values (scaled and centered)
)
# inspect names
names(data_fit)

# define initial values for latent variables and sensitive parameters to speed up convergence
# (all other parameters are initiated randomly using Stan standards)
init_latent <- list(
  phi      = rep(0.001, 12),                                 # misclassification rates
  nsc_true = (with(data_fit, tapply(nsc_obs, J_nsc, mean))), # true relative sugar values
  hsm_true = (data_fit$hsm_obs),                             # true safety margins
  sd_nsc   = .1,                                             # SD of within-species variability in relative sugar 
  c = .8                                                     # scaling parameter for neighbourhood effects
)


# 4. Fit full model in Stan ---------------------------------------------------
# the code for the model is in the "stan" subfolder of the repository
# Warning - model fitting takes 6-8 hours on a fast desktop computer
mod <- stan(
  file = "stan/neighbourhood_mortality_model.stan",
  data   = data_fit,
  iter   = 5000,
  chains = 4,
  # list of stored model parameters
  pars   = c("b_Intercept", "b_hsm", "b_nsc", "b", "gamma0",
             "sd_nsc", 
             "sd_1", "sd_2", "sd_3", "cor_2", "cor_neighbor",
             "nsc_true", "hsm_true", "phi", "c",
             "r_1_1", paste0("r_2_", 1:data_fit$M_2),  paste0("r_3_", 1:12),
             "log_lik"), 
  init   = function() init_latent,
  seed   = "1234567",
  control = list(max_treedepth = 15,
                 adapt_delta   = 0.95),
  verbose = TRUE,
  refresh = 10
)

# save model object
saveRDS(mod, file = "output/neighbourhood_mortality_model.rds")

# save exact image of the datasets used for modeling
save(list = c("data", "data_fit", "hsm1", "relsug"), 
     file = "output/neighbourhood_mortality_model_data.RData")

# 5. Model inspection via shinystan -------------------------------------------
# launch shinystan to inspect convergence characteristics, trace plots etc.
launch_shinystan(mod)

# 6. Approximate leave-one-out cross validation via the loo package -----------
# get log likelihood and effective sample size
log_lik <- as.array(mod, pars = "log_lik")
r_eff   <- relative_eff(log_lik)

# compute LOO CV
loo <- loo(x = log_lik, r_eff = r_eff, cores = 6) 
loo # no observations with Pareto k > 0.7, almost all < 0.5 

# 7. Export model estimates ---------------------------------------------------
# parameters for export
pars <-  c("b_Intercept", "b_hsm", "b_nsc", "b", "gamma0",
           "c", "sd_nsc", 
           "sd_1", "sd_2", "sd_3", "cor_neighbor", "cor_2",
           "phi", "nsc_true", "hsm_true", 
           paste0("r_2_", 1:data_fit$M_2))

# Get model summary
mcmcsumm <- summary(
  mod, 
  pars =pars
)$summary %>%
  as.data.frame() %>% 
  mutate(par = rownames(.))

# compute highest posterior density intervals
hdi <- posterior_samples(mod, pars = pars) %>% 
  gather(par, value) %>%
  group_by(par) %>% 
  summarize(hdi = list(mean_hdci(value))) %>% 
  unnest(hdi) %>% 
  ungroup %>% 
  select(par, hdi_min = ymin, hdi_max = ymax)

# join 
full_output <- full_join(hdi, mcmcsumm) %>% 
  select(par, mean:`97.5%`, hdi_min, hdi_max, n_eff, Rhat)
full_output

# get correct names for correlation parameters
names <- c("Int", "H", "Beet.neigh", "N", "NP", "P")

# (following the parameterization of the correlation matrix for the species-wise effects)
cors <- vector(mode = "character")
for (k in 1:6) {
  for (j in 1:(k - 1)) {
    cors[choose(k - 1, 2) + j] = paste0("cor(", names[j], "~", names[k], ")")
  }
}
cors

# replace correlation parameter names
full_output$par[grepl("cor_2", full_output$par)]
full_output$par[grepl("cor_2", full_output$par)] <- cors[c(1, 10:15, 2:9)]

# define new order for correlation parameters
c_order <- cors[c(1, 2, 4, 11, 7, 3, 5, 12, 8, 6, 13, 9, 14, 10, 15)]
c_order

# reorder and export model summary
mutate(
  full_output, 
  par = factor(par, 
               levels = c(mcmcsumm$par[1:32], 
                          c_order, 
                          mcmcsumm$par[48:155]),
               ordered = TRUE)) %>% 
  arrange(par) %>% 
  write_csv("output/full_model_output.csv")

