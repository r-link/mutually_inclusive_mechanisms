###############################################################################
#
#    Plotting model output
#
#    Script: Roman Link  (roman.link@plant-ecology.de)	
#
###############################################################################

# This script generates all figures in the publication based on the model output
# from Script 05.

# 0.1 Load packages -----------------------------------------------------------
# if pacman is not installed, install via install.packages("pacman").
# the package manager will take care of the installation of all subsequent 
# packages.
pacman::p_load(tidyverse, cowplot, brms, rstan, tidybayes, ggrepel)

# 0.2 Read model objects and data used for fitting ----------------------------
mod <- readRDS("output/neighbourhood_mortality_model.rds")
load("output/neighbourhood_mortality_model_data.RData")

# 0.3 Get posterior draws and split into meaningful units ---------------------
# species keys
spec <- sort(unique(data$species))

# response
Y <- as.numeric(data_fit$Y)

# posterior sample  (if your memory is limited, take a subsample, the
# calculations based on the posterior are very memory intensive)
post <- posterior_samples(mod) # %>% sample_n(1000)
dim(post)

# fixed effects parameter vector
beta <- select(post, b_Intercept, contains("b["), b_hsm, b_nsc) %>% 
  set_names(colnames(data_fit$X), "HSM", "NSC") 
dim(beta)

# main effects
main      <- select(beta, Intercept, HSM, NSC, height:fertP) %>% 
  set_names(c("Intercept", "HSM", "NSC", "Height", "Beetles" ,
              "Neigh. w/ beetles",  "N", "NP", "P")) %>% 
  select(everything(), "NP")

# latent variables from measurement models
nsc_true <- select(post, contains("nsc_true"))  %>% set_names(spec)
hsm_true <- select(post, contains("hsm_true"))  %>% set_names(spec)

nsc_true0 <- with(data_fit, tapply(nsc_obs, J_nsc, mean))
hsm_true0 <- data_fit$hsm_obs 

# species-wise HSM + NSC effects
hsm_eff  <- sweep(hsm_true, 1, beta$HSM, "*")
nsc_eff  <- sweep(nsc_true, 1, beta$NSC, "*")

# random plot effects
RE_plot <- select(post, contains("r_1_1"))
names(RE_plot) <-sort(unique(data$plotspec))

# species-wise random intercepts
RE_int <- select(post, contains("r_2_1"))  %>% set_names(spec)

# species-wise random height effects
RE_height <- select(post, contains("r_2_2")) %>% set_names(spec)

# species-wise random effects of beetle-effected neighbour trees
RE_beet <- select(post, contains("r_2_3"))  %>% set_names(spec)

# species-wise random N effects
RE_N <- select(post, contains("r_2_4"))  %>% set_names(spec)

# species-wise random NP effects
RE_NP <- select(post, contains("r_2_5"))  %>% set_names(spec)

# species-wise random N effects
RE_P <- select(post, contains("r_2_6"))  %>% set_names(spec)


# average neighbourhood effects
neigh_eff0 <- select(post, contains("gamma0"))
# neighbourhood effect scaling factor
c <- select(post, c) 

# species-wise neighbourhood effects
RE_neigh <- select(post, contains("r_3_"), -contains("r_3_0")) %>% 
  # get iteration indicator to enable reshaping to wide format
  mutate(iter = 1:nrow(.)) %>% 
  # make longtable
  gather(parm, est, -iter) %>% 
  # get indicators for species and neighbour species
  mutate(neighbour = str_match(parm, "r_3_(.*?)\\[")[,2],
         species  = str_match(parm, "\\[(.*?)\\]")[,2],
         neighbour = spec[as.numeric(neighbour)],
         species  = spec[as.numeric(species)]) %>% 
  # remove unnecessary variable
  select(-parm) %>% 
  # reshape to long form
  spread(neighbour, est) %>% 
  arrange(species, iter)
RE_neigh

# get full neighbourhood effects
# ...start with random component
neigh_eff <- RE_neigh

# ...add main effect
for (i in 1:12){
  neigh_eff[,i + 2] <- neigh_eff[,i + 2] + rep(neigh_eff0[, i], 12)
}

# get mean neighbourhood effects matrix
neigh_eff_mean <- neigh_eff %>% 
  gather(neighbour, est, ACPL:QURU) %>% 
  group_by(species, neighbour) %>% 
  summarize(mhdi = list(mean_hdi(est)[1:3])) %>% 
  unnest(cols = mhdi) %>% 
  ungroup() %>% 
  mutate(cred = ymin > 0 | ymax < 0)
neigh_eff_mean


# get one inflation parameters
phi_raw <- posterior_samples(mod, "phi") %>% as.matrix
# summarize one inflation parameters
phi <- as.data.frame(phi_raw) %>% 
  set_names(spec) %>% 
  gather(species, value) %>%
  group_by(species) %>% 
  summarize(mhdi = list(mean_hdci(value)[1, 1:3])) %>% 
  unnest(mhdi)  %>% 
  mutate(
    dir = case_when(ymax < 0 ~ "neg",
                    ymin > 0 ~ "pos",
                    TRUE ~ "0"),
    dir = factor(dir, levels = c("neg", 0, "pos"), ordered = TRUE),
    species = factor(species, levels = ord, ordered = TRUE),
    group = "Misclassification rates")

# 0.4 Predicted average mortality rates ----------------------------------------
# get full conditional predictions from log likelihood
# -- Explanation: the log_lik has to be exported to be able to analyze the model
# with approximate leave-one-out crossvalidation via loo().
# As the data are simulated as a Bernoulli process, it is easy to compute the
# likelihood for each observation & to marginalize out the misclassification rate
# after the fact to reduce the number of exported parameters.
lik <- exp(extract_log_lik(mod)) 

# get corresponding one inflation parameters for each observation
phimat <- phi_raw[,data_fit$J_2]

# get true mortality rate probability
Yhat0 <- matrix(nrow = nrow(lik), ncol = ncol(lik)) 
Yhat0[, Y == 0]  <- 1 - lik[, Y == 0] / (1 - phimat[, Y == 0])
Yhat0[, Y == 1]  <- (lik[, Y == 1] -  phimat[, Y == 1]) / (1 - phimat[, Y == 1])

# get probability of observing dead tree
Yhat1 <- matrix(nrow = nrow(lik), ncol = ncol(lik)) 
Yhat1[, Y == 0]  <- 1 - lik[, Y == 0]
Yhat1[, Y == 1]  <- lik[, Y == 1]

# compute species average predictions for figure 1
# WARNING: this step takes a lot of time and memory!
specmort_prob <- t(Yhat0) %>%
  data.frame() %>% 
  mutate(species = data$species) %>% 
  group_by(species) %>% 
  summarize_all(mean) %>% 
  gather(key, value, -species) %>% 
  group_by(species) %>% 
  summarize(hdi = list(mean_hdci(value)[1:3])) %>% 
  unnest(hdi)

# get order for species indicators
ord <- arrange(specmort_prob, y)$species

# 0.5 Pseudo-Rsquared ---------------------------------------------------------
# function for explained variance following Gelman et al. (2018)
rsq_bayes <- function(pred, Y){
  res <- -sweep(pred, 2, Y)
  var_y   <- apply(pred, 1, var)
  var_res <- apply(res, 1, var)
  var_y / (var_y + var_res)
}

# compute explained variance
(Rsq_cond <- mean_hdci(rsq_bayes(Yhat1, Y)))  # after removing misclassification .514 
(Rsq_cond1 <- mean_hdci(rsq_bayes(Yhat0, Y))) # with misclassification .515 - practically identical

 # baseline: variance explained by species-averages alone
data %>%  
  group_by(species) %>% 
  mutate(mmort = mean(mort1, na.rm = TRUE),
         res   = mort1 - mmort) %>% 
  ungroup %>% 
  with(var(mmort) / (var(mmort) + var(res)))


# Fig. 1. Overall species predictions -----------------------------------------
# get data for the overall average mortality probability
avgmort <- mean_hdci(rowMeans(Yhat0)) %>%  mutate(xmin = -Inf, xmax = Inf)

# get observed average mortality rates
obsmort <- data %>% 
  group_by(species) %>% 
  summarize(m = mean(mort1)) %>% 
  mutate(species = factor(species, levels = ord, ordered = TRUE)) 

# combine in one single plot
mort_plot <- specmort_prob  %>% 
  mutate(dir = case_when(ymax < .3 ~ "neg",
                         ymin > .39 ~ "pos",
                         TRUE ~ "0"),
         dir = factor(dir, levels = c("neg", 0, "pos"), ordered = TRUE),
         group = "Overall species-specific mortality effects",
         species = factor(species, levels = ord, ordered = TRUE)) %>% 
  ggplot(aes(x = species)) +
  geom_point(data = obsmort, aes(y = m), pch = 1, size = 1.8) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
            data = avgmort, fill = "grey", alpha = .6, inherit.aes = FALSE) +
  geom_hline(yintercept = avgmort$y, lty = 2, size = .25) +
  geom_pointrange(aes(y = y , ymin = ymin, ymax = ymax, col = dir), fatten = -.5, size = .25) +
  geom_text(aes(x, y, label = label), data = tibble(x =  ord[10], y = avgmort$y - .07, label = "Overall avg. mortality rate"),
            inherit.aes = FALSE, size = 2.8)  +
  scale_color_manual(values = c("blue", "grey50", "red")) +
  facet_wrap(~group) +
  theme_linedraw(base_size = 9) +
  labs(y = "Estimated mortality rate") +
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1), units = "mm"),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(angle = 90, size = 6, vjust = .5)) 
mort_plot

# Fig. 2. Observed vs. predicted ----------------------------------------------
# compute species-wise plot average posterior predictions
ovp <-tibble(
  pred = colMeans(Yhat1), 
  species = data$species, 
  plot = data$plot_cn, 
  block = data$block, 
  mort = data$mort) %>% 
  group_by(species, plot) %>% 
  summarize_if(is.numeric, mean) 

# plot observed versus predicted species-wise plot averages 
ovp_plot <- ovp %>%
  mutate(group = "Observed and predicted plot-level mortality") %>% 
  ggplot(aes(x = mort, y = pred)) + 
  geom_point(alpha = 0.3, col = 4, size = 1) +
  #geom_smooth(method = "lm", col = 1, lty = 2, se = F) +
  geom_abline()  +
  facet_wrap(~group) +
  theme_linedraw(base_size = 9) +
  labs(x = "Observed mortality rate", y = "Predicted mortality rate") +
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1), units = "mm")) 
ovp_plot

# plot level (raw) rsq
with(ovp, var(pred) / (var(pred) + var(res1))) 
#  > 0.9 - high proportion of variance explained at plot level

# Fig. 3. HSM and NSC effects -------------------------------------------------
## .....a. HSM -----
hsm

# function for retransformation of scaled and centred HSM
retrans_hsm <- function(x) x * sd(hsm$HSM) + mean(hsm$HSM)

# compute values for smooth curve
hsms <- seq(min(hsm$HSM) - .1, max(hsm$HSM) + .1, .01)
hsms_trans <- ((hsms) - mean(hsm$HSM)) / sd(hsm$HSM)

# compute smooth predictions
hsm_pred <- hsms_trans %*% t(beta$HSM) %>% 
  sweep(2, beta$Intercept, "+") %>% 
  as.data.frame %>% 
  mutate(hsm = hsms, hsm_trans = hsms_trans) %>% 
  gather(var, value, -hsm, -hsm_trans) %>% 
  group_by(hsm, hsm_trans) %>% 
  summarize(mhdi = list(mean_hdci(value)[1, 1:3])) %>% 
  unnest(mhdi)
hsm_pred

# get true values
hsm_true_summary <- hsm_true %>% 
  gather(species, value) %>% 
  mutate(value = retrans_hsm(value)) %>% 
  group_by(species) %>% 
  summarize(hdi = list(mean_hdci(value))) %>% 
  unnest(hdi) %>% 
  set_names(gsub("y", "hsm", names(.)))
hsm_true_summary

# plot output
hsm_plot <- sweep(hsm_eff + RE_int,1, beta$Intercept, "+" ) %>% 
  as.data.frame %>% 
  gather(species, value) %>% 
  group_by(species) %>% 
  summarize(mhdi = list(mean_hdci(value)[1, 1:3])) %>% 
  unnest(mhdi) %>% 
  mutate(species = factor(species, levels = ord, ordered = TRUE),
         hsm0 = hsm$HSM,
         main = "HSM effect"
  ) %>% 
  full_join(hsm_true_summary) %>% 
  ggplot(aes(x = hsm)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), data = hsm_pred, fill = "grey") +
  geom_line(aes( y = y), data = hsm_pred) +
  geom_linerange(aes(xmin = hsmmin, xmax = hsmmax, y = y)) +
  geom_pointrange(aes(y = y, ymin = ymin, ymax = ymax)) +
  geom_point(aes(x = hsm0, y = y), col = 4) +
  geom_label_repel(aes(x = hsm, y = y,  label = species), size = 2.5) +
  theme_linedraw() +
  facet_wrap(~main) +
  labs(x = "Hydraulic safety margin (MPa)", 
       y = "Partial effect on mortality (logit scale)")
hsm_plot

## .....b. NSC -----
# function for retransformation of the standardized relative sugar content changes
retrans_nsc <- function(x) (x * sd(relsug$relsug) + mean(relsug$relsug))

# average relative sugar content changes
rs <- with(relsug, tapply(relsug, species, mean))

# compute values for smooth curve
nsc <- seq(min(rs) - .04, max(rs) + 0.04, .01)
nsc_trans <- ((nsc) - mean(relsug$relsug)) / sd(relsug$relsug)

# compute smooth predictions
nsc_pred <- nsc_trans %*% t(beta$NSC) %>% 
  sweep(2, beta$Intercept, "+") %>% 
  as.data.frame %>% 
  mutate(nsc = nsc * 100, nsc_trans = nsc_trans) %>% 
  gather(var, value, -nsc, -nsc_trans) %>% 
  group_by(nsc, nsc_trans) %>% 
  summarize(mhdi = list(mean_hdci(value)[1, 1:3])) %>% 
  unnest(mhdi)
nsc_pred

# get true values
nsc_true_summary <- nsc_true %>% 
  gather(species, value) %>% 
  mutate(value = retrans_nsc(value) * 100) %>% 
  group_by(species) %>% 
  summarize(hdi = list(mean_hdci(value))) %>% 
  unnest(hdi) %>% 
  set_names(gsub("y", "nsc", names(.)))
nsc_true_summary

# prepare data for plot
nsc_plot <- sweep(nsc_eff + RE_int, 1, beta$Intercept, "+" ) %>% 
  as.data.frame %>% 
  gather(species, value) %>% 
  group_by(species) %>% 
  summarize(mhdi = list(mean_hdci(value)[1, 1:3])) %>% 
  unnest(mhdi) %>% 
  mutate(species = factor(species, levels = ord, ordered = TRUE),
         nsc0 = with(relsug, tapply(relsug, species, mean)) * 100,
         main = "NSC effect") %>% 
  full_join(nsc_true_summary) %>% 
  ggplot(aes(x = nsc)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), data = nsc_pred, fill = "grey") +
  geom_line(aes( y = y), data = nsc_pred) +
  geom_linerange(aes(xmin = nscmin, xmax = nscmax, y = y)) +
  geom_pointrange(aes(y = y, ymin = ymin, ymax = ymax)) +
  geom_point(aes(x = nsc0, y = y), col = 4) +
  geom_label_repel(aes(x = nsc, y = y,  label = species), size = 2.5) +
  theme_linedraw() +
  facet_wrap(~main) +
  labs(x = "Change in relative leaf sugar content (%)", 
       y = "Partial effect on mortality (logit scale)")
nsc_plot

# ..... combine plots -----
hsmnsc <- plot_grid(hsm_plot, nsc_plot, 
                    labels = c("a)", "b)"), 
                    label_size = 10, 
                    nrow = 1)
hsmnsc

# Fig. 4. Main effects --------------------------------------------------------
# compute and plot main effects
main_effs <- main %>% 
  gather() %>%
  mutate(key = factor(key, levels = unique(key), ordered = TRUE)) %>% 
  group_by(key) %>% 
  summarize(mhdi = list(mean_hdci(value)[1, 1:3])) %>% 
  unnest(mhdi)  %>% 
  mutate(dir = case_when(ymax < 0 ~ "neg",
                         ymin > 0 ~ "pos",
                         TRUE ~ "0"),
         dir = factor(dir, levels = c("neg", 0, "pos"), ordered = TRUE),
         key   = fct_relevel(key, "NP", after = Inf),
         group = "Main effects") %>% 
  ggplot(aes(x = key, y , ymin = ymin, ymax = ymax, col = dir)) +
  geom_pointrange(fatten = 0, size = .3) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_color_manual(values = c("blue", "grey50", "red")) +
  facet_wrap(~group) +
  theme_linedraw(base_size = 9) +
  labs(y = "Estimate (logit scale)") +
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1), units = "mm"),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(angle = 90, size = 6, vjust = .5, hjust =1)) 

main_effs

# Fig. 5. Neighbourhood effect heatmap ----------------------------------------
# .....a. average neighbourhood effects ------
neigh_effs_av <- neigh_eff0 %>% 
  set_names(spec) %>% 
  gather() %>%
  mutate(key = factor(key, levels = unique(key), ordered = TRUE)) %>% 
  group_by(key) %>% 
  summarize(mhdi = list(mean_hdci(value)[1, 1:3])) %>% 
  unnest(mhdi)  %>% 
  mutate(dir = case_when(ymax < 0 ~ "neg",
                         ymin > 0 ~ "pos",
                         TRUE ~ "0"),
         dir = factor(dir, levels = c("neg", 0, "pos"), ordered = TRUE),
         group = "Average neighbourhood effects",
         key = factor(key, levels = ord, ordered = TRUE)) %>% 
  ggplot(aes(x = key, y , ymin = ymin, ymax = ymax, col = dir)) +
  geom_pointrange(fatten = 0, size = .3) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_color_manual(values = c("blue", "grey50", "red")) +
  facet_wrap(~group) +
  theme_linedraw(base_size = 9) +
  labs(y = "Est. (logit scale)") +
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1), units = "mm"),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(angle = 90, size = 6, vjust = .5)) 

neigh_effs_av

# .....b. full matrix of interactions ------
neigh_mat <- neigh_eff_mean %>%
  mutate(species  = factor(species, levels = rev(ord), ordered = TRUE),
         neighbour = factor(neighbour, levels = ord, ordered = TRUE)) %>% 
  mutate(fac = "Neighbourhood interaction matrix") %>% 
  ggplot(aes(x = neighbour, y = species)) +
  geom_tile(aes(fill = y), show.legend = TRUE, col = 1) + 
  geom_tile(data = function(x) filter(x, species == neighbour), 
            fill = NA, col = 1, size = .5) + 
  geom_point(aes(shape = cred), col = "lightgrey", size = 2.5) +
  scale_fill_gradient2(low = "blue", mid = "grey95", high = "red", 
                       guide = "colourbar") +
  scale_shape_manual(values = c("", "\U2605")) +
  theme_linedraw(base_size = 9) +
  theme(legend.position = "bottom",
        aspect.ratio = 1,
        plot.margin = unit(c(1,1,1,1), units = "mm"),
        axis.text   = element_text(size = 6),
        axis.text.x  = element_text(angle = 90, vjust = .5)) +
  facet_wrap(~fac) +
  coord_cartesian(expand = FALSE) +
  labs(x = "Species identity of neighbour tree", y = "Species identity of focal tree", fill = "Neighbour effect") +
  guides(shape = "none")
neigh_mat

# ..... combine plots -----
neigh_plot <- plot_grid(
  neigh_effs_av, neigh_mat, 
  nrow = 2, align = "hv", 
  labels = c("a)", "b)"),
  label_size = 10,
  axis = "tlr",
  rel_heights = c(0.3, 1))
neigh_plot                         


# Fig. S2.1. Misclassification parameter --------------------------------------
ones <- phi %>% 
  ggplot(aes(x = species, y , ymin = ymin, ymax = ymax, col = dir)) +
  geom_pointrange(fatten = 0, size = .3) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_color_manual(values = c("blue", "grey50", "red")) +
  facet_wrap(~group) +
  theme_linedraw(base_size = 9) +
  labs(y = "Estimate") +
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1), units = "mm"),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(angle = 90, size = 6, vjust = .5, hjust =1)) 
ones

# Fig. S2.2. Species level effects --------------------------------------------
# prepare dataset of species-wide effects (main effect + random species effect)
specwise_dat <- bind_rows(
  mutate(RE_height + beta$height, group = "Height effect"),
  mutate(RE_beet + beta$beetles_neigh, group = "Beetle neig. eff."),
  mutate(RE_int, group = "Unexplained species diffs."),
  mutate(RE_N + beta$fertN, group = "N effect"),
  mutate(RE_P + beta$fertP, group = "P effect"),
  mutate(RE_NP + beta$fertNP, group = "NP effect")
) %>% 
  mutate(group = factor(group, levels = unique(group), ordered = TRUE)) %>% 
  gather(species, est, -group) %>% 
  group_by(group, species) %>%   
  summarize(mhdi = list(mean_hdci(est)[1, 1:3])) %>% 
  unnest(mhdi)  %>% 
  ungroup

specwise_dat1 <- specwise_dat %>%  
  mutate(group = factor(group, levels = unique(group), ordered = TRUE)) %>% 
  mutate(dir = case_when(ymax < 0 ~ "neg",
                         ymin > 0 ~ "pos",
                         TRUE ~ "0"),
         dir = factor(dir, levels = c("neg", 0, "pos"), ordered = TRUE),
         species = factor(species, levels = ord, ordered = TRUE)
  ) %>% 
  # remove fertilizer effects of unfertilized species (information-free samples from prior)
  filter(!(group %in% c("N effect", "P effect", "NP effect") & 
             species %in% c("BEPA", "LALA", "PIST", "PIGL", "QURU", "ACSA")))

# plot species-wise effects
specwise_effs <- specwise_dat1 %>% 
  ggplot(aes(x = y, xmin = ymin, xmax = ymax, y = species, col = dir)) +
  geom_pointrange(fatten = 0, size = .3) +
  geom_vline(xintercept = 0, lty = 2) +
  facet_wrap(~group, scales = "free_x", nrow = 2) +
  scale_color_manual(values = c("blue", "grey50", "red")) +
  theme_linedraw(base_size = 9) +
  labs(x = "Estimate (logit scale)") +
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1), units = "mm"),
        axis.title.y = element_blank()
  )
specwise_effs

# Fig. S2.3 Variance parameters for appendix ----------------------------------
# get correct variable names
names <- c("Int", "H", "Beet.neigh", "N", "NP", "P")
# (following the parameterization of the correlation matrix for the species-wise effects)
newnames <- vector(mode = "character")
for (k in 1:6) {
  for (j in 1:(k - 1)) {
    newnames[choose(k - 1, 2) + j] = paste0("rho[beta]*('", names[j], ",", names[k], "')")
  }
}
newnames
# add the other names
newnames1 <- c(newnames,
               "rho[neigh]~('Neighbour pair corr.')", 
               "tau[delta]~('SD of spec.wise plot effs.')", 
               paste0("tau[beta[", 1:6,"]]~('SD for ", names, "')"),
               "tau[gamma]~('SD of neigh. effects')", 
               "sigma['Δsugar']~('Δsugar meas. error')")

# extract parameters and reshape
var_pars <- select(post, contains("cor"), contains("sd")) %>% 
  set_names(newnames1) %>% 
  gather() %>% 
  group_by(key) %>% 
  summarize(hdi = list(mean_hdci(value))) %>% 
  unnest(hdi) %>% 
  mutate(group = ifelse(!grepl("rho", key), "Standard deviations", "Correlation parameters"),
         group = factor(group, levels = rev(unique(group)), ordered = TRUE),
         key   = factor(key, levels = rev(newnames1[c(1, 2, 4, 11, 7, 3, 5, 12, 8, 6, 13, 9, 14, 10, 15,
                                                      16, 18:24, 17, 25)]),
                        ordered = TRUE)
  )

# create plot
varplot <- var_pars %>% 
  ggplot(aes(y, y = key, xmin = ymin, xmax = ymax)) +
  geom_vline(lty = 2, col = "grey", xintercept = 0) +
  geom_pointrange(fatten = 0, size = .3) +
  geom_hline(yintercept = 0, lty = 2) +
  facet_wrap(~group, scales = "free") +
  theme_linedraw(base_size = 9) +
  labs(x = "Estimate (logit scale)") +
  scale_y_discrete(labels = function(x) parse(text = x)) +
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1), units = "mm"),
        axis.title.y = element_blank()) 

varplot


# Fig. S2.4. Direction of neighbourhood effects -------------------------------
head(neigh_eff)
# calculate average neighbourhood effects
test <- neigh_eff %>% 
  gather(neighbour, value, -species, -iter)
test 

# get dataset with reverse effects
test1 <- rename(test, species = neighbour, neighbour = species, 
                value_rev = value)
test1

# merge datasets
test2 <- full_join(test, test1)  %>% 
  group_by(species, neighbour) %>% 
  summarize(dir1 = list(mean_hdci(value)[, 1:3]), 
            dir2 = list(mean_hdci(value_rev)[, 1:3]),
            diff = list(mean_hdci(abs(value - value_rev))[, 1:3]),
            sum  = list(mean_hdci((value + value_rev))[, 1:3]))   %>% 
  mutate(dir2 = map(dir2, ~set_names(.x, paste0(names(.x), "_rev"))),
         diff = map(diff, ~set_names(.x, paste0(names(.x), "_diff"))),
         sum  = map(sum,  ~set_names(.x, paste0(names(.x), "_sum")))) %>% 
  unnest(cols = c(dir1, dir2, diff, sum)) %>% 
  mutate(pair = apply(cbind(species, neighbour), 1, function(x) paste(sort(x), collapse = "-"))) %>% 
  ungroup %>% 
  filter(!duplicated(pair))
test2

# plot of the effects for the pairings with the strongest disparity
disp <- arrange(test2, -y_diff)[1:5, ] %>% 
  mutate(group = "Pairs w/ strongest disparities",
         pair = factor(pair, levels = unique(pair), ordered = TRUE)
  ) %>% 
  ggplot(aes(y = pair)) +
  geom_pointrange(aes(x = y, xmin = ymin, xmax = ymax, col = "Original"), fatten = 0) +
  geom_pointrange(aes(x = y_rev, xmin = ymin_rev, xmax = ymax_rev, col = "Reverse"), fatten = 0) +
  geom_vline(aes(xintercept = 0), lty = 2) +
  facet_wrap(~group) +
  scale_color_viridis_d(begin = .07, end = .91) +
  labs(x = "Estimate (logit scale)") +
  theme_linedraw(base_size = 9) + 
  theme(legend.position = "none", 
        axis.title.y = element_blank())

# plot of the effects for the pairings with the strongest beneficial effects
benef <- arrange(test2, y_sum)[1:5, ] %>% 
  mutate(group = "Most beneficial pairings",
         pair = factor(pair, levels = unique(pair), ordered = TRUE)
  ) %>% 
  ggplot(aes(y = pair)) +
  geom_pointrange(aes(x = y, xmin = ymin, xmax = ymax, col = "Original"), fatten = 0) +
  geom_pointrange(aes(x = y_rev, xmin = ymin_rev, xmax = ymax_rev, col = "Reverse"), fatten = 0) +
  geom_vline(aes(xintercept = 0), lty = 2) +
  facet_wrap(~group) +
  labs(x = "Estimate (logit scale)") +
  scale_color_viridis_d(begin = .07, end = .91) +
  theme_linedraw(base_size = 9) + 
  theme(legend.position = "none", 
        axis.title.y = element_blank())

# plot of the effects for the pairings with the most detrimental effects
detri <- arrange(test2, -y_sum)[1:5, ] %>% 
  mutate(group = "Most detrimental pairings",
         pair = factor(pair, levels = unique(pair), ordered = TRUE)
  ) %>% 
  ggplot(aes(y = pair)) +
  geom_pointrange(aes(x = y, xmin = ymin, xmax = ymax, col = "Original"), fatten = 0) +
  geom_pointrange(aes(x = y_rev, xmin = ymin_rev, xmax = ymax_rev, col = "Reverse"), fatten = 0) +
  geom_vline(aes(xintercept = 0), lty = 2) +
  facet_wrap(~group) +
  labs(x = "Estimate (logit scale)") +
  scale_color_viridis_d( name = "Direction",begin = .07, end = .91) +
  theme_linedraw(base_size = 9) + 
  theme(axis.title.y = element_blank())

# combine the three subplots
neigh_pairs <- plot_grid(disp, benef, detri + theme(legend.position = "none"),
                         get_legend(detri), nrow = 1, rel_widths = c(1,1,1,.4))

neigh_pairs


# Fig. S2.5. Plot effects -----------------------------------------------------
# plot effects sorted by posterior mean
plot_effs <- RE_plot %>% 
  gather() %>%
  mutate(key = factor(key, levels = unique(key), ordered = TRUE)) %>% 
  group_by(key) %>% 
  summarize(mhdi = list(mean_hdci(value)[1, 1:3])) %>% 
  unnest(mhdi)  %>% 
  separate(key, into = c("plot", "species"), sep = "_", remove = FALSE) %>% 
  mutate(dir = case_when(ymax < 0 ~ "neg",
                         ymin > 0 ~ "pos",
                         TRUE ~ "0"),
         key = fct_reorder(key, y),
         dir = factor(dir, levels = c("neg", 0, "pos"), ordered = TRUE),
         group = "Varying species-wise plot effects") %>% 
  ggplot(aes(x = key, y =  y , ymin = ymin, ymax = ymax, col = dir)) +
  geom_pointrange(fatten = 0, size = .3, alpha = 0.7) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_color_manual(values = c("blue", "grey50", "red")) +
  facet_wrap(~group) +
  theme_linedraw(base_size = 9) +
  labs(y = "Estimate (logit scale)", x = "Plot/species combination") +
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1), units = "mm"),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )
plot_effs 
# plot identity is not supposed to be visible - this plot is about the distribution
# of the plot effects

# 1.0. Export output ----------------------------------------------------------
# Main text
ggsave("figures/fig1_overall_mortality.png", mort_plot, width = 11.25, 
       height = 5, units = "cm")
ggsave("figures/fig2_observed_predicted.png", ovp_plot, width = 8.5, 
       height = 6, units = "cm")
ggsave("figures/fig3_hsm_nsc.png", hsmnsc, width = 18, 
       height = 8, units = "cm")
ggsave("figures/fig4_main_effects.png", main_effs, width = 11.25, 
       height = 5.5, units = "cm")
ggsave("figures/fig5_neighbourhood_effects.png", neigh_plot, width = 11.25, 
       height = 18, units = "cm")


# Supplementary material
ggsave("figures/S2-1_one_inflation.png", ones, width = 11.25, 
       height = 5, units = "cm")
ggsave("figures/S2-2_species_level_effects.png", specwise_effs, width = 16, 
       height = 12, units = "cm")
ggsave("figures/S2-3_variance_parameters.png", varplot, width = 14.5, 
       height = 16, units = "cm")
ggsave("figures/S2-4_neighbour_pairs.png", neigh_pairs, width = 18, 
       height = 8, units = "cm")
ggsave("figures/S2-5_plot_effects.png", plot_effs, width = 18, 
       height = 5, units = "cm")
