Mutually inclusive mechanisms of drought-induced tree mortality
================

## Description

The present repository is part of Hajek, Link et al. (2021) and contains
the code necessary for preparing the dataset and fitting the
hierarchical Bayesian model described in the paper.

## Structure

The documents in the root folder (`/`) perform the basic steps of the
data analysis. The raw data are in the `/data` folder, the `/stan`
folder contains the Stan code required for model fitting, and
intermediate results and model output are stored in `/output`. Visual
output generated after running script 06 is stored in `/figures`.

    /           top-level scripts:
                  01_vulnerability_curves.R
                  02_hydraulic_safety_margins.R
                  03_impute_plot_margin_heights.R
                  04_compute_neighbourhood_matrices.R
                  05_stan_model_fitting.R
                  06_plot_model_output.R           
    /data       dataset used for model fitting
    /figures    visual output
    /output     saved model objects, model output etc.
    /stan       model code for the Stan model
