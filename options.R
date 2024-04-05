###########################################################
# OPTIONS
#
# Set key options for all things model related. The output
# of this function, o (a list), lives in the global 
# environment, so can be referenced throughout the pipeline.
#
###########################################################

# ---------------------------------------------------------
# Set model options and assumptions
# ---------------------------------------------------------
set_options = function(do_step = NA) {

  message("* Setting options")

  # Several global R settings to make life easier
  default_R_options()  # See auxiliary.R

  # Initiate options list
  o = list(do_step = do_step)

  # Prepare output directory system
  o = prepare_dirs(o)  # See directories.R
  
  # ---- Analysis settings ----
  
  # Disease burden metrics to report
  o$metrics = c("deaths", "dalys")
  
  # ---- Time settings ----
  
  # Years to analyse
  o$years = 1974 : 2024
  
  # Ages modelled
  o$ages = 0 : 100
  
  # Annual of 5-year bins for population data
  o$pop_bin = 1  # OPTIONS: 1 or 5
  
  # Model final year only partially
  #
  # NOTE: This represents results up to May 2024
  o$partial_final_year = 5 / 12
  
  # ---- Data and coverage settings ----
  
  # Bound coverage values below x%
  o$max_coverage = 0.995
  
  # Year that high-income countries switch to acellular pertussis vaccine
  o$wholecell_acellular_switch = 1995
  
  # Smooth static model coverage to better align with GBD burden estimates
  #
  # NOTE: Set to NULL to turn smoothing off
  o$gbd_coverage_smoother = "kernel"  # OPTIONS: "kernel" or "spline"
  
  # Define method for extrapolating GBD burden post 2019
  o$gbd_extrap = "constant"  # OPTIONS: "trend" or "constant"
  
  # Flag to force (re)download Gapminder data from github
  o$force_download_gapminder = FALSE
  
  # ---- External models ----
  
  # Directly simulate Dynamice model
  #
  # NOTE: If false, DynaMICE results must be otherwise available to user
  o$simulate_dynamice = TRUE
  
  # GitHub repo for simulating DynaMICE model for EPI50 analysis
  o$github_dynamice = "ashattock/dynamice"
  
  # ---- Global Burden of Disease settings ----
  
  # Use GBD estimates starting from a given year
  o$gbd_estimate_years = 1990 : 2019
  
  # ---- Regression settings ----
  
  # Country must have at least n data points for fitting
  o$min_data_requirement = 5
  
  # Length of each 'period' in years
  o$period_length = 10

  # ---- Impact function settings ----
  
  # Impact function model selection metric
  #
  # OPTIONS:
  #  aicc := Akaike information criterion score
  #    ll := Log likelihood score
  o$selection_metric = "aicc"
  
  # Multiply impact when fitting for more consistent FVP-impact scales
  o$impact_scaler = 1000
  
  # Number of times to repeat optimisation
  o$n_optim = 10
  
  # Prior characteristics (assuming Gaussian)
  o$prior_weight = 1
  o$prior_sd     = 0.05
  
  # MCMC run time settings
  o$mcmc_burnin  = 100
  o$mcmc_iter    = 1000
  
  # MCMC posterior samples (ensure greater than required uncertainty_samples)
  o$mcmc_samples = 100

  # Defalt x scale for evaluating impact functions
  o$eval_x_scale = 5  # Not a critical value - often overwritten with actual FVPs
  
  # Number of initial years to average over for back-projecting impact ratio
  o$init_impact_years = 3
  
  # ---- Uncertainty settings ----

  # Number of posterior samples for generating uncertainty
  # 
  # NOTE: Bounded by value of mcmc_samples
  o$uncertainty_samples = 10
  
  # Quantiles for prediction intervals
  o$quantiles = c(0.025, 0.975)  # Represents 95% bounds
  
  # ---- Parallelisation settings ----
  
  # Use multiple cores to speed up several processes
  o$parallel = list(
    interp  = TRUE,  # NOTE: Occurs in two places in preparation step
    impute  = TRUE, 
    impact  = FALSE,  # NOTE: Having issues with shared memory
    history = TRUE)

  # Detect number of cores available to this user
  o$n_cores = detectCores()

  # ---- Plotting flags ----

  # Turn figures on or off
  o$plot_inputs     = FALSE
  o$plot_static     = FALSE
  o$plot_imputation = FALSE
  o$plot_impact     = FALSE
  o$plot_history    = TRUE

  # ---- Plotting settings ----

  # Saved figure size
  o$save_width  = 14
  o$save_height = 10

  # Units of figures sizes
  o$save_units = "in"

  # Plotting resolution (in dpi)
  o$save_resolution = 300

  # Image format for saving figure
  #
  # NOTE: Use a character vector to save with multiple formats at once
  o$figure_format = c("png", "svg")  # Lancet requires: "eps" or "svg"

  return(o)
}

