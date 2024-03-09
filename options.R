###########################################################
# OPTIONS
#
# Set key options for all things model related. The output
# of this function, o (a list), lives in the global environment,
# so can be referenced throughout the pipeline.
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
  
  # ---- Time settings ----
  
  # Years to analyse
  o$years = 1974 : 2024
  
  # Ages modelled
  o$ages = 0 : 95
  
  # Annual of 5-year bins for population data
  o$pop_bin = 1  # OPTIONS: 1 or 5
  
  # ---- Data and coverage settings ----
  
  # Bound coverage values below x%
  o$max_coverage = 0.99
  
  # Year that high-income countries switch to acellular pertussis vaccine
  o$wholecell_acellular_switch = 1995
  
  # Smooth static model coverage to better align with GBD burden estimates
  #
  # NOTE: Set to NULL to turn smoothing off
  o$gbd_coverage_smoother = "kernel"  # OPTIONS: "kernel" or "spline"
  
  # Standard deviation around maternal mean age at birth
  o$birth_age_sd = 3
  
  # ---- External models ----
  
  # Create/recreate dummy polio results
  o$dummy_polio = FALSE  # TODO: Remove when polio results are available
  
  # Directly simulate Dynamice model
  #
  # NOTE: If false, DynaMICE results must be otherwise available to user
  o$simulate_dynamice = FALSE
  
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
  o$selection_metric = "ll"
  
  # Multiply impact when fitting for more consistent FVP-impact scales
  o$impact_scaler = 1000

  # Defalt x scale for evaluating impact functions
  o$eval_x_scale = 3  # Not a critical value - often overwritten with actual FVPs
  
  # Number of initial years to average over for back-projecting impact ratio
  o$init_impact_years = 3
  
  # ---- Uncertainty settings ----

  # Flag for reproducible uncertainty draws (consistent randomly sampled seeds)
  o$uncertainty_reproducible = TRUE  # TODO: Implement this

  # Number of draws to sample
  o$n_draws = 200

  # Parameter bounds for fitting beta distribution to GBD disease vaccine efficacy
  o$par_lower = log(1)
  o$par_upper = log(10)

  # Statistical summary to use for 'best estimate' projection
  #
  # OPTIONS:
  #  "median" := Median of uncertainty simulations
  #    "mean" := Mean of uncertainty simulations
  # o$best_estimate_simulation = "mean"

  # Quantiles for credibility intervals
  # o$quantiles = c(0.025, 0.975)
  
  # ---- Parallelisation settings ----
  
  # Use multiple cores to speed up several fitting processes
  o$parallel = list(
    impute = FALSE, 
    impact = FALSE) 

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
  
  # Colour packages and palettes (see colour_scheme in auxiliary.R)
  o$palette = list(
    disease = "pals::kovesi.rainbow",  # 14 values needed
    region  = "brewer::paired",        # 6 values needed
    income  = "brewer::dark2")         # 4 values needed
  
  # WHO colour scheme (vibrant blue, navy blue)
  o$palette_who = c("#009CDE", "#001F58")

  # Font sizes: title, axis, tick, strip, legend, key
  o$font_size = c(34, 28, 16, 24, 20, 18)

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
  o$figure_format = "png" # Classic options: "png", "pdf", or "svg"

  return(o)
}

