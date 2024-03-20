###########################################################
# IMPACT
#
# An alternative to impact factor for representing impact as
# a function of vaccine coverage.
#
###########################################################

# ---------------------------------------------------------
# Parent function for calculating non-linear impact
# ---------------------------------------------------------
run_impact = function() {
  
  # Only continue if specified by do_step
  if (!is.element(5, o$do_step)) return()
  
  message("* Fitting impact functions")
  
  # ---- FVPs and impact estimates ----

  message(" - Preparing FVP-impact data")

  # Prepare impact-FVP data to fit to
  data_dt = get_impact_data()

  # Exploratory plots of data used to fit impact functions
  plot_impact_data()

  # Plot all-time impact per FVPs
  plot_impact_fvps(scope = "all_time")

  # ---- Model fitting ----

  message(" - Evaluating impact functions")

  # Display number of cores if running in parallel
  if (o$parallel$impact)
    message("  > Using ", o$n_cores, " cores")
  
  # Country-disease-vaccine-activity combinations
  run_dt = data_dt %>%
    select(d_v_a_id, country) %>%
    unique() %>%
    mutate(d_v_a_str = str_pad(d_v_a_id, 2, pad = "0"),
           run_id    = paste1(country, d_v_a_str)) %>%
    select(-d_v_a_str)

  # Initiate progress bar
  pb = start_progress_bar(nrow(run_dt))

  # Run get_best_model in parallel
  if (o$parallel$impact)
    mclapply(
      X    = run_dt$run_id,
      FUN  = get_best_model,
      run  = run_dt,
      data = data_dt,
      mc.cores = n_cores,
      mc.preschedule = FALSE)

  # Run get_best_model consecutively
  if (!o$parallel$impact)
    lapply(
      X    = run_dt$run_id,
      FUN  = get_best_model,
      run  = run_dt,
      data = data_dt,
      pb   = pb)
  
  # Close progress bar
  close(pb)

  # ---- Model selection ----

  # Compile all results
  compile_results(run_dt)

  # Select best function for each country-d_v_a combination
  model_selection()

  # ---- Plot results ----

  # Plot function selection statistics
  plot_model_selection()
  
  # Plot impact vs coverage by vaccine, income, and decade 
  # plot_impact_coverage()
  
  # Plot impact function evaluation
  plot_model_fits()
}

# ---------------------------------------------------------
# Prepare impact-FVP data to fit to
# ---------------------------------------------------------
get_impact_data = function() {
  
  # Population size of each country over time
  pop_dt = table("wpp_pop") %>%
    lazy_dt() %>%
    group_by(country, year) %>%
    summarise(pop = sum(pop)) %>%
    ungroup() %>%
    as.data.table()
  
  # Load impact estimates from VIMC (inlcuding imputed)
  #
  # NOTE: Result of imputation is already in cumulative form
  vimc_dt = read_rds("impute", "impute_result", err = FALSE)
  
  # Load static model impact estimates
  static_dt = 
    read_rds("static", "deaths_averted_vaccine", err = FALSE) %>%
    lazy_dt() %>%
    # Scale results to per capita...
    left_join(y  = pop_dt, 
              by = c("country", "year")) %>%
    mutate(fvps   = fvps / pop, 
           impact = impact / pop) %>%
    select(-pop) %>%
    # Convert to cumulative FVP and impact...
    group_by(country, d_v_a_id) %>%
    mutate(fvps   = cumsum(fvps), 
           impact = cumsum(impact)) %>%
    ungroup() %>%
    as.data.table()
  
  # Impact estimates per capita from all sources
  data_dt = rbind(vimc_dt, static_dt) %>%
    mutate(impact_fvp = impact / fvps)
  
  # Save to file
  save_rds(data_dt, "impact", "data") 
  
  return(data_dt)
}

# ---------------------------------------------------------
# Set of functions to fit - we'll determine the 'best'
# ---------------------------------------------------------
fn_set = function(dict = FALSE) {
  
  # Set of statistical models / functions we want to test
  out = list(
    lin = function(x, a)       y = a * x,
    log = function(x, a, b)    y = logarithmic_growth(x, a, b),
    exp = function(x, a, b)    y = exponential_growth(x, a, b),
    sig = function(x, a, b, c) y = sigmoidal_growth(x, a, b, c))
  
  # Alternative functionality - return dictionary
  if (dict == TRUE)
    out = c(
      lin = "Linear gradiant (1 parameter)", 
      log = "Logarithmic growth (2 parameters)",
      exp = "Exponential growth (2 parameters)", 
      sig = "Sigmoidal growth (3 parameters)")
  
  return(out)
}

# ---------------------------------------------------------
# Parent function to determine best fitting function
# ---------------------------------------------------------
get_best_model = function(id, run, data, pb = NULL) {
  
  # Details of this run
  this_run = run[run_id == id]
  
  # Reduce data down to what we're interested in
  fit_data_dt = data %>%
    lazy_dt() %>%
    filter(country  == this_run$country,
           d_v_a_id == this_run$d_v_a_id) %>%
    select(x = fvps,
           y = impact) %>%
    # Multiply impact for more consistent x-y scales...
    mutate(y = y * o$impact_scaler) %>%
    as.data.table()
  
  # Append the origin (zero vaccine, zero impact)
  fit_data_dt = data.table(x = 1e-12, y = 0) %>%
    rbind(fit_data_dt)
  
  # Do not fit if insufficient data
  if (nrow(fit_data_dt) > 1) {
    
    # Declare x and y values for which we want to determine a relationship
    x = fit_data_dt$x
    y = fit_data_dt$y
    
    # Functions we'll attempt to fit with
    fns = fn_set()
    
    # Use optim algorithm to get good starting point for MLE
    start = credible_start(fns, x, y)
    
    # Run MLE from this starting point
    fit = run_mle(fns, start, x, y)
    
    # Determine AICc value for model suitability
    model_quality(fns, fit, x, y, id)
  }
  
  # Update progress bar
  if (!is.null(pb))
    setTxtProgressBar(pb, which(run$run_id == id))
}

# ---------------------------------------------------------
# Determine credible starting points for MLE - it needs it
# ---------------------------------------------------------
credible_start = function(fns, x, y, plot = FALSE) {
  
  # Initialise starting point for sigma in likelihood function
  s0 = 1  # This is essentially a placeholder until run_mle 
  
  # Let's start with any old points
  start = list(
    lin = list(s = s0, a = 1),
    log = list(s = s0, a = 1, b = 1), 
    exp = list(s = s0, a = 1, b = 1), 
    sig = list(s = s0, a = 1, b = 1, c = 1))
  
  # Define an objective function to minimise - sum of squares
  obj_fn = function(fn, ...) {
    
    # Squared difference
    diff_sq = (y - fns[[fn]](x, ...)) ^ 2
    
    # The sum of the squared difference
    obj_val = list(y = sum(diff_sq))
    
    return(obj_val)
  }
  
  # Define model-specific calls to objective function
  asd_fn = list(
    lin = function(p, args) obj_fn("lin", a = p[1]),
    log = function(p, args) obj_fn("log", a = p[1], b = p[2]),
    exp = function(p, args) obj_fn("exp", a = p[1], b = p[2]),
    sig = function(p, args) obj_fn("sig", a = p[1], b = p[2], c = p[3]))
  
  # Initiate list for storing plotting datatables
  plot_list = list()
  
  # Points to evaluate when plotting
  x_eval = seq(0, max(x), length.out = 101)
  
  # Iterate through stats models
  for (fn in names(fns)) {
    
    # Number of parameters for this model
    n_pars = length(start[[fn]]) - 1
    
    # Run ASD optimisation algorithm
    optim = asd(
      fn   = asd_fn[[fn]],
      x0   = rep(1, n_pars),
      lb   = rep(1e-10, n_pars),
      ub   = rep(Inf, n_pars),
      max_iters = 1e3)
    
    # Overwrite starting point with optimal parameters
    start[[fn]] = c(s0, optim$x) %>%
      setNames(names(start[[fn]])) %>%
      as.list()
    
    # Store diagnostic plotting data if desired
    if (plot == TRUE) {
      
      # Construct plotting datatable for this model
      plot_args = c(list(x_eval), as.list(optim$x))
      plot_list[[fn]] = data.table(
        x = x_eval,
        y = do.call(fns[[fn]], as.list(plot_args)), 
        fn = fn)
    }
  }
  
  # Create diagnostic plot if desired
  if (plot == TRUE) {
    
    # Plot the quality of fit for each model starting point
    g_asd = ggplot(rbindlist(plot_list)) +
      aes(x = x, y = y) +
      geom_line(aes(colour = fn)) + 
      geom_point(data = data.table(x = x, y = y),
                 colour = "black")
  }
  
  return(start)
}

# ---------------------------------------------------------
# Fit MLE for each fn using previously determined start point
# ---------------------------------------------------------
run_mle = function(fns, start, x, y, plot = FALSE) {
  
  # Log likelihood function to maximise
  likelihood = function(s, y_pred)
    ll = -sum(dnorm(x = y, mean = y_pred, sd = s, log = TRUE))
  
  # Annoyingly we need to name likelihood inputs, hence wrapper functions
  model = list(
    lin = function(s, a)       likelihood(s, fns$lin(x, a)),
    log = function(s, a, b)    likelihood(s, fns$log(x, a, b)), 
    exp = function(s, a, b)    likelihood(s, fns$exp(x, a, b)), 
    sig = function(s, a, b, c) likelihood(s, fns$sig(x, a, b, c)))
  
  # Initiate list for storing plotting datatables
  fit = plot_list = list()
  
  # Points to evaluate when plotting
  x_eval = seq(0, max(x), length.out = 101)
  
  # Iterate through stats models
  for (fn in names(fns)) {
    
    # Vector of lower bounds
    n_pars = length(start[[fn]]) - 1  # Zero for function coefficients
    l_bnds = c(1e-1, rep(0, n_pars))  # Small but non-zero for sigma
    
    # Catch errors when running in parallel
    if (!o$parallel$impact) {
      
      # Attempt to fit MLE model using prevriously determined start point
      fit_result = tryCatch(
        expr = mle(
          minuslogl = model[[fn]],
          start     = start[[fn]],
          lower     = l_bnds,
          nobs      = length(y)),
        
        # Return trivial if MLE failed
        error   = function(e) return(NULL),
        warning = function(w) return(NULL))
    }
    
    # Do not catch errors when running in parallel
    if (o$parallel$impact) {
      
      # Attempt to fit MLE model using prevriously determined start point
      fit_result = mle(
        minuslogl = model[[fn]], 
        start     = start[[fn]], 
        lower     = l_bnds,
        nobs      = length(y))
    }
    
    # Store results if we've had success
    if (!is.null(fit_result)) {
      
      # Calculate final log likelihood
      ll = do.call(model[[fn]], as.list(fit_result@coef))
      
      # Store fit with associated log likelihood 
      fit[[fn]] = list(mle = fit_result, ll = ll)
      
      # Store diagnostic plotting data if desired
      if (plot == TRUE) {
        
        # Extract coefficients for plotting
        coef = fit_result@coef[-1]
        
        # Construct plotting datatable for this model
        plot_args = c(list(x_eval), as.list(coef))
        plot_list[[fn]] = data.table(
          x = x_eval,
          y = do.call(fns[[fn]], as.list(plot_args)), 
          fn = fn)
      }
    }
  }
  
  # Create diagnostic plot if desired
  if (plot == TRUE) {
    
    # Plot the quality of fit for each model
    g_mle = ggplot(rbindlist(plot_list)) +
      aes(x = x, y = y) +
      geom_line(aes(colour = fn)) + 
      geom_point(data = data.table(x = x, y = y),
                 colour = "black")
  }
  
  return(fit)
}

# ---------------------------------------------------------
# Determine model quality - primarily this is via AICc
# ---------------------------------------------------------
model_quality = function(fns, fit, x, y, run_id) {
  
  # Return out if no fits succesful 
  if (length(fit) == 0)
    return()
  
  # Coefficients for each successful model
  coef = unlist(lapply(fit, function(a) a$mle@coef[-1]))
  coef = tibble(var = names(coef), value = coef) %>%
    separate(var, c("fn", "coef")) %>%
    mutate(run_id = run_id, 
           .before = 1) %>%
    as.data.table()
  
  # Save to file
  save_name = paste1(run_id, "coef")
  save_rds(coef, "runs", save_name)
  
  # Extract log likelihood...
  ll = lapply(fit, function(a) a$ll) %>%
    as.data.table() %>%
    pivot_longer(cols = everything(), 
                 values_to = "ll") %>%
    as.data.table()
  
  # Calculate AIC - adjusted for sample size
  aicc = sapply(fit, AICc) %>%
    as.list() %>%
    as.data.table() %>%
    pivot_longer(cols = everything(), 
                 values_to = "aicc") %>%
    # Append log likelihood...
    left_join(ll, by = "name") %>%
    rename(fn = name) %>%
    mutate(run_id = run_id, 
           .before = 1) %>%
    as.data.table()
  
  # Save to file
  save_name = paste1(run_id, "aicc")
  save_rds(aicc, "runs", save_name) 
}

# ---------------------------------------------------------
# Use AICc rather than AIC to reduce overfitting
# ---------------------------------------------------------
AICc = function(x) {
  
  # See en.wikipedia.org/wiki/Akaike_information_criterion
  
  # Sample size
  n = x$mle@nobs
  
  # Number of parameters
  k = length(x$mle@details$par) - 1
  
  # The usual AIC term
  aic_term = stats::AIC(x$mle)
  
  # An additional penalty term for small sample size
  pen_term = (2*k^2 + 2*k) / (n - k - 1)
  
  # Sum these terms
  aicc = aic_term # + pen_term
  
  return(aicc)
}

# ---------------------------------------------------------
# Compile all results into full datatables
# ---------------------------------------------------------
compile_results = function(run_dt) {
  
  message(" - Compiling all results")
  
  # Repeat for each type of output
  for (type in c("coef", "aicc")) {
    
    # All file names for this type of result
    names = paste1(run_dt$run_id, type)
    files = paste0(o$pth$runs, names, ".rds")
    
    # Files that exist 
    #
    # NOTE: Files with insufficient data are not saved
    result_files = files[file.exists(files)]
    
    # Load all results and squash into single datatable
    results_dt = result_files %>%
      lapply(readRDS) %>%
      rbindlist(fill = TRUE)
    
    # Append run details
    compiled_dt = run_dt %>%
      left_join(y  = results_dt, 
                by = "run_id") %>%
      select(-run_id)
    
    # Save to file
    save_rds(compiled_dt, "impact", type)
  }
}

# ---------------------------------------------------------
# Select best function considering complexity
# ---------------------------------------------------------
model_selection = function() {
  
  message(" - Selecting best functions")
  
  # Extract best fitting function based on AICc
  best_dt = read_rds("impact", "aicc") %>%
    lazy_dt() %>%
    rename(value = !!o$selection_metric) %>%
    # Select function with lowest AICc (or LL)...
    group_by(d_v_a_id, country) %>%
    slice_min(value, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    # Tidy up...
    select(d_v_a_id, country, fn) %>%
    as.data.table()
  
  # Save to file
  save_rds(best_dt, "impact", "best_model")
}

