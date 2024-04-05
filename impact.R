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
run_impact = function(metric) {
  
  # Only continue if specified by do_step
  if (!is.element(5, o$do_step)) return()
  
  message("* Fitting impact functions: ", metric)
  
  # ---- FVPs and impact estimates ----

  message(" > Preparing FVP-impact data")

  # Prepare impact-FVP data to fit to
  data_dt = get_impact_data(metric)

  # Exploratory plots of data used to fit impact functions
  plot_impact_data(metric)

  # Plot all-time impact per FVPs
  # plot_impact_fvps(metric, scope = "all_time")

  # ---- Model fitting ----

  message(" > Evaluating impact functions")
  
  # Country-disease-vaccine-activity combinations
  run_dt = data_dt %>%
    select(d_v_a_id, country) %>%
    unique() %>%
    mutate(run_id = paste1(d_v_a_id, country))
  
  # Iterate through d-v-a one by one
  for (id in unique(run_dt$d_v_a_id)) {
    
    # Details of this d_v_a
    d_v_a_name = data.table(d_v_a_id = id) %>%
      format_d_v_a_name() %>%
      pull(d_v_a_name)
    
    # Display progress message to user
    message("  - ", d_v_a_name)
    
    # Subset what to run, and data to use
    run  = run_dt[d_v_a_id == id]
    data = data_dt[d_v_a_id == id]

    # Initiate progress bar
    pb = start_progress_bar(nrow(run))

    # Run get_best_model
    results_list = lapply(
      X    = run$run_id,
      FUN  = get_best_model,
      run  = run,
      data = data,
      pb   = pb)

    # Squash results into single datatable
    results_dt = rbindlist(results_list)

    # Save to file
    save_rds(results_dt, "impact", "impact", metric, id)

    # Close connections opened by sink
    closeAllConnections()
  }
  
  # ---- Model selection ----

  # Select best function for each country-d_v_a combination
  model_selection(run_dt, metric)

  # ---- Plot results ----

  # Plot function selection statistics
  plot_model_selection(metric)

  # Plot impact function evaluation
  plot_model_fits(metric)
  
  # Plot impact vs coverage by vaccine, income, and decade
  # plot_impact_coverage(metric)
}

# ---------------------------------------------------------
# Prepare impact-FVP data to fit to
# ---------------------------------------------------------
get_impact_data = function(metric) {
  
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
  vimc_dt = read_rds("impute", "impute", metric, "result", err = FALSE)
  
  # Load static model impact estimates
  static_dt = read_rds("static", metric, "averted_vaccine", err = FALSE) %>%
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
  save_rds(data_dt, "impact", "impact", metric, "data") 
  
  return(data_dt)
}

# ---------------------------------------------------------
# Set of functions to fit - we'll determine the 'best'
# ---------------------------------------------------------
fn_set = function(params = FALSE, dict = FALSE) {
  
  # Set of statistical models / functions we want to test
  out = list(
    lin = function(x, p) y = x * p[1],
    log = function(x, p) y = logarithmic_growth(x, p[1], p[2]),
    exp = function(x, p) y = exponential_growth(x, p[1], p[2]),
    sig = function(x, p) y = sigmoidal_growth(x, p[1], p[2], p[3]))
  
  # Alternative functionality - return number of params
  if (params == TRUE)
    out = c(lin = 1, log = 2, exp = 2, sig = 3)
  
  # Alternative functionality - return dictionary
  if (dict == TRUE)
    out = c(
      lin = "Linear gradient (1 parameter)", 
      log = "Logarithmic growth (2 parameters)",
      exp = "Exponential growth (2 parameters)", 
      sig = "Sigmoidal growth (3 parameters)")
  
  return(out)
}

# ---------------------------------------------------------
# Parent function to determine best fitting function
# ---------------------------------------------------------
get_best_model = function(id, run, data, pb) {
  
  # Initiate trivial output
  result = NULL
  
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
  
  # Number of genuine data points (aside from the origin)
  n_data = nrow(fit_data_dt) - 1 
  
  # Do not fit if insufficient data
  if (n_data >= 1) {

    # Declare x and y values for which we want to determine a relationship
    x = fit_data_dt$x
    y = fit_data_dt$y

    # Functions we'll attempt to fit with
    fns = fn_set()[fn_set(params = TRUE) <= n_data]

    # Attempt to determine global minimum for each function
    optim = run_optim(fns, x, y)
    
    # Apply MCMC using assumed global optimum as strong prior
    fit = run_mcmc(fns, optim, x, y)
    
    # Determine AICc value for model suitability
    result = model_quality(fns, fit, x, y, id)
  }
  
  # Update progress bar
  pb$tick()
  
  # ---- Diagnostic plot ----
  
  # data_dt   = data.table(x = x, value = y)
  # models_dt = data.table(x = x)
  # 
  # for (fn in names(fit))
  #   models_dt[[fn]] = fns[[fn]](x, fit[[fn]]$coef)
  # 
  # models_dt %<>%
  #   pivot_longer(cols = -x) %>%
  #   select(fn = name, x, value) %>%
  #   arrange(fn, x) %>%
  #   as.data.table()
  # 
  # plot_dt = result %>%
  #   filter(param == "ll") %>%
  #   mutate(lab = paste0(
  #     fn, "\nll = ",
  #     round(value, 2))) %>%
  #   select(fn, lab) %>%
  #   left_join(y  = models_dt,
  #             by = "fn")
  # 
  # g = ggplot(plot_dt) +
  #   aes(x = x, y = value) +
  #   geom_line(
  #     mapping = aes(colour = lab)) +
  #   geom_point(
  #     data   = data_dt,
  #     colour = "black")
  
  return(result)
}

# ---------------------------------------------------------
# Attempt to determine global minimum for each function
# ---------------------------------------------------------
run_optim = function(fns, x, y) {
  
  # Define an objective function to minimise - sum of squares
  obj_fn = function(p, fn) {
    
    # Squared difference
    diff_sq = (y - fns[[fn]](x, p)) ^ 2
    
    # The sum of the squared difference
    obj_val = list(y = sum(diff_sq))
    
    return(obj_val)
  }
  
  # Initiate optimal results list
  optim = list()
  
  # Iterate through stats models
  for (fn in names(fns)) {
    
    # Number of parameters for this model
    n_pars  = fn_set(params = TRUE)[[fn]]
    par_ref = letters[1 : n_pars]
    
    # Set lower and upper parameter bounds
    lb = rep(1e-10, n_pars)
    ub = rep(1e3, n_pars)
    
    # Inititae list to store results
    asd_results = list()
    
    # Repeat optimisation several times
    for (i in 1 : o$n_optim) {
      
      # Different starting point each time
      x0 = pmax(runif(n_pars), lb)

      # Run ASD optimisation algorithm
      asd_result = asd(
        fn   = obj_fn,
        args = fn,
        x0   = x0,
        lb   = lb,
        ub   = ub,
        iters = 200)
      
      # Store result and optimal parameters
      asd_results[[i]] = c(
        asd_result$y, 
        asd_result$x)
    }
    
    # Select best fitting parameters for this function
    optim[[fn]] = do.call(rbind, asd_results) %>%
      as_named_dt(c("y", par_ref)) %>%
      # Sort by objective function value...
      arrange(y) %>%
      slice_head(n = 1) %>%
      select(-y) %>%
      as.list()
  }
  
  return(optim)
}

# ---------------------------------------------------------
# Apply MCMC using assumed global optimum as strong prior
# ---------------------------------------------------------
run_mcmc = function(fns, optim, x, y) {
  
  # Log-likelihood function
  likelihood_fn = function(p) {
    
    # Set poor likelihood when parameter bounds are violated
    if (any(p < 1e-10)) 
      return(-1e6)
    
    # Evaluate model emulator for given parameters
    y_pred = fns[[fn]](x, p)
    
    # Calculate the log-likelihood
    ll = dnorm(
      x    = y_pred, 
      mean = y, 
      sd   = sd(y - y_pred), 
      log  = TRUE)
    
    # Calculate log-prior for all parameters
    lp = dnorm(
      x    = p, 
      mean = x0, 
      sd   = x0 * o$prior_sd, 
      log  = TRUE)
    
    # Weighting to be applied to priors
    #
    # NOTE: Dividing by number of parameters such that more complex
    #       models are not double punished when computing AICc
    prior_weight = o$prior_weight / length(p)
    
    # Sum and appply weighting to priors
    likelihood = sum(ll) + sum(lp) * prior_weight
    
    return(likelihood)
  }
  
  # Wrapper function for MCMC call
  mcmc_fn = function() {
    
    # Call Metropolis-Hasting algorithm
    mcmc_result = MCMCmetrop1R(
      fun    = likelihood_fn, 
      burnin = o$mcmc_burnin,
      mcmc   = o$mcmc_iter,
      thin   = o$mcmc_iter / o$mcmc_samples,
      tune   = 1.5,
      seed   = 1,
      theta.init   = x0, 
      optim.method = "L-BFGS-B",
      optim.lower  = 1e-10)
    
    return(mcmc_result)
  }
  
  # We'll send noisy output to a null file
  sink(nullfile())
  
  # Initiate results list
  fit = list()
  
  # Iterate through stats models
  for (fn in names(fns)) {
    
    # Parameter reference
    par_ref = names(optim[[fn]])
    
    # Assumed global minimum
    x0 = unlist(optim[[fn]])
    
    # Wrap MCMC call in try catch in case of errors
    mcmc_result = tryCatch(
      expr  = suppressWarnings(mcmc_fn()),
      error = function(e) return())
    
    # Store unless null result
    if (!is.null(mcmc_result)) {
      
      # Format resulting chain (burn-in already discarded)
      mcmc_chain = mcmc_result %>%
        as_named_dt(names(optim[[fn]]))
      
      # Take mean of posteriors as best fitting coefficients
      coef = colMeans(mcmc_result) %>%
        setNames(par_ref)
      
      # Store fit with associated log likelihood
      fit[[fn]] = list(
        chain = mcmc_chain,
        coef  = coef,
        ll    = likelihood_fn(coef))
    }
  }
  
  # Sink the output
  sink()
  
  return(fit)
}

# ---------------------------------------------------------
# Determine model quality - primarily this is via AICc
# ---------------------------------------------------------
model_quality = function(fns, fit, x, y, run_id) {
  
  # Return out if no fits succesful 
  if (length(fit) == 0)
    return()
  
  # ---- Model selection metrics ----
  
  # Calculate AIC - adjusted for sample size
  aicc = sapply(fit, aicc, n = length(y)) %>%
    as.list() %>%
    as.data.table() %>%
    pivot_longer(cols = everything(), 
                 names_to  = "fn") %>%
    mutate(param = "aicc") %>%
    as.data.table()
  
  # Extract log likelihood
  ll = lapply(fit, function(a) a$ll) %>%
    as.data.table() %>%
    pivot_longer(cols = everything(), 
                 names_to  = "fn") %>%
    mutate(param = "ll") %>%
    as.data.table()
  
  # ---- Model parameters ----
  
  # Coefficients for each model
  coef = unlist(lapply(fit, function(a) a$coef))
  coef = tibble(var = names(coef), value = coef) %>%
    separate(var, c("fn", "param")) %>%
    as.data.table()
  
  # ---- Parameter posteriors ----
  
  # MCMC chains for each model
  chains = lapply(fit, function(a) a$chain) %>%
    as.data.table() %>%
    mutate(iter = 1 : n()) %>%
    pivot_longer(cols = -iter) %>%
    separate(col  = "name", 
             into = c("fn", "param"), 
             fill = "right") %>%
    replace_na(list(param = "a")) %>%
    select(fn, param, iter, value) %>%
    arrange(fn, param, iter) %>%
    as.data.table()
  
  # ---- Concatenate output ----
  
  # Squash all details into single datatable()
  quality_dt = bind_rows(aicc, ll, coef, chains) %>%
    mutate(run_id = run_id) %>%
    select(run_id, fn, param, iter, value)
  
  return(quality_dt)
}

# ---------------------------------------------------------
# Use AICc rather than AIC to reduce overfitting
# ---------------------------------------------------------
aicc = function(x, n) {
  
  # See en.wikipedia.org/wiki/Akaike_information_criterion
  
  # Number of parameters
  k = length(x$coef)
  
  # Log likelihood associated with these parameters
  l = x$ll
  
  # The usual AIC term
  aic_term = 2*k - 2*l
  
  # An additional penalty term for small sample size
  pen_term = (2*k^2 + 2*k) / (n - k - 1)
  
  # Sum these terms
  aicc = aic_term + pen_term
  
  return(aicc)
}

# ---------------------------------------------------------
# Select best function considering complexity
# ---------------------------------------------------------
model_selection = function(run_dt, metric) {
  
  message(" > Selecting best functions")
  
  # ---- Extract results ----
  
  # All d-v-a combinations considered
  d_v_a = unique(run_dt$d_v_a_id)
  
  # Construct paths to results files
  names = paste1("impact", metric, d_v_a)
  files = paste0(o$pth$impact, names, ".rds")
  
  # Extract best fitting function based on AICc
  results_dt = lapply(files, read_rds) %>%
    rbindlist() %>%
    left_join(y  = run_dt, 
              by = "run_id")
  
  # ---- Model selection ----
  
  # Select best model based on AICc or LL
  selection_dt = results_dt %>%
    filter(param %in% c("aicc", "ll")) %>%
    # Transform log-likelihood so we search for the lowest...
    mutate(value = ifelse(param == "ll", -value, value)) %>%
    # Select models according to AICc and LL...
    group_by(d_v_a_id, country, param) %>%
    slice_min(value, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    # Model selection according to o$selection_metric...
    filter(param == o$selection_metric) %>%
    # Tidy up...
    select(d_v_a_id, country, fn) %>%
    as.data.table()
  
  # Save to file
  save_rds(selection_dt, "impact", "model_choice", metric)
  
  # ---- Posterior chains ----
  
  # Select best model based on AICc or LL
  posteriors_dt = results_dt %>%
    filter(iter > 0) %>%
    inner_join(y  = selection_dt, 
               by = c("d_v_a_id", "country", "fn")) %>%
    select(d_v_a_id, country, fn, param, iter, value)
  
  # Save to file
  save_rds(posteriors_dt, "impact", "posteriors", metric)
}

