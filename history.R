###########################################################
# HISTORY
#
# Bring everything together: evaluate non-linear impact 
# functions to determine vaccine impact. Back project using
# ratio of first few years of impact estimates.
#
###########################################################

# ---------------------------------------------------------
# Use impact functions to calculate historical impact
# ---------------------------------------------------------
run_history = function() {
  
  # Only continue if specified by do_step
  if (!is.element(6, o$do_step)) return()
  
  message("* Calculating impact of historical coverage")
  
  # ---- Extract FVP data to be evaluated ----
  
  message(" - Preparing historical coverage")
  
  # Population size of each country over time
  pop_dt = table("wpp_pop") %>%
    lazy_dt() %>%
    group_by(country, year) %>%
    summarise(pop = sum(pop)) %>%
    ungroup() %>%
    as.data.table()
  
  # Extract FVPs over time
  #
  # NOTE: We do not need cumulative FVPs yet, these are only for evaluating
  #       impact fns, and we first need to subset what is to be evaluated.
  fvps_dt = table("coverage") %>%
    lazy_dt() %>%
    # Summarise over age...
    group_by(d_v_a_id, country, year) %>%
    summarise(fvps_abs = sum(fvps)) %>%
    ungroup() %>%
    # Scale results to per capita...
    left_join(y  = pop_dt, 
              by = c("country", "year")) %>%
    select(d_v_a_id, country, year, pop, fvps_abs) %>%
    mutate(fvps_rel = fvps_abs / pop) %>%
    as.data.table()
  
  # From which years have impact functions been fit from
  start_fit_dt = read_rds("impact", "data") %>%
    lazy_dt() %>%
    select(d_v_a_id, country, year) %>%
    group_by(d_v_a_id, country) %>%
    slice_min(year, with_ties = FALSE) %>%
    ungroup() %>%
    rename(start_year = year) %>%
    as.data.table()
  
  # FVPs to evaulate using impact functions 
  eval_dt = 
    # Begin for full factorial of points...
    expand_grid(
      d_v_a_id  = table("d_v_a")$d_v_a_id, 
      country   = all_countries(), 
      eval_year = o$years) %>%
    # Remove years prior to fit start year...
    inner_join(y  = start_fit_dt, 
               by = c("d_v_a_id", "country")) %>%
    filter(eval_year >= start_year) %>%
    select(d_v_a_id, country, year = eval_year) %>%
    # Append FVPs...
    inner_join(y  = fvps_dt, 
               by = c("d_v_a_id", "country", "year")) %>%
    # Cumulative sum FVPs...
    arrange(d_v_a_id, country, year) %>%
    group_by(d_v_a_id, country) %>%
    mutate(fvps_cum = cumsum(fvps_rel)) %>%
    ungroup() %>%
    as.data.table()
  
  # ---- Evaluate impact functions -----
  
  message(" - Evaluating impact functions")
  
  # Evaluate impact of relevant FVPs
  result_fit_dt = eval_dt %>%
    # Rename solely for use in evaluation fn...
    rename(fvps = fvps_cum) %>%
    evaluate_impact_function() %>%
    # Convert back to meaningful names...
    rename(fvps_cum   = fvps, 
           impact_cum = impact) %>%
    # Reverse cumsum to derive annual relative impact...
    lazy_dt() %>%
    group_by(d_v_a_id, country) %>%
    mutate(impact_rel = rev_cumsum(impact_cum)) %>%
    ungroup() %>%
    # Rescale back to population scale...
    mutate(impact_abs = impact_rel * pop) %>%
    as.data.table()
  
  # ---- Back project using initial impact ratios -----
  
  message(" - Back projecting")
  
  # Initial impact per FVPs - used to back project
  #
  # NOTE: Idea behind init_impact_years is to smooth out any 
  #       initially extreme or jumpy impact ratios
  initial_ratio_dt = result_fit_dt %>%
    lazy_dt() %>%
    group_by(d_v_a_id, country) %>%
    # Take the first init_impact_years years...
    slice_min(order_by  = year, 
              n         = o$init_impact_years, 
              with_ties = FALSE) %>%
    # Take the mean over these first years...
    summarise(impact_mean = mean(impact_abs), 
              fvps_mean   = mean(fvps_abs)) %>%
    ungroup() %>%
    # Calculate the mean initial ratio...
    mutate(initial_ratio = impact_mean / fvps_mean) %>%
    select(d_v_a_id, country, initial_ratio) %>%
    as.data.table()
  
  # Save initial ratio to file for diagnostic plotting
  save_rds(initial_ratio_dt, "impact", "initial_ratio") 
  
  # Back project by applying initial ratio 
  back_project_dt = result_fit_dt %>%
    lazy_dt() %>%
    select(d_v_a_id, country, year, impact_abs) %>%
    # Join with full FVPs data...
    full_join(y  = fvps_dt, 
              by = c("d_v_a_id", "country", "year")) %>%
    select(d_v_a_id, country, year, 
           fvps   = fvps_abs, 
           impact = impact_abs) %>%
    arrange(country, d_v_a_id, year) %>%
    # Append impact ratio...
    left_join(y  = initial_ratio_dt, 
              by = c("d_v_a_id", "country")) %>%
    replace_na(list(initial_ratio = 0)) %>%
    # Ratio of past impact assumed consistent with initial years...
    mutate(impact = ifelse(
      test = is.na(impact), 
      yes  = fvps * initial_ratio, 
      no   = impact)) %>%
    select(-initial_ratio) %>%
    as.data.table()
  
  # ---- Append external models ----
  
  message(" - Appending external models")
  
  # Load up results from external models
  extern_dt = table("extern_estimates") %>%
    lazy_dt() %>%
    filter(d_v_a_id %in% table("d_v_a")$d_v_a_id) %>%
    # Summarise results over age...
    group_by(d_v_a_id, country, year) %>%
    summarise(impact = sum(deaths_averted)) %>%
    ungroup() %>%
    # TODO: Update placeholder with actual values
    mutate(fvps = 1, .before = impact) %>%
    as.data.table()
  
  # Concatenate results from external models
  result_dt = back_project_dt %>%
    rbind(extern_dt)
  
  # Save results to file
  save_rds(result_dt, "history", "deaths_averted") 
  
  # ---- Plot results ----
  
  # TODO: Add new plots here...
  
  # Plot inital impact ratios used to back project
  plot_impact_fvps(scope = "initial")
  
  # Plot primary results figure: historical impact over time
  plot_historical_impact()
  
  # Plot change in child mortality rates over time
  plot_child_mortality()
}

# ---------------------------------------------------------
# Evaluate impact function given FVPs
# ---------------------------------------------------------
evaluate_impact_function = function(eval_dt = NULL) {
  
  # ---- Load stuff ----
  
  # Load results: best fit functions and associtaed coefficients
  best_dt = read_rds("impact", "best_model")
  coef_dt = read_rds("impact", "coef")
  
  # Best coefficients
  best_coef = coef_dt %>%
    lazy_dt() %>%
    inner_join(y  = best_dt, 
               by = c("d_v_a_id", "country", "fn")) %>%
    mutate(par = as.list(setNames(value, coef))) %>% 
    group_by(d_v_a_id, country, fn) %>% 
    summarise(par = list(par)) %>% 
    ungroup() %>% 
    arrange(d_v_a_id, country) %>% 
    as.data.table()
  
  # ---- Interpret trivial argument ----
  
  # Countrues to evaluate
  if (is.null(eval_dt)) {
    
    # Default points at which to evaluate
    x_eval = seq(0, o$eval_x_scale, length.out = 101)
    
    # Construct evaluation datatable
    eval_dt = best_coef %>%
      select(d_v_a_id, country) %>%
      expand_grid(fvps = x_eval) %>%
      as.data.table()
  }
  
  # ---- Evaluate best fit model and coefficients ----
  
  # Function to valuate best coefficients
  eval_fn = function(i) {
    
    # message(paste(c_d_v_a[i, ], collapse = ", "))
    
    # Exract FVPs for this c_d_v_a
    data = c_d_v_a[i, ] %>%
      inner_join(y  = eval_dt, 
                 by = c("d_v_a_id", "country"))
    
    # Fitted function and parameters
    pars = c_d_v_a[i, ] %>%
      inner_join(y  = best_coef, 
                 by = c("d_v_a_id", "country"))
    
    # Load fitted function
    fn = fn_set()[[pars$fn]]
    
    # Call function with fitted coefficients
    impact = do.call(fn, c(list(x = data$fvps), pars$par[[1]]))
    
    # Output result in datatable form
    result = cbind(data, impact)
    
    return(result)
  }
  
  # All country - dva combos to evaluate
  c_d_v_a = eval_dt %>%
    inner_join(y  = best_coef,
               by = c("d_v_a_id", "country")) %>%
    select(d_v_a_id, country) %>%
    unique()
  
  # Apply the evaluation funtion
  result_dt = seq_row(c_d_v_a) %>%
    lapply(eval_fn) %>%
    rbindlist() %>%
    # Transform impact to real scale...
    mutate(impact = impact / o$impact_scaler)
  
  return(result_dt)
}

# ---------------------------------------------------------
# Calculate child mortality rates in vaccine and no vaccine scenarios
# ---------------------------------------------------------
mortality_rates = function(age_bound = 5, grouping = "none") {
  
  # NOTE: Options for 'grouping' argument: "none", "region", or "income"
  
  # Construct grouping datatable to be joined to results
  grouping_dt = table("country") %>%
    left_join(y  = table("income_status"), 
              by = "country") %>%
    mutate(none = "none") %>%
    select(country, year, group = !!grouping)
  
  # Child deaths as recorded by WPP
  deaths_dt = table("wpp_deaths") %>%
    filter(age <= age_bound) %>%
    left_join(y  = grouping_dt, 
              by = c("country", "year")) %>%
    group_by(group, year) %>%
    summarise(deaths = sum(deaths)) %>%
    ungroup() %>%
    arrange(group, year) %>%
    as.data.table()
  
  # Vaccine impact disaggregated by age
  age_effect = impact_age_multiplier()
  
  # Estimated child deaths averted by vaccination
  averted_dt = read_rds("history", "deaths_averted") %>%
    expand_grid(age_effect) %>%
    filter(age <= age_bound) %>%
    left_join(y  = grouping_dt, 
              by = c("country", "year")) %>%
    group_by(group, year) %>%
    summarise(averted = sum(impact * scaler)) %>%
    ungroup() %>%
    arrange(group, year) %>%
    as.data.table()
  
  # Population as per WPP - needed to convert to rates
  pop_dt = table("wpp_pop") %>%
    filter(age <= age_bound) %>%
    left_join(y  = grouping_dt, 
              by = c("country", "year")) %>%
    group_by(group, year) %>%
    summarise(pop = sum(pop)) %>%
    ungroup() %>%
    arrange(group, year) %>%
    as.data.table()
  
  # Calculate mortality rates in each scenario
  mortality_dt = pop_dt %>%
    # Join child death estimates...
    left_join(y  = deaths_dt, 
              by = c("group", "year")) %>%
    left_join(y  = averted_dt, 
              by = c("group", "year")) %>%
    replace_na(list(averted = 0)) %>%
    # Calculate child mortality rates...
    mutate(no_vaccine = (deaths + averted) / pop, 
           vaccine    = deaths / pop) # %>%
    # select(group, year, vaccine, no_vaccine)
  
  return(mortality_dt)
}

# ---------------------------------------------------------
# Scaler to estimate vaccination impact disaggregated by age
# ---------------------------------------------------------
impact_age_multiplier = function() {
  
  # TODO: Extract impact distribution by age for each d-v-a
  #       using original impact estimates
  
  # TEMP: Assume impact is highest for infants and decays exponentially
  age_effect = data.table(age = o$ages) %>%
    mutate(scaler = exp(-sqrt(2 * age)), 
           scaler = scaler / sum(scaler)) 
  
  return(age_effect)
}

