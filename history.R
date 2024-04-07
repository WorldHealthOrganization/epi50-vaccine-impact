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
run_history = function(metric) {
  
  # Only continue if specified by do_step
  if (!is.element(6, o$do_step)) return()
  
  message("* Calculating impact of historical coverage: ", metric)
  
  # Eausy reference external d-v-a IDs
  extern_id = table("d_v_a")[source == "extern", d_v_a_id]
  
  # ---- Extract FVP data to be evaluated ----
  
  message(" > Preparing historical coverage")
  
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
  start_fit_dt = read_rds("impact", "impact", metric, "data") %>%
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
  
  message(" > Evaluating impact functions")
  
  # Parent function to evaluate impact functions
  evaluate_fn = function(data, uncert) {
    
    # Evaluate FVPs using impact functions
    result = data %>%
      select(d_v_a_id, country, year, pop, fvps = fvps_cum) %>%
      # Rename solely for use in evaluation fn...
      evaluate_impact_function(
        metric = metric, 
        uncert = uncert) %>%
      # Convert back to meaningful names...
      rename(impact_cum = impact) %>%
      # Reverse cumsum to derive annual relative impact...
      lazy_dt() %>%
      group_by(d_v_a_id, country, sample) %>%
      mutate(impact_rel = rev_cumsum(impact_cum)) %>%
      ungroup() %>%
      # Rescale back to population scale...
      mutate(impact = impact_rel * pop) %>%
      select(d_v_a_id, country, year, sample, impact) %>%
      as.data.table()
    
    return(result)
  }
  
  message("  - Best estimate coefficients")
  
  # Evaluate impact of relevant FVPs - best fit coefficients
  fit_best_dt = eval_dt %>%
    evaluate_fn(uncert = FALSE) %>%
    mutate(sample = "best")
  
  message("  - Posterior coefficients (", 
          o$uncertainty_samples, " samples)")
  
  # Evaluate impact of relevant FVPs - uncertainty samples
  fit_uncert_dt = eval_dt %>%
    evaluate_fn(uncert = TRUE)
    
  # Concatenate results
  fit_dt = rbind(fit_best_dt, fit_uncert_dt)
  
  # ---- Back project using initial impact ratios -----
  
  message(" > Back projecting")
  
  # Expand all FVP datatable for each uncertainty sample
  all_fvps_dt = fvps_dt %>%
    select(d_v_a_id, country, year, fvps = fvps_abs) %>%
    expand_grid(sample = unique(fit_dt$sample)) %>%
    as.data.table()
  
  # Initial impact per FVPs - used to back project
  #
  # NOTE: Idea behind init_impact_years is to smooth out any 
  #       initially extreme or jumpy impact ratios
  initial_ratio_dt = fit_dt %>%
    lazy_dt() %>%
    # Join associated FVPs...
    inner_join(y  = all_fvps_dt, 
               by = qc(d_v_a_id, country, year, sample)) %>%
    group_by(d_v_a_id, country, sample) %>%
    # Take the first init_impact_years years...
    slice_min(order_by  = year,
              n         = o$init_impact_years,
              with_ties = FALSE) %>%
    # Take the mean over these first years...
    summarise(impact_mean = mean(impact),
              fvps_mean   = mean(fvps)) %>%
    ungroup() %>%
    # Calculate the mean initial ratio...
    mutate(initial_ratio = impact_mean / fvps_mean) %>%
    select(d_v_a_id, country, sample, initial_ratio) %>%
    as.data.table()
  
  # Save initial ratio to file for diagnostic plotting
  save_rds(initial_ratio_dt, "history", "initial_ratio", metric) 
  
  # Back project by applying initial ratio 
  back_project_dt = fit_dt %>%
    lazy_dt() %>%
    # Join associated FVPs...
    full_join(y  = all_fvps_dt, 
              by = c("d_v_a_id", "country", "year", "sample")) %>%
    arrange(d_v_a_id, country, year, sample) %>%
    filter(!d_v_a_id %in% extern_id) %>%
    # Append impact ratio...
    left_join(y  = initial_ratio_dt, 
              by = c("d_v_a_id", "country", "sample")) %>%
    replace_na(list(initial_ratio = 0)) %>%
    # Ratio of past impact assumed consistent with initial years...
    mutate(impact = ifelse(
      test = is.na(impact), 
      yes  = fvps * initial_ratio, 
      no   = impact)) %>%
    select(d_v_a_id, country, year, sample, impact) %>%
    as.data.table()
  
  # ---- Append external models ----
  
  message(" > Appending external models")
  
  # Uncertainty sample from external models
  extern_uncertainty = table(paste1("extern_uncertainty", metric))
  
  # Concatenate for full set of samples
  all_samples = back_project_dt %>%
    rbind(extern_uncertainty)
  
  # Save all samples - allows cumulative summing in any direction
  save_rds(all_samples, "history", "all_samples", metric) 
  
  # ---- Summarise uncertainty ----
  
  message(" > Summarising uncertainty bounds")
  
  # Determine uncertainty bounds for temporal results
  result_time_dt = all_samples %>%
    summarise_uncertainty(cumulative = FALSE)  # See uncertainty.R
  
  # Equivalent uncertainty bounds for cumulative results
  #
  # NOTE: This special case is needed as it is not legitimate 
  #       to simply cumulatively sum temporal bounds
  # result_cum_dt = all_samples %>%
  #   summarise_uncertainty(cumulative = TRUE)  # See uncertainty.R
  
  # Save results to file
  save_rds(result_time_dt, "history", "burden_averted", metric) 
  # save_rds(result_cum_dt,  "history", "cumulative_averted", metric) 
  
  # ---- Supporting results ----
  
  # Only approporiate/necessary when computing deaths averted
  if (metric == "deaths") {
    
    # ---- Age structure of impact ----
    
    message(" > Summarising age structure")
    
    # Load all results and summarise age at impact
    age_impact_dt = table("vimc_estimates") %>%
      # Append all results...
      bind_rows(table("static_estimates")) %>%
      bind_rows(table("extern_estimates")) %>%
      select(d_v_a_id, country, year, age, impact = deaths_averted) %>%
      # Pregnancy vaccines are a special case...
      left_join(y  = table("d_v_a"), 
                by = "d_v_a_id") %>%
      mutate(impact = ifelse(
        test = grepl("_px$", vaccine) & age > 0, 
        yes  = 0, 
        no   = impact)) %>%
      # Attribute all impact to infants...
      mutate(impact = ifelse(
        test = grepl("_px$", vaccine) & age == 0, 
        yes  = 1, 
        no   = impact)) %>%
      # Summarise deaths over space and time...
      lazy_dt() %>%
      group_by(d_v_a_id, year, age) %>%
      summarise(value = sum(abs(impact))) %>%
      ungroup() %>%
      # Normalise absolute numbers...
      group_by(d_v_a_id, year) %>%
      mutate(scaler = value / sum(value)) %>%
      ungroup() %>%
      replace_na(list(scaler = 0)) %>%
      # Mean over time...
      group_by(d_v_a_id, age) %>%
      summarise(scaler = mean(scaler)) %>%
      ungroup() %>%
      # Normalise means...
      group_by(d_v_a_id) %>%
      mutate(scaler = scaler / sum(scaler)) %>%
      ungroup() %>%
      as.data.table()
    
    # Expand for each age
    age_effect_dt = table("d_v_a") %>%
      select(d_v_a_id) %>%
      expand_grid(age = o$ages) %>%
      left_join(y  = age_impact_dt, 
                by = c("d_v_a_id", "age")) %>%
      replace_na(list(scaler = 0)) %>%
      as.data.table()
    
    # g = age_effect_dt %>%
    #   group_by(d_v_a_id) %>%
    #   mutate(scaler = cumsum(scaler)) %>%
    #   ungroup() %>%
    #   format_d_v_a_name() %>%
    #   ggplot() +
    #   aes(x = age,
    #       y = scaler,
    #       colour = d_v_a_name) +
    #   geom_line() +
    #   facet_wrap(~d_v_a_name) +
    #   xlim(0, 50)
    
    # Save to tables cache
    save_table(age_effect_dt, "impact_age_multiplier")
    
    # ---- Convert deaths to YLL ----
    
    message(" > Calculating years of life lost")
    
    # Apply age structure and calculate YLL from deaths
    yll_samples_dt = all_samples %>%
      lazy_dt() %>%
      # Append life expectancy...
      left_join(y  = table("wpp_life_exp"), 
                by = c("country", "year")) %>%
      select(-age) %>%
      # Apply age structure...
      left_join(y  = table("impact_age_multiplier"), 
                by = "d_v_a_id", 
                relationship = "many-to-many") %>%
      mutate(deaths = impact * scaler) %>%
      # Calculate years of life lost...
      mutate(yll = deaths * pmax(0, life_exp - age)) %>%
      group_by(d_v_a_id, country, year, sample) %>%
      summarise(impact = sum(yll)) %>%
      ungroup() %>%
      as.data.table()
    
    # Save all samples - allows cumulative summing in any direction
    save_rds(yll_samples_dt, "history", "all_samples_yll") 
    
    message("  - Summarising YLL uncertainty")
    
    # Determine uncertainty bounds for temporal results
    yll_time_dt = yll_samples_dt %>%
      summarise_uncertainty(cumulative = FALSE)  # See uncertainty.R
    
    # Equivalent uncertainty bounds for cumulative results
    #
    # NOTE: This special case is needed as it is not legitimate 
    #       to simply cumulatively sum temporal bounds
    # yll_cum_dt = yll_samples_dt %>%
    #   summarise_uncertainty(cumulative = TRUE)  # See uncertainty.R
    
    # Save results to file
    save_rds(yll_time_dt, "history", "burden_averted_yll") 
    # save_rds(yll_cum_dt,  "history", "cumulative_averted_yll") 
  }
  
  # ---- Plot outcomes ----
  
  # Plot inital impact ratios used to back project
  plot_impact_fvps(metric, scope = "initial")
}

# ---------------------------------------------------------
# Evaluate impact function given FVPs
# ---------------------------------------------------------
evaluate_impact_function = function(data, metric, uncert = TRUE) {
  
  # ---- Evaluation functions ----
  
  # Function to valuate best coefficients
  eval_fn = function(i, sets, data, coef, fns) {
    
    # Index this set
    set = sets[i, ]
    
    # message(paste(id, collapse = ", "))
    
    # Exract FVPs for this ID
    data %<>% inner_join(set, by = qc(d_v_a_id, country))
    
    # Fitted function and parameters
    coef %<>% inner_join(set, by = qc(d_v_a_id, country, sample))
    
    # Load fitted function
    fn = fns[[unique(coef$fn)]]
    
    # Call function with fitted coefficients
    impact = fn(x = data$fvps, p = coef$value)
    
    # Output in datatable form
    result = cbind(data, impact)
    
    return(result)
  }
  
  # Wrapper function for evaluating all models for given country
  set_fn = function(sets, data, coef, fns) {
    
    # Country ID of this set
    id = unique(sets$country)
    
    # Apply evaluation function this set
    result_list = lapply(
      X    = seq_row(sets), 
      FUN  = eval_fn,
      sets = sets,
      data = data[country == id], 
      coef = coef[country == id], 
      fns  = fns)
    
    # Concatenate results
    result_dt = rbindlist(result_list)
    
    return(result_dt)
  }
  
  # ---- Sample coefficients ----
  
  # Load impact function posteriors
  posteriors = read_rds("impact", "posteriors", metric)
  
  # Generate uncertainty: sample posteriors
  if (uncert == TRUE) {
    
    # Sample uncertainty_samples sets
    samples = sample_vec(
      x    = unique(posteriors$iter), 
      size = o$uncertainty_samples)
    
    # Select associated coefficients
    coef = posteriors %>% 
      filter(iter %in% samples) %>%
      rename(sample = iter)
  }
  
  # No uncertainty: mean coefficients
  if (uncert == FALSE) {
    
    # Only one sample needed
    samples = 1
    
    # Take the mean of each coefficient
    coef = posteriors %>% 
      lazy_dt() %>%
      group_by(d_v_a_id, country, fn, param) %>% 
      summarise(value = mean(value)) %>% 
      ungroup() %>% 
      mutate(sample = samples, 
             .before = value) %>%
      arrange(d_v_a_id, country, param) %>% 
      as.data.table()
  }
  
  # Convert sample numbers to sample IDs
  sample_dict = get_sample_ids(samples)
  
  # ---- Perform evaluations ----
  
  # Load set of functions that may be evaluated
  fns = fn_set()
  
  # All sample-country-dva combos to evaluate
  sets = data %>%
    select(d_v_a_id, country) %>%
    unique() %>%
    expand_grid(sample = samples) %>%
    as.data.table() %>%
    split(.$country)
  
  # Apply evaluations in parallel
  if (o$parallel$history)
    result_list = mclapply(
      X    = sets, 
      FUN  = set_fn, 
      data = data, 
      coef = coef,
      fns  = fns,
      mc.cores = o$n_cores,
      mc.preschedule = FALSE)
  
  # Apply evaluations consecutively
  if (!o$parallel$history)
    result_list = lapply(
      X    = sets, 
      FUN  = set_fn,
      data = data, 
      coef = coef, 
      fns  = fns)
  
  # Squash results into single datatable
  result_dt = rbindlist(result_list) %>%
    # Transform impact to real scale...
    mutate(impact = impact / o$impact_scaler) %>%
    # Recode sample names for readability...
    mutate(sample = recode(sample, !!!sample_dict))
  
  return(result_dt)
}

# ---------------------------------------------------------
# Calculate child mortality rates in vaccine and no vaccine scenarios
# ---------------------------------------------------------
mortality_rates = function(age_bound = 0, grouping = "none") {
  
  # NOTE: Options for 'grouping' argument: "none", "region", or "income"
  
  # ---- Demography ----
  
  # Population as per WPP - needed to convert to rates
  pop_dt = table("wpp_pop") %>%
    filter(age <= age_bound) %>%
    group_by(country, year) %>%
    summarise(pop = sum(pop)) %>%
    ungroup() %>%
    as.data.table()
  
  # Child deaths as recorded by WPP
  deaths_dt = table("wpp_deaths") %>%
    filter(age <= age_bound) %>%
    group_by(country, year) %>%
    summarise(deaths = sum(deaths)) %>%
    ungroup() %>%
    as.data.table()
  
  # ---- Age-structured deaths averted ----

  # Vaccine impact disaggregated by age
  age_effect = table("impact_age_multiplier") %>%
    filter(age <= age_bound + 1) %>%
    group_by(d_v_a_id) %>%
    summarise(scaler = sum(scaler)) %>%
    ungroup() %>%
    as.data.table()
  
  # g = ggplot(age_effect %>%
  #              format_d_v_a_name()) +
  #   aes(fill = d_v_a_name,
  #       x    = d_v_a_name,
  #       y    = scaler) +
  #   geom_col() +
  #   theme(axis.text.x = element_text(
  #     angle = 50, hjust = 1))
  
  # Estimated child deaths averted by vaccination
  averted_dt = read_rds("history", "burden_averted_deaths") %>%
    left_join(y  = age_effect, 
              by = "d_v_a_id") %>%
    mutate(impact_age = impact * scaler) %>%
    # Summarise over d-v-a...
    group_by(country, year) %>%
    summarise(averted = sum(impact_age)) %>%
    ungroup() %>%
    as.data.table()
  
  # ---- Mortality rates ----
  
  # Construct grouping datatable to be joined to results
  grouping_dt = table("country") %>%
    left_join(y  = table("income_status"), 
              by = "country") %>%
    mutate(none = "none") %>%
    select(group = !!grouping, country, year)
  
  # Calculate mortality rates in each scenario
  mortality_dt = grouping_dt %>%
    left_join(y  = pop_dt, 
              by = c("country", "year")) %>%
    # Join child death estimates...
    left_join(y  = deaths_dt, 
              by = c("country", "year")) %>%
    left_join(y  = averted_dt, 
              by = c("country", "year")) %>%
    # Summarise over group...
    group_by(group, year) %>%
    summarise(pop     = sum(pop), 
              deaths  = sum(deaths), 
              averted = sum(averted)) %>%
    # Calculate child mortality rates...
    mutate(deaths_alt1 = deaths + averted, 
           rate        = deaths      / pop, 
           rate_alt1   = deaths_alt1 / pop) %>%
    group_by(group) %>%
    # Alternative scenario: no improvement in mortality rate...
    mutate(rate_alt2 = rate_alt1[1]) %>%
    ungroup() %>%
    mutate(deaths_alt2 = pop * rate_alt2, 
           .after = deaths_alt1) %>%
    select(-pop, -averted) %>%
    as.data.table()
  
# --- Double-counting check ---
  # Estimated child deaths averted by vaccination
  averted_dva_dt = read_rds("history", "burden_averted_deaths") %>%
    left_join(y  = age_effect, 
              by = "d_v_a_id") %>%
    mutate(impact_age = impact * scaler) %>%
    group_by(d_v_a_id, country, year) %>%
    summarise(averted = sum(impact_age)) %>%
    ungroup() %>%
    as.data.table()
  
  # Calculate mortality rates in each scenario
  mortality_dva_dt = grouping_dt %>%
    left_join(y  = pop_dt, 
              by = c("country", "year")) %>%
    # Join child death estimates...
    left_join(y  = deaths_dt, 
              by = c("country", "year")) %>%
    left_join(y  = averted_dva_dt, 
              by = c("country", "year")) %>%
    # Summarise over group...
    group_by(d_v_a_id, group) %>%
    summarise(pop     = sum(pop), 
              deaths  = sum(deaths), 
              averted = sum(averted)) %>%
    # Calculate child mortality rates...
    mutate(deaths_alt1 = deaths + averted, 
           rate        = deaths      / pop, 
           rate_alt1   = deaths_alt1 / pop) %>%
    select(-averted) %>%
    as.data.table()
  
  uncorrected = mortality_dva_dt %>%
                 mutate(diff = rate_alt1 - rate) %>%
                 select(diff) %>%
                 sum()
  
  bernoulli = mortality_dva_dt %>%
               mutate(diff = rate_alt1 - rate) %>%
               mutate(inverse = 1- diff) %>%
               select(inverse) %>%
               cumprod()
  
  bernoulli = 1-min(bernoulli)
  
  lower_bound_error = uncorrected-bernoulli
  lower_bound_deaths = mortality_dva_dt %>%
                        summarise(pop = mean(pop)) %>%
                        mutate(value = pop * lower_bound_error) %>%
                        select(-pop)
                       
  # ---- Format output ----
  
  # Use more descriptive scenario names
  scenarios = qc(vaccine, no_vaccine, no_other)
  col_names = c("group", "year", scenarios)
  
  # Function to format output datatable
  mortality_format_fn = function(dt, metric) {

    # Selection of metrics and melt to long format
    formated_dt = dt %>%
      select(group, year, starts_with(metric)) %>%
      as_named_dt(col_names) %>%
      pivot_longer(cols = -c(group, year), 
                   names_to = "scenario") %>%
      select(scenario, group, year, value) %>%
      arrange(scenario, group, year) %>%
      as.data.table()
    
    return(formated_dt)
  }
  
  # Store results in list
  mortality = list(
    rate  = mortality_format_fn(mortality_dt, "rate"),
    value = mortality_format_fn(mortality_dt, "deaths"))
  
  return(mortality)
}

