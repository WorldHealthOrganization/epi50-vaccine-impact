###########################################################
# UNCERTAINTY
#
# Functions for generating and summarising model uncertainty.
#
###########################################################

# ---------------------------------------------------------
# Generate samples from extern model results
# ---------------------------------------------------------
extern_uncertainty = function() {
  
  message(" > Generating uncertainty (", 
          o$uncertainty_samples, " samples)")
  
  # ---- External model uncertainty summary ----
  
  # Function to load external model uncertainty
  load_fn = function(model) {
    
    # Load external model uncertainty where available
    model_uncert_dt = try_load(
      pth  = o$pth$extern, 
      file = paste1("epi50", model, "uncertainty"), 
      throw_error = FALSE) 
    
    return(model_uncert_dt)
  }
  
  # All external models
  model_dt = table("extern_models")
  
  # Average uncertainty bounds per d-v-a per year
  bounds_dt = model_dt %>%
    pull(model) %>%
    lapply(load_fn) %>%
    rbindlist() %>%
    left_join(y  = model_dt, 
              by = "model") %>%
    left_join(y  = table("d_v_a"), 
              by = "disease") %>%
    group_by(d_v_a_id, year) %>%
    summarise(lb = pmin(1, mean(lower)), 
              ub = pmax(1, mean(upper))) %>%
    ungroup() %>%
    as.data.table()
  
  # ---- Generate uncertianty samples -----
  
  # IDs of samples to generate
  sample_ids = get_sample_ids(1 : o$uncertainty_samples)
  
  # Load best estimate outcomes (all metrics)
  best_dt = table("extern_estimates")
  
  # Consider one metric at a time
  for (metric in o$metrics) {
    
    # Best estimate for this metric
    metric_dt = best_dt %>%
      rename(value = !!paste1(metric, "averted")) %>%
      lazy_dt() %>%
      # Summarise over age...
      group_by(d_v_a_id, country, year) %>%
      summarise(impact = sum(value)) %>%
      ungroup() %>%
      as.data.table()
    
    # Apply bounds and generate uncertainty samples
    samples_dt = metric_dt %>%
      left_join(y  = bounds_dt, 
                by = c("d_v_a_id", "year")) %>%
      # Assume trivial uncertainty if undefined...
      pivot_longer(cols = c(lb, ub), 
                   names_to  = "direction", 
                   values_to = "scaler") %>%
      replace_na(list(scaler = 1)) %>%
      # Apply uncertainty bounds...
      mutate(bound = impact * scaler, 
             diff  = abs(impact - bound)) %>%
      group_by(d_v_a_id, country, year) %>%
      slice_min(diff, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      # Interpret bound represents 3 standard deviations...
      mutate(sd = diff / 3) %>%
      select(d_v_a_id, country, year, mean = impact, sd) %>%
      # Sample from Gaussian uncertainty_samples times...
      expand_grid(sample = sample_ids) %>%
      mutate(impact = rnorm(n(), mean, sd)) %>%
      select(d_v_a_id, country, year, sample, impact) %>%
      as.data.table()
    
    # Concatenate with best estimate
    uncert_dt = metric_dt %>%
      mutate(sample = "best", 
             .before = impact) %>%
      rbind(samples_dt) %>%
      arrange(d_v_a_id, country, year, sample)
    
    # Save in tables cache
    save_table(uncert_dt, paste1("extern_uncertainty", metric))
  }
}

# ---------------------------------------------------------
# Summarise over posterior samples for lower and upper bounds
# ---------------------------------------------------------
summarise_uncertainty = function(data, cumulative) {
  
  # Grouping - may or may not include country
  grouping = intersect(
    x = names(data), 
    y = qc(d_v_a_id, disease, region, country, year))
  
  # Cumulative summing must be done before sample summary
  if (cumulative == TRUE) {
    
    # Grouping for cumsumming over time
    group_cum = c(setdiff(grouping, "year"), "sample")
    
    # Cumulatively sum over time first
    data %<>%
      lazy_dt() %>%
      group_by(across(all_of(group_cum))) %>%
      mutate(impact = cumsum(impact)) %>%
      ungroup() %>%
      as.data.table()
  }
  
  # Compute quantiles of all samples to obtain bounds
  bounds_dt = data %>%
    lazy_dt() %>%
    filter(sample != "best") %>%
    group_by(across(all_of(grouping))) %>%
    summarise(lower = quantile(impact, o$quantiles[1]),
              upper = quantile(impact, o$quantiles[2])) %>%
    ungroup() %>%
    as.data.table()
  
  # Append bounds to best estimate results
  summary_dt = data %>%
    filter(sample == "best") %>%
    select(-sample) %>%
    left_join(y  = bounds_dt,
              by = grouping)
  
  return(summary_dt)
}

# ---------------------------------------------------------
# Convert sample number to sample ID
# ---------------------------------------------------------
get_sample_ids = function(samples) {
  
  # For readability, we'll convert from sample numbers to IDs
  sample_names = sprintf("s%04i", seq_along(samples))
  
  # Construct named vector dictionary for recoding
  sample_dict = setNames(sample_names, samples)
  
  return(sample_dict)
}

