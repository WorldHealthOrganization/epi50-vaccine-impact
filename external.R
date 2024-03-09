###########################################################
# EXTERNAL
#
# Prepare, simulate (if appropriate), and extract results
# for external polio and measles models.
#
###########################################################

# ---------------------------------------------------------
# Parent function for extracting results for all external models
# ---------------------------------------------------------
run_external = function() {
  
  # Only continue if specified by do_step
  if (!is.element(2, o$do_step)) return()
  
  message("* Preparing external models")
  
  # TODO: Split impact across each vaccine...
  
  # # Create templates for measles and polio models
  # template_measles()
  # template_polio()
  # 
  # # Simulate DynaMICE measles model
  # simulate_dynamice()

  # Read in external polio results (only required once)
  read_polio()  # TODO: Remove when polio results available
  read_penn_state()

  # Format external modelling results for EPI50 use
  format_measles()
  format_polio()

  # Extract results from all external models
  extract_extern_results()
  
  # ---- Data visualisation plots ----
  
  # Plot total number of FVP over time
  plot_total_fvps()
  
  # Plot coverage density by disease
  plot_coverage()
  
  # Coverage data density by age
  plot_coverage_age_density()
  
  # Missing coverage data by country
  plot_missing_data()
}

# ---------------------------------------------------------
# Create results template for measles models
# ---------------------------------------------------------
template_measles = function() {
  
  message(" - Creating results template: measles")
  
  # All metrics of intrerest (epi outcomes and number of doses)
  metrics = qc(deaths, dalys, MCV1_doses, MCV2_doses, SIA_doses)
  
  # Full factorial of factors: scenario, country, year, age, metric
  template_dt = 
    expand_grid(
      scenario = c("no_vaccine", "vaccine"), 
      country  = all_countries(), 
      year     = o$years, 
      age      = o$ages,
      metric   = metrics) %>%
    # Append random result placeholders...
    mutate(value = runif(n())) %>%
    as.data.table()
  
  # Write to file
  file = paste0(o$pth$extern, "template_measles.csv")
  fwrite(template_dt, file = file)
}

# ---------------------------------------------------------
# Create results template for polio model
# ---------------------------------------------------------
template_polio = function() {
  
  message(" - Creating results template: polio")
  
  # All metrics of intrerest (epi outcomes and number of doses)
  metrics = qc(paralytic_cases, deaths, dalys, OPV_doses, IPV_doses)
  
  # Two different stratifications of setting
  geo = list(
    region = sort(unique(table("country")$region)), 
    income = sort(unique(table("income_status")$income)))
  
  # Repeat process for setting stratification
  for (setting in names(geo)) {
    
    # Full factorial of factors: scenario, setting, year, age, metric
    template_dt = 
      expand_grid(
        scenario  = c("no_vaccine", "vaccine"), 
        setting   = geo[[setting]], 
        year      = o$years, 
        age_group = paste1("age_group", 1 : 7),
        metric    = sort(metrics)) %>%
      # Append random result placeholders...
      mutate(value = runif(n())) %>%
      rename(!!setting := setting) %>%
      as.data.table()
    
    # Write to file
    file = paste0(o$pth$extern, "template_polio_", setting, ".csv")
    fwrite(template_dt, file = file)
  }
}

# ---------------------------------------------------------
# Prepare and simulate DynaMICE measles model
# ---------------------------------------------------------
simulate_dynamice = function() {
  
  # TODO: How important are 'mid_day' and 'coverage_subnat' 
  #       and how should they be derived?
  
  # Return out now if direct simulation not required
  if (!o$simulate_dynamice)
    return()
  
  # Extract path of local DynaMICE repo
  repo_path = repo_exists("dynamice")
  
  # Throw error if repo doesn't exist locally
  if (is.null(repo_path)) {
    
    # Construct error message
    err_msg = paste0(
      "In order to simulate the DynaMICE model, you must: ", 
      "\n  1) Clone the repo '", o$github_dynamice, "'", 
      "\n  2) Have access to a SLURM-queued cluster")
    
    stop(err_msg)
  }
  
  message("\n----- Simulating DynaMICE -----\n")
  
  # ---- Dynamice coverage inputs ----
  
  message("* Model set up")
  
  # Convert EPI50-DynaMICE vaccine references
  dynamice_dict = c(
    MCV1    = "MCV1", 
    MCV2    = "MCV2", 
    Measles = "SIA")
  
  # Load EPI50 coverage details
  coverage_dt = table("coverage_everything") %>%
    inner_join(y  = table("d_v_a_extern"), 
               by = "d_v_a_id") %>%
    filter(vaccine %in% names(dynamice_dict)) %>%
    mutate(vaccine = recode(vaccine, !!!dynamice_dict)) %>%
    select(vaccine, country, year, age, coverage)
  
  # Routine coverage - no age disaggregation required
  routine_dt = coverage_dt %>%
    filter(vaccine != "SIA") %>%
    select(vaccine, country, year, coverage) %>%
    # Summarise over age groups (noting that coverage the same)...
    group_by(vaccine, country, year) %>%
    summarise(coverage = mean(coverage)) %>%
    ungroup() %>%
    as.data.table()
  
  # Non-routine coverage - requires additonal details
  sia_dt = coverage_dt %>%
    filter(vaccine == "SIA") %>%
    # Convert coverage to char to enable pivot...
    mutate(coverage = round(coverage, 8), 
           coverage = as.character(coverage)) %>%
    # Group by all but age to find age bounds per campaign...
    group_by(vaccine, country, year, coverage) %>%
    mutate(age_first = min(age), 
           age_last  = max(age)) %>%
    ungroup() %>%
    # Reduce down to individual campaigns...
    select(-age) %>%
    unique() %>%
    # Convert coverage back to numeric...
    mutate(coverage = as.numeric(coverage)) %>%
    # Append any other additional variables needed...
    mutate(mid_day = 180) %>%             # NOTE: A placeholder assumption
    # Assume trivial subntaional coverage...
    mutate(coverage_subnat = coverage,    # NOTE: A placeholder assumption
           .after = coverage) %>%  
    as.data.table()
  
  # ---- Save input files to DynaMICE repo ----
  
  # Concatenate routine and non-routine coverage data
  data_dt = bind_rows(routine_dt, sia_dt) %>%
    arrange(vaccine, country, year, 
            age_first, age_last)
  
  # Inputs for 'all vaccine' and 'no vaccine' scenarios
  data_list = list(
    mcv1_mcv2_sia = data_dt, 
    nomcv         = data_dt[coverage == 0])
  
  # Iterate through scenarios
  for (scenario in names(data_list)) {
    
    # Construct file path to save to (in DynaMICE repo)
    save_path = file.path(repo_path, "input", "coverage", "coverage")
    save_file = paste0(paste1(save_path, scenario), ".csv")
    
    # Save data as a csv
    fwrite(data_list[[scenario]], file = save_file)
  }
  
  # ---- Other config files ----
  
  # Save full EPI50 country list
  country_dt   = data.table(country = all_countries())
  country_file = file.path(repo_path, "config", "countries.csv")
  
  # Also save associated regions
  region_dt   = table("country")[, .(country, region)]
  region_file = file.path(repo_path, "config", "regions.csv")
  
  # Write to config folder
  fwrite(country_dt, file = country_file)
  fwrite(region_dt,  file = region_file)
  
  # ---- Simulate model ----
  
  # Set working directory to DynaMICE repo
  setwd(repo_path)
  
  # Launch to model
  #
  # NOTE: For full functionality, step should be set to 1 : 3 in DynaMICE repo
  system("sh launch.sh")
  
  # Once we're done, reset working directory to EPI50 repo
  setwd(o$pth$code)
  
  message("\n----- DynaMICE complete -----\n")
  
  # Name of DynaMICE results file that should have been produced
  results_name = "epi50_dynamice_results.rds"
  results_file = file.path(repo_path, "output", results_name)
  
  # Copy results file from DynaMICE repo to EPI50 repo
  invisible(file.copy(results_file, o$pth$extern, overwrite = TRUE))
}

# ---------------------------------------------------------
# TEMP: Read in external polio results
# ---------------------------------------------------------
read_polio = function() {
  
  # Check whether dummy results already exist
  results_exist = file.exists(paste0(o$pth$extern, "epi50_polio_results.rds"))
  
  # Return out if process not necessary
  if (results_exist && o$dummy_polio == FALSE)
    return()
  
  message(" - Reading polio outcomes")
  
  # Load template of regional results
  template_file = "template_polio_region.csv"
  read_dt = fread(paste0(o$pth$extern, template_file)) 
  
  # Save EPI50-formatted polio results
  save_rds(read_dt, "extern", "epi50_polio_results")

}

# ---------------------------------------------------------
# TEMP: Read in Ferrari's measles results
# ---------------------------------------------------------
read_penn_state = function() {
 
  # Check whether results already exist
   results_exist = file.exists(paste0(o$pth$extern, "epi50_penn_state_results.rds"))
  
  # Return out if process not necessary
  if (results_exist && o$dummy_penn_state == FALSE)
   return()
  
  message(" - Reading Penn State measles outcomes")
  
  # Load input files
  no_vaccination_file = "cases_and_deaths_estimates_no_vaccination_2022.csv"
  no_vaccines_dt = fread(paste0(o$pth$extern, no_vaccination_file)) %>%
                     select(country, year, deaths) %>%
                     mutate(scenario = "no_vaccine")
  
  vaccination_file = "cases_and_deaths_estimates_2022.csv"
  vaccines_dt = fread(paste0(o$pth$extern, vaccination_file)) %>%
                 select(country, year, deaths) %>%
                 mutate(scenario = "vaccine")
  
  summary_dt = bind_rows(no_vaccines_dt, vaccines_dt) 
  
  # Save summary of Penn State measles results
  save_rds(summary_dt, "extern", "epi50_penn_state_results.rds")
  
}

# ---------------------------------------------------------
# Format measles modelling results for EPI50 use
# ---------------------------------------------------------
format_measles = function() {
  
  message(" - Appending to measles outcomes")
  
  # Dictionary for converting dose names
  dose_dict = c(
    MCV1_doses = "MCV1",
    MCV2_doses = "MCV2",
    SIA_doses  = "Measles")
  
  # Load template of regional results
  template_file = "template_measles.csv"
  template_dt = fread(paste0(o$pth$extern, template_file)) %>%
    mutate(metric = recode(metric, !!!dose_dict)) %>%
    select(-value)
  
  # Load measles vaccine coverage by vaccine
  fvps_dt = table("coverage_everything") %>%
    inner_join(y  = table("d_v_a_extern"), 
               by = "d_v_a_id") %>%
    filter(disease == "Measles") %>%
    mutate(scenario = "vaccine") %>%
    select(scenario, country, year, age, 
           metric = vaccine, value = fvps)
  
  # All measles models to append to
  measles_models = table("extern_models") %>%
    filter(disease == "Measles") %>%
    pull(model)
  
  # Iterate through measles models
  for (model in measles_models) {
    
    # File names for raw and formatted results
    raw_name   = paste1("epi50",  model, "results")
    table_name = paste1("extern", model, "results")

    # Load raw results, removing any appended FVPs info
    raw_dt = read_rds("extern", raw_name) %>%
      filter(!metric %in% dose_dict)
    
    # Append original FVPs from coverage table
    model_dt = template_dt %>%
      left_join(y  = rbind(raw_dt, fvps_dt),
                by = names(template_dt)) %>%
      replace_na(list(value = 0))
    
    # Save EPI50-formatted results for this model 
    save_table(model_dt, table_name)
  }
}

# ---------------------------------------------------------
# Format polio modelling results for EPI50 use
# ---------------------------------------------------------
format_polio = function() {
  
  message(" - Interpolating polio outcomes")
  
  # TODO: Convert doses into FVPs by dividing through...
  
  # Load raw polio results
  raw_dt = read_rds("extern", "epi50_polio_results")
  
  # ---- Expand age groups in single years ----
  
  # Age groupings as defined in polio results
  age_groups   = sort(unique(raw_dt$age_group))
  age_group_dt = data.table(age_group = age_groups) %>% 
    mutate(age = 2 ^ (seq_along(age_groups) - 1),
           age = pmin(age, max(o$ages)))
  
  # Construct age datatable to expand age bins to single years
  age_dt = data.table(age = o$ages) %>%
    left_join(y  = age_group_dt, 
              by = "age") %>%
    fill(age_group, .direction = "downup") %>%
    group_by(age_group) %>%
    add_count(age_group) %>%
    ungroup() %>%
    as.data.table()
  
  # ---- Crudely expand regions to countries ----
  
  # We'll (very crudely) disaggregate results into countries
  #
  # NOTE: This is simply to have consistent format with other diseases
  setting_dt = table("country") %>%
    select(region, country) %>%
    # Append population size in most recent year...
    mutate(year = max(o$years)) %>%
    left_join(y  = table("wpp_pop"), 
              by = c("country", "year")) %>%
    # Summarise over all ages...
    group_by(region, country) %>%
    summarise(pop = sum(pop)) %>%
    ungroup() %>%
    # Country population share by region...
    group_by(region) %>%
    mutate(pop_share = pop / sum(pop)) %>%
    ungroup() %>%
    select(-pop) %>%
    as.data.table()
  
  # ---- Disaggregate raw results by country and age ----
  
  # Bring it all together to disaggregate raw results...
  polio_dt = raw_dt %>%
    # Expand age groups to all ages...
    full_join(age_dt, by = "age_group", 
              relationship = "many-to-many") %>%
    select(-age_group) %>%
    # Expand regions to countries...
    full_join(setting_dt, by = "region", 
              relationship = "many-to-many") %>%
    select(-region) %>%
    # Divide results through for each country and age ...
    mutate(value = (value * pop_share) / n) %>%
    select(scenario, country, year, age, metric, value) %>%
    arrange(scenario, country, year, age)
  
  # ---- Sanity checks ----
  
  # Function to compute total outcomes
  total_fn = function(dt, name) {
    
    # Total outcomes by scenario and metric
    total_dt = dt %>%
      group_by(scenario, metric) %>%
      summarise(value = sum(value)) %>%
      ungroup() %>%
      rename(!!name := value) %>%
      as.data.table()
    
    return(total_dt)
  }
  
  # Compare raw with formatted model outcomes
  check_dt = total_fn(polio_dt, "clean") %>%
    left_join(y  = total_fn(raw_dt, "raw"), 
              by = c("scenario", "metric")) %>%
    mutate(diff = abs(clean - raw), 
           err  = diff > 1e-6)

  # Throw an error if any differences are identified
  if (any(check_dt$err))
    stop("Error in country or age polio results disaggregation")
  
  # ---- Finally, convery doses to FVPs ----
  
  # Divide doses through to get FVPs
  polio_dt %<>%
    mutate(metric = str_remove(metric, "_doses$")) %>%
    left_join(y  = table("regimen"), 
              by = c("metric" = "vaccine")) %>%
    replace_na(list(schedule = 1)) %>%
    mutate(value = value / schedule) %>%
    select(-schedule)
  
  # Save EPI50-formatted polio results 
  save_table(polio_dt, "extern_polio_results")
}

# ---------------------------------------------------------
# Extract results from all extern models
# ---------------------------------------------------------
extract_extern_results = function() {
  
  message(" - Extracting results from all external models")
  
  # ---- Extract outcomes ----
  
  # Function for extracting model outcomes
  extract_fn = function(model) {
    
    message("   ~ ", model)
    
    # Name of formatted table for this model
    model_table = paste1("extern", model, "results")
    
    # Extract death estimates from vaccine and no vaccine scenarios
    model_dt = table(model_table) %>%
      mutate(model = model, .before = 1)
    
    return(model_dt)
  }
  
  # All extern models and associated disease name
  all_models = table("extern_models") %>%
    pivot_wider(
      names_from  = model, 
      values_from = disease) %>%
    unlist()
  
  # Load historical outcomes from all models
  historical_dt = names(all_models) %>%
    lapply(extract_fn) %>%
    rbindlist() %>%
    lazy_dt() %>%
    # Define d_v_a classification...
    mutate(d_v_a_name = all_models[model]) %>%
    left_join(y  = table("d_v_a"), 
              by = "d_v_a_name") %>%
    # Summarise by d_v_a (mean across all models)...
    group_by(d_v_a_id, scenario, country, year, age, metric) %>%
    summarise(value = mean(value, na.rm = TRUE)) %>%
    ungroup() %>%
    as.data.table()
  
  # ---- Historical deaths and DALYs ----
  
  message("  > Summarising historical estimates")
  
  # Historical deaths in each scenario
  #
  # NOTE: Used for final plotting purposes
  extern_deaths_dt = historical_dt %>%
    filter(metric == "deaths") %>%
    pivot_wider(names_from  = scenario,
                values_from = value) %>%
    as.data.table()
  
  # Save in table cache
  save_table(extern_deaths_dt, "extern_deaths")
  
  # ---- Deaths and DALYs averted ----
  
  message("  > Calculating deaths and DALYs averted")
  
  # Extract deaths and DALYs averted
  extern_averted_dt = historical_dt %>%
    filter(metric %in% c("deaths", "dalys")) %>%
    # Burden in baseline minus burden in vaccine scenario...
    lazy_dt() %>%
    group_by(d_v_a_id, country, year, age, metric) %>%
    mutate(value = value[scenario == "no_vaccine"] - value) %>%
    ungroup() %>%
    mutate(value = pmax(value, 0)) %>%  # TEMP: Experimenting
    # Remove reference to baseline...
    filter(scenario != "no_vaccine") %>%
    select(-scenario) %>%
    # Spread to wide format...
    mutate(metric = paste1(metric, "averted")) %>%
    pivot_wider(names_from  = metric,
                values_from = value) %>%
    as.data.table()
  
  # Save in table cache
  save_table(extern_averted_dt, "extern_estimates")
  
  # TODO: Migrate to plotting.R
  
  # plot_dt = extern_averted_dt %>%
  #   pivot_longer(cols = c(dalys_averted, deaths_averted),
  #                names_to = "metric") %>%
  #   group_by(d_v_a_id, metric, year) %>%
  #   summarise(value = sum(value)) %>%
  #   ungroup() %>%
  #   group_by(d_v_a_id, metric) %>%
  #   mutate(cum_value = cumsum(value)) %>%
  #   ungroup() %>%
  #   format_d_v_a_name() %>%
  #   as.data.table()
  # 
  # g = ggplot(plot_dt) +
  #   aes(x = year, y = cum_value, colour = d_v_a_name) +
  #   geom_line() +
  #   facet_wrap(~metric, scales = "free_y") +
  #   scale_y_continuous(labels = comma)
  
  # ---- Update coverage estimates using model outputs ----
  
  message("  > Extracting vaccine coverage")
  
  # Coverage data prior to appending external pathogen coverage
  base_coverage_dt = table("coverage") %>%
    filter(!d_v_a_id %in% unique(historical_dt$d_v_a_id))
  
  # Extract vaccine coverage from external model outcomes
  extern_coverage_dt = historical_dt %>%
    # Reduce down to non-trival dose estimates...
    filter(scenario == "vaccine",
           metric %in% table("d_v_a_extern")$vaccine, 
           value > 0) %>%
    # Divide doses through by regimen...
    left_join(y  = table("regimen"), 
              by = c("metric" = "vaccine")) %>%
    mutate(fvps = value / schedule) %>%
    # Append cohort and calculate coverage...
    left_join(y  = table("wpp_pop"), 
              by = c("country", "year", "age")) %>%
    rename(cohort = pop) %>%
    mutate(fvps = pmin(fvps, cohort * o$max_coverage),
           coverage = fvps / cohort) %>%
    # Tidy up...
    select(all_names(base_coverage_dt)) %>%
    arrange(d_v_a_id, country, year, age)
  
  # Append external coverage to coverage table
  base_coverage_dt %>%
    rbind(extern_coverage_dt) %>%
    arrange(d_v_a_id, country, year, age) %>%
    save_table("coverage")
}

# ---------------------------------------------------------
# Accredit impact across vaccine types
# ---------------------------------------------------------
split_impact = function() {
  
  browser() # TODO
}

# ---------------------------------------------------------
# Determine if specific repo exists locally
# ---------------------------------------------------------
repo_exists = function(repo) {
  
  # Path for the parent directory of this EPI50 repository
  parent_path = str_remove(o$pth$code, "[a-z,A-Z,0-9]+/$")
  
  # Path to the repo in question
  repo_path = paste0(parent_path, repo)
  
  # If repo exists, return path
  if (dir.exists(repo_path))
    return(repo_path)
  
  # If it doesn't exist, return trivial
  if (!dir.exists(repo_path))
    return(NULL)
}

