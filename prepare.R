###########################################################
# PREPARE
#
# Prepare various data sources for use throughout the analysis.
# The idea is that this process needs to be done only once, 
# and that the prepared inputs are streamlined for quick loading.
#
###########################################################

# ---------------------------------------------------------
# Parent function for preparing model inputs from various data sources
# ---------------------------------------------------------
run_prepare = function() {
  
  # Only continue if specified by do_step
  if (!is.element(1, o$do_step)) return()
  
  message("* Preparing input data")
  
  # Convert config yaml files to datatables
  prepare_config_tables()

  # Streamline VIMC impact estimates for quick loading
  prepare_vimc_estimates()

  # Prepare GBD estimates of deaths for non-VIMC pathogens
  prepare_gbd_estimates()

  # Parse vaccine efficacy profile for non-VIMC pathogens
  prepare_vaccine_efficacy()

  # Prepare country income status classification over time
  prepare_income_status()

  # Prepare demography-related estimates from WPP
  prepare_demography()

  # Prepare all covariates for regression modelling
  prepare_covariates()  # See covariates.R

  # Prepare historical vaccine coverage
  prepare_coverage()  # See coverage.R
}

# ---------------------------------------------------------
# Convert config yaml files to datatables 
# ---------------------------------------------------------
prepare_config_tables = function() {
  
  message(" > Config files")
  
  # NOTE: Convert from yaml (/config) to rds (/tables) for fast loading
  
  # List of config yaml files to convert
  config_files = o$pth$config %>%
    list.files(pattern = ".+\\.yaml$") %>%
    str_remove(".yaml$")
  
  # Iterate through these files
  for (file in config_files) {
    
    # Load the yaml file
    yaml_file = paste0("config/", file, ".yaml")
    yaml_data = read_yaml(yaml_file)$table
    
    # Convert to datatable
    config_dt = yaml_data %>%
      lapply(as.data.table) %>%
      rbindlist(fill = TRUE)
    
    # Save in tables cache
    save_table(config_dt, file)
  }
}

# ---------------------------------------------------------
# Streamline VIMC impact estimates for quick loading
# ---------------------------------------------------------
prepare_vimc_estimates = function() {
  
  message(" > VIMC estimates")
  
  # All diseases to load VIMC outcomes for
  vimc_info = table("d_v_a") %>%
    filter(source == "vimc") %>%
    select(disease) %>%
    unique() %>%
    left_join(y  = table("disease_name"), 
              by = "disease")
  
  # Initiate list to store outcomes
  vimc_list = list()
  
  # Iterate through diseases
  for (i in seq_row(vimc_info)) {
    
    # Disease ID and associated full name
    id   = vimc_info[i]$disease
    name = vimc_info[i]$disease_name

    message("  - ", name)
    
    # Load VIMC impact estimates for this disease
    vimc_list[[id]] = read_rds("vimc", id) %>%
      lazy_dt() %>%
      pivot_longer(cols = ends_with("impact"),
                   names_to = "vaccine") %>%
      replace_na(list(value = 0)) %>%
      # Intrepet disease, vaccine, and activity...
      mutate(disease  = tolower(disease),
             vaccine  = str_remove(vaccine, "_impact"),
             activity = ifelse(
               test = vaccine %in% c("routine", "campaign"),
               yes  = vaccine,
               no   = "routine")) %>%
      # Tidy up...
      select(disease, vaccine, activity, country,
             year, age, metric = outcome, value) %>%
      as.data.table()
  }
  
  # Squash results into single datatable
  vimc_dt = rbindlist(vimc_list) %>%
    # Interpret activity...
    mutate(vaccine = ifelse(
      test = vaccine %in% c("routine", "campaign"),
      yes  = disease,
      no   = vaccine)) %>%
    # Deal with rubella special case...
    mutate(is_all   = vaccine == "combined",
           vaccine  = ifelse(is_all, disease, vaccine),
           activity = ifelse(is_all, "all", activity)) %>%
    # Wide format of metrics...
    mutate(metric = paste1(metric, "averted")) %>%
    pivot_wider(names_from = metric) %>%
    # Append d-v-a ID...
    left_join(y  = table("d_v_a"),
              by = c("disease", "vaccine", "activity")) %>%
    # Tidy up...
    select(d_v_a_id, country, year, age,
           deaths_averted, dalys_averted) %>%
    arrange(d_v_a_id, country, year, age) %>%
    as.data.table()
  
  # Save in tables cache
  save_table(vimc_dt, "vimc_estimates")
}

# ---------------------------------------------------------
# Prepare GBD estimates of deaths for non-modelled pathogens
# ---------------------------------------------------------
prepare_gbd_estimates = function() {
  
  message(" > GBD estimates")
  
  # Dictionary of GBD disease names
  gbd_dict = table("gbd_dict") %>%
    rename(name  = gbd_name, 
           value = disease) %>%
    pivot_wider() %>%
    as.list()
  
  # ---- Age mapping ----
  
  # Parse specific age strings
  age_dict = c(
    "<28 days"     = "-1", 
    "28..364 days" = "0..1", 
    "80+ years"    = "80..100")
  
  # Age bins in data before transformation
  age_bins = c(-1, 0, 1, seq(5, 80, by = 5))
  
  # Construct age datatable to expand age bins to single years
  age_dt = data.table(age = c(-1, o$ages)) %>%
    mutate(age_bin = ifelse(age %in% age_bins, age, NA)) %>%
    fill(age_bin, .direction = "down") %>%
    group_by(age_bin) %>%
    add_count(age_bin) %>%
    ungroup() %>%
    as.data.table()
  
  # ---- Load and format data ----
  
  message("  - Loading data")
  
  # Initiate list to store burden results
  burden_list = list()
  
  # Iterate through burden metrics to load
  for (metric in o$metrics) {
    
    # File path to GBD burden file 
    file_name = paste1("gbd19", metric)
    file_path = paste0(o$pth$input, file_name, ".csv")
    
    # Load GBD burden estimates for relevant diseases
    burden_list[[metric]] = fread(file_path) %>%
      # Parse disease and countries...
      mutate(disease = recode(cause, !!!gbd_dict), 
             country = countrycode(
               sourcevar   = location,
               origin      = "country.name", 
               destination = "iso3c")) %>%
      # Retain only what we're interesting in...
      filter(disease %in% table("d_v_a")$disease, 
             country %in% all_countries()) %>%
      # Parse age groups...
      mutate(age = str_replace(age, "-", ".."), 
             age = recode(age, !!!age_dict), 
             age_bin = str_extract(age, "^-*[0-9]+"), 
             age_bin = as.numeric(age_bin)) %>%
      # Reduce down to variables of interest...
      select(disease, country, year, age_bin, value = val) %>%
      arrange(disease, country, year, age_bin) %>%
      mutate(metric = !!metric)
  }
  
  # Squash then split by disease-metric-age
  gbd_list = burden_list %>%
    rbindlist() %>%
    split(f = list(
      .$disease, 
      .$metric, 
      .$age_bin))
  
  # Remove any trivial splits
  gbd_list[lapply(gbd_list, nrow) == 0] = NULL
  
  # ---- Extrapolate recent years ----
  
  message("  - Extrapolating trends")
  
  # Function for extrapolating trends for post-2019 period
  extrap_fn = function(data) {
    
    # Constant method...
    if (o$gbd_extrap == "constant") {
      
      # Expand temporal scope and extrapolate
      extrap_data = data %>%
        # Expand out to recent years without data...
        complete(country, disease, age_bin, metric,
                 year = min(year) : max(o$years)) %>%
        # Extrapolate most recent value...
        fill(value, .direction = "down") %>%
        # filter(!is.na(value)) %>%
        select(all_names(data)) %>%
        as.data.table()
    }
    
    # Using time series method...
    if (o$gbd_extrap == "trend") {
      
      # Identifiers of this split
      id_dt = data %>%
        select(disease, metric, age_bin) %>%
        unique()
      
      # Expand temporal scope and extrapolate
      extrap_data = data %>%
        # Expand out to recent years without data...
        complete(country, year = min(year) : max(o$years)) %>%
        # Extrapolate recent trends...
        as_tsibble(index = year,
                   key   = country) %>%
        interp_ts_trend() %>%
        as.data.table() %>%
        # Tidy up...
        cbind(id_dt) %>%
        select(all_names(data))
    }
    
    return(extrap_data)
  }
  
  # Interpolate metrics in parallel
  if (o$parallel$interp)
    extrap_list = mclapply(
      X   = gbd_list,
      FUN = extrap_fn, 
      mc.cores = o$n_cores,
      mc.preschedule = FALSE)
  
  # Interpolate metrics consecutively
  if (!o$parallel$interp)
    extrap_list = lapply(
      X   = gbd_list,
      FUN = extrap_fn)
  
  # Squash everything back together
  extrap_dt = extrap_list %>%
    rbindlist() %>%
    arrange(metric, disease, country, year, age_bin)
  
  # ---- Expand to singel age bins ----
  
  message("  - Expanding age groups")
  
  # Expand to all ages
  gbd_dt = extrap_dt %>%
    full_join(y  = age_dt, 
              by = "age_bin", 
              relationship = "many-to-many") %>%
    mutate(value = value / n) %>%
    select(disease, country, year, age, value, metric)
  
  # Squash into single datatable
  gbd_dt %>%
    mutate(metric = paste1(metric, "disease")) %>%
    # Pivot metrics to wide format...
    pivot_wider(names_from = metric) %>%
    replace_na(list(
      deaths_disease = 0, 
      dalys_disease  = 0)) %>%
    arrange(disease, country, year, age) %>%
    as.data.table() %>%
    # Save in tables cache
    save_table("gbd_estimates")
  
  # Plot GBD death estimates by age
  plot_gbd_estimates()
}

# ---------------------------------------------------------
# Parse vaccine efficacy profile for non-modelled pathogens
# ---------------------------------------------------------
prepare_vaccine_efficacy = function() {
  
  message(" > Vaccine efficacy")
  
  # ---- Optimisation functions ----
  
  # Function to determine optimal immunity profiles parameters
  optimisation_fn = function(vaccine) {
    
    # Efficacy details (incl data) for this vaccine
    efficacy_info = table("vaccine_efficacy") %>%
      filter(vaccine == !!vaccine)
    
    # Extract the data points (efficacy, year)
    data = efficacy_info %>%
      pull(data) %>%
      unlist() %>%
      matrix(nrow = 2) %>%
      t()
    
    # Extract user-defined functional form for vaccine efficacy
    fn = eval_str(unique(efficacy_info$fn))
    
    # Number of function input arguments (without default values)
    #
    # NOTE: These are the set of values to be optimised
    n_args = sum(!unlist(lapply(formals(fn), is.numeric)))
    
    # Repeat call to optimisation algorithm multiple times
    optim_list = lapply(
      X      = 1 : o$n_optim, 
      FUN    = asd_fn, 
      data   = data, 
      fn     = fn, 
      n_args = n_args)
    
    # Extract parameters from best fitting result
    optim_pars = optim_list %>%
      rbindlist() %>%
      pivot_wider() %>%
      slice_min(y, n = 1, with_ties = FALSE) %>%
      select(-id, -y) %>%
      unlist() %>%
      unname() %>%
      as.list()
    
    # Convert into datatable
    pars_dt = optim_pars %>%
      as_named_dt(letters[seq_along(optim_pars)]) %>%
      mutate(vaccine = vaccine) %>%
      pivot_longer(cols = -vaccine, 
                   names_to = "var") %>%
      as.data.table()
    
    # Evaluate function using optimal parameters
    profile = do.call(fn, optim_pars)
    
    # Form profile into a datatable
    profile_dt = data.table(
      vaccine = vaccine,
      var     = t,
      value   = profile)
    
    # Bind optimal parameters and optimal profile for single output
    output_dt = rbind(pars_dt, profile_dt)
    
    return(output_dt)
  }
  
  # Objective algorithm
  asd_fn = function(i, data, fn, n_args) {
    
    # Reset random number generator
    set.seed(i)
    
    # Fit all required parameters to the data available
    optim = asd(
      fn = obj_fn,
      x0 = runif(n_args),
      lb = 1e-6, 
      ub = 1e6,
      iters = 1e3, 
      args  = list(
        data   = data, 
        fn     = fn, 
        n_args = n_args))
    
    # Convert result to datatable
    result = optim[qc(x, y)] %>%
      unlist() %>%
      enframe() %>%
      mutate(id = i) %>%
      as.data.table()
    
    return(result)
  }
  
  # Objective function to minimise
  obj_fn = function(x, args) {
    
    # Evalulate immunity function
    y = do.call(args$fn, as.list(x))
    
    # Data points we want to hit
    data_x = args$data[, 2] + 1
    data_y = args$data[, 1]
    
    # Calculate sum of squared error
    obj_val = sum((y[data_x] - data_y) ^ 2)
    
    return(list(y = obj_val))
  }
  
  # ---- Perform optimisation for each vaccine ----
  
  # Points at which to evaluate efficacy functions
  t = seq_along(o$years) - 1  # Immunity in the years following vaccination
  
  # Vaccines we want efficacy profiles for (all static modelled vaccines)
  vaccines = table("d_v_a")[source == "static", vaccine]
  
  # Apply optimisation to determine optimal immunity parameters
  optim_results = vaccines %>%
    lapply(optimisation_fn) %>%
    rbindlist()
  
  # Extract optimal profiles
  profile_dt = optim_results %>%
    filter(grepl("^[0-9]+$", var)) %>%
    mutate(time = as.integer(var)) %>%
    left_join(y  = table("d_v_a"), 
              by = "vaccine") %>%
    select(disease, vaccine, time, profile = value)
  
  # Extract optimal parameters
  pars_dt = optim_results %>%
    filter(grepl("^[a-z]+$", var)) %>%
    rename(parameter = var)
  
  # Save both in tables cache
  save_table(profile_dt, "vaccine_efficacy_profiles")
  save_table(pars_dt,    "vaccine_efficacy_parameters")
  
  # Plot these profiles
  plot_vaccine_efficacy()
}

# ---------------------------------------------------------
# Prepare country income status classification over time
# ---------------------------------------------------------
prepare_income_status = function() {
  
  message(" > Income status")
  
  # Path to data file
  #
  # SOURCE: https://datacatalogfiles.worldbank.org/ddh-published/0037712/
  #         DR0090755/CLASS.xlsx?versionId=2023-11-16T18:35:30.5758473Z
  #
  # Alternatively, download 'Historical classification by income' Excel file from: 
  # datacatalog.worldbank.org/search/dataset/0037712/World-Development-Indicators
  file = paste0(o$pth$input, "worldbank_income_status.csv")
  
  # Full country-year combination
  full_dt = expand_grid(
    country = all_countries(), 
    year    = o$years) %>%
    as.data.table()
  
  # Load and format country income status over time
  income_dt = fread(file, header = TRUE) %>%
    # Countries of interest...
    filter(country %in% all_countries()) %>%
    select(-country_name) %>%
    # Convert to tidy format...
    pivot_longer(cols = -country, 
                 names_to  = "year", 
                 values_to = "income") %>%
    # Country with all full country-year combo...
    mutate(year = as.integer(year)) %>%
    full_join(y  = full_dt, 
              by = c("country", "year")) %>%
    arrange(country, year) %>%
    # Fill missing data with pro/preceding value...
    mutate(income = ifelse(income == "", NA, income)) %>%
    group_by(country) %>%
    fill(income, .direction = "downup") %>%
    ungroup() %>%
    # Niue and Cook Islands missing, both are HIC...
    replace_na(list(income = "H")) %>%
    mutate(income = paste0(tolower(income), "ic")) %>%
    as.data.table()
  
  # Save in tables cache
  save_table(income_dt, "income_status")
}

# ---------------------------------------------------------
# Prepare demography-related estimates from WPP
# ---------------------------------------------------------
prepare_demography = function() {
  
  message(" > Demography data")
  
  # Function to apply element-wise scaler to data
  scaler_fn = function(m) {
    
    # Population scaling: by country, year, and age
    if (grepl("pop", m$scale))
      scaler_dt = setnames(
        x   = table("wpp_pop"), 
        old = "pop", 
        new = "scaler")
    
    # Numeric values: simple repitition
    if (grepl("^[0-9,\\.]+$", m$scale))
      scaler_dt = expand_grid(
        country = all_countries(), 
        year    = o$years,
        age     = if (m$age) o$ages else NA, 
        scaler  = as.numeric(m$scale)) %>%
        as.data.table()
    
    return(scaler_dt)
  }
  
  # Details of WPP metrics to load
  wpp_metrics = table("wpp_dict") %>%
    mutate(metric = fct_inorder(metric)) %>%
    split(.$metric)
  
  # Iterate through metrics to load
  for (metric in names(wpp_metrics)) {
    m = as.list(wpp_metrics[[metric]])
    
    message("  - ", metric)
    
    # Past and future in separate data sets
    if (m$proj == TRUE) {
      
      # Age-disaggregation specified in data set file name
      age = ifelse(m$age, "Age", "")
      
      # Names of WPP2022 data files to load
      past   = paste0(m$file,         age, o$pop_bin, "dt")
      future = paste0(m$file, "proj", age, o$pop_bin, "dt") 
      
      # Load pop data from WPP github package
      data_list = data_package(past, future, package = "wpp2022")
    }
    
    # Past and future combined into single data set
    if (m$proj == FALSE) {
      
      # Name of WPP2022 data file - history and projection in one
      all_time = paste0(m$file, o$pop_bin, "dt")
      
      # Load data from WPP github package
      data_list = data_package(all_time, package = "wpp2022")
    }
    
    # Combine (extended) past and future data
    data_dt = rbindlist(data_list, fill = TRUE) %>%
      {if (!m$age) mutate(., age = NA) else .} %>%
      # Select countries of interest...
      inner_join(y  = table("country"),  
                 by = "country_code") %>%
      select(country, year, age, value = !!m$var) %>%
      # Shift year by one (see github.com/PPgp/wpp2022 for details)...
      mutate(year = as.integer(year) + 1) %>%
      filter(year %in% o$years) %>%
      # Scale metrics...
      left_join(y = scaler_fn(m), 
                by = c("country", "year", "age")) %>%
      mutate(value = value * scaler) %>%
      select(-scaler) %>%
      # Tidy up...
      rename(!!metric := value) %>%
      arrange(country, year, age)
    
    # Save in tables cache
    save_table(data_dt, paste1("wpp", metric))
  }
}

# ---------------------------------------------------------
# Simple wrapper to load all countries
# ---------------------------------------------------------
all_countries = function(as_dt = FALSE) {
  
  # Pull all countries defined in config file
  countries = table("country")$country
  
  # Convert to simple datatable if desired
  if (as_dt == TRUE)
    countries = data.table(country = countries)
  
  return(countries)
}

# ---------------------------------------------------------
# Simple wrapper to load all regions
# ---------------------------------------------------------
all_regions = function() {
  
  # Pull all regions defined in config file
  regions = table("region_dict")$region
  
  return(regions)
}

# ---------------------------------------------------------
# Save table in cache directory for quick loading
# ---------------------------------------------------------
save_table = function(x, table) {
  
  # Save table in tables cache directory
  save_rds(x, "tables", table, "table")
}

# ---------------------------------------------------------
# Load and return cached datatable
# ---------------------------------------------------------
table = function(table) {
  
  # TODO: If 'table' contains 'extern' report step 2 back to user 
  
  # Construct file path
  file = paste0(o$pth$tables, table, "_table.rds")
  
  # Throw an error if this file doesn't exist
  if (!file.exists(file))
    stop("Table ", table, " has not been cached - have you run step 1?")
  
  # Load rds file
  y = read_rds(file)
  
  return(y)
}

