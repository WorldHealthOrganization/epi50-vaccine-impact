###########################################################
# COVARIATES
#
# Prepare covariates from various sources for use in regression
# modelling. That is, for geographical imputation, and also
# for inferring drivers of impact.
#
###########################################################

# ---------------------------------------------------------
# Parent function for preparing all covariate data
# ---------------------------------------------------------
prepare_covariates = function() {
  
  message(" > Regression covariates")
  
  # Download Gapminder data (if necessary)
  download_gapminder()
  
  # Format all covariates from the various sources
  c1 = covariates_gapminder()
  c2 = covariates_unicef()
 
  # Concatenate and interpolate
  covariates_dt = c(c1, c2) %>%
    # Expand to temporal scope and interpolate trends...
    interpolate_covariates() %>%
    rbindlist() %>%
    # Tidy up...
    select(metric, country, year, value) %>%
    arrange(metric, country, year)
    
  # Save in tables cache
  save_table(covariates_dt, "regression_covariates")
}

# ---------------------------------------------------------
# Download all Gapminder data from github
# ---------------------------------------------------------
download_gapminder = function() {
  
  # Details of Gapminder data we're interested in
  gapminder_dict = table("gapminder_dict")
  
  # Destination of downloaded files - may or may not exist
  files = paste0(o$pth$gapminder, gapminder_dict$var, ".rds")
  
  # Skip process if all files have already been downloaded
  if (all(file.exists(files)) && !o$force_download_gapminder)
    return()
  
  # ---- Download all files ----
  
  message("  - Downloading data")
  
  # Construct base URL for Gapminder github data files
  gapminder_url = paste0(
    "https://raw.githubusercontent.com/open-numbers/", 
    "ddf--gapminder--systema_globalis/", 
    "master/countries-etc-datapoints/", 
    "ddf--datapoints--[name]--by--geo--time.csv")
  
  # Iterate through metrics of interest
  for (i in seq_row(gapminder_dict)) {
    
    # File name (Gapminder convention)
    data_var = gapminder_dict$file[i]
    
    # Adapt URL for this specific metric
    data_url = str_replace(
      string  = gapminder_url, 
      pattern = "\\[name\\]", 
      replacement = data_var)
    
    # Load the data and briefly format
    data = read_csv(data_url, show_col_types = FALSE) %>%
      rename(value = !!data_var) %>%
      mutate(metric = data_var) %>%
      as.data.table()
    
    # Save rds file locally for easy re-loading
    saveRDS(data, file = files[i])
  }
}

# ---------------------------------------------------------
# Prepare covariates from Gapminder
# ---------------------------------------------------------
covariates_gapminder = function() {
  
  # Gapminder dictionary
  gapminder_dict = table("gapminder_dict")
  
  # Convert to named vector
  vars = setNames(
    gapminder_dict$var, 
    gapminder_dict$file)
  
  # All previously downloaded gapminder files
  files = paste0(o$pth$gapminder, vars, ".rds")
  
  # Load and format gapminder data
  covariates = lapply(files, readRDS) %>%
    rbindlist() %>%
    # Recode variables...
    mutate(metric = recode(metric, !!!vars)) %>%
    # Recode countries
    mutate(country = toupper(geo)) %>%
    filter(country %in% all_countries()) %>%
    # Data ten years prior EPI to capture historical effect...
    mutate(year = as.integer(time)) %>%
    filter(year >= min(o$years) - 10, 
           year <= max(o$years)) %>%
    # Tidy up...
    select(country, year, value, metric) %>%
    split(f = .$metric)
  
  return(covariates)
}

# ---------------------------------------------------------
# Prepare covariates from UNICEF
# ---------------------------------------------------------
covariates_unicef = function() {
  
  # Initiate covariates list
  covariates = list()
  
  # ---- Stunting ----
  
  # Read in UNICEF stunting data
  covariates$stunting = 
    fread(paste0(o$pth$input, "unicef_stunting.csv")) %>%
    # Format column names...
    setnames(names(.), tolower(names(.))) %>% 
    rename(country = "iso code") %>%
    # Reduce down to values of interest...
    filter(country %in% all_countries(), 
           indicator == "Stunting", 
           estimate  == "Point Estimate") %>%
    select(country, matches("^[0-9]+")) %>% 
    # Melt to long format...
    pivot_longer(cols = -country,
                 names_to = "year") %>% 
    # Remove NAs...
    replace_with_na(list(value = "-")) %>% 
    filter(!is.na(value)) %>% 
    # Format values...
    mutate(value  = as.numeric(value), 
           year   = as.integer(year), 
           metric = "stunting") %>%
    as.data.table()

  return(covariates)
}

# ---------------------------------------------------------
# Interpolate timeseries trends for all metrics
# ---------------------------------------------------------
interpolate_covariates = function(covariates) {
  
  message("  - Interpolating timeseries trends")
  
  # Function to perform interpolation
  interpolate_fn = function(data) {
    
    # Expand years to data limit and interpolate trends
    interp_data = data %>%
      # Normalise...
      mutate(value = value / max(value)) %>%
      # Expand to complete temporal scope...
      complete(country, year = min(year) : max(o$years)) %>%
      # Interpolate timeseries trends...
      as_tsibble(index = year, key = country) %>%
      interp_ts_trend() %>%
      # Tidy up...
      mutate(metric = unique(data$metric)) %>%
      as.data.table()
    
    return(interp_data)
  }
  
  # Interpolate metrics in parallel
  if (o$parallel$interp)
    interp_list = mclapply(
      X   = covariates,
      FUN = interpolate_fn, 
      mc.cores = o$n_cores,
      mc.preschedule = FALSE)
  
  # Interpolate metrics consecutively
  if (!o$parallel$interp)
    interp_list = lapply(
      X   = covariates,
      FUN = interpolate_fn)

  return(interp_list)
}

