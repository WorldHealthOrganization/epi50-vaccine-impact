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
 # prepare_config_tables()

  # Streamline VIMC impact estimates for quick loading
  #prepare_vimc_estimates()

  # Parse vaccine efficacy profile for non-VIMC pathogens
  #prepare_vaccine_efficacy()

  # Prepare GBD estimates of deaths for non-VIMC pathogens
  #prepare_gbd_estimates()

  # Prepare GBD covariates for extrapolating to non-VIMC countries
  #prepare_gbd_covariates()

  # Prepare UNICEF covariates for imputing non-VIMC countries
  prepare_unicef()
  
  # Prepare World Bank covariates for imputing non-VIMC countries
  prepare_world_bank()

  # Prepare Gapminder covariates for imputing non-VIMC countries
  prepare_gapminder()

  # Prepare country income status classification over time
  prepare_income_status()

  # Prepare demography-related estimates from WPP
  prepare_demography()

  # Prepare age at birth by country and year
  prepare_birth_age()
  
  # Prepare historical vaccine coverage
  prepare_coverage()  # See coverage.R
}

# ---------------------------------------------------------
# Convert config yaml files to datatables 
# ---------------------------------------------------------
prepare_config_tables = function() {
  
  message(" - Config files")
  
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
  
  message(" - VIMC estimates")
  
  # TODO: In raw form, we could instead use vimc_estimates.csv
  
  # Prepare VIMC vaccine impact estimates
  read_rds("input", "vimc_estimates") %>%
    # Years and countries of interest...
    filter(country %in% all_countries(), 
           year    %in% o$years) %>%
    # Disease, vaccines, and activities of interest...
    left_join(y  = table("d_v_a"), 
              by = c("disease", "vaccine", "activity")) %>%
    filter(!is.na(d_v_a_id), 
           source == "vimc") %>%
    # Tidy up...
    select(d_v_a_id, country, year, age, deaths_averted) %>%
    arrange(d_v_a_id, country, year, age) %>%
    save_table("vimc_estimates")
  
  # Simply store VIMC in it's current form
  read_rds("input", "vimc_uncertainty") %>%
    save_table("vimc_uncertainty")
}

# ---------------------------------------------------------
# Parse vaccine efficacy profile for non-modelled pathogens
# ---------------------------------------------------------
prepare_vaccine_efficacy = function() {
  
  message(" - Vaccine efficacy")
  
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
    
    # Fit all required parameters to the data available
    optim = asd(
      fn   = obj_fn,
      x0   = runif(n_args),
      args = list(
        data   = data, 
        fn     = fn, 
        n_args = n_args),
      lb   = 1e-6, 
      ub   = 1e6,
      max_iters = 1e3)
    
    # Evaluate function using optimal parameters
    profile = do.call(fn, as.list(optim$x))
    
    # Form profile into a datatable
    profile_dt = data.table(
      vaccine = vaccine,
      time    = t,
      profile = profile)
    
    return(profile_dt)
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
  lapply(vaccines, optimisation_fn) %>%
    rbindlist() %>%
    left_join(y  = table("d_v_a"), 
              by = "vaccine") %>%
    select(disease, vaccine, time, profile) %>%
    save_table("vaccine_efficacy_profiles")
  
  # Plot these profiles
  plot_vaccine_efficacy()
}

# ---------------------------------------------------------
# Prepare GBD estimates of deaths for non-modelled pathogens
# ---------------------------------------------------------
prepare_gbd_estimates = function() {
  
  message(" - GBD estimates")
  
  # Parse specific age strings
  age_dict = c(
    "<28 days"     = "-1", 
    "28..364 days" = "0..1", 
    "80+ years"    = "80..95")
  
  # Dictionary of GBD disease names
  gbd_dict = table("gbd_dict") %>%
    rename(name  = gbd_name, 
           value = disease) %>%
    pivot_wider() %>%
    as.list()
  
  # Load GBD estimates of deaths for relevant diseases
  deaths_dt = fread(paste0(o$pth$input, "gbd19_deaths.csv")) %>%
    filter(year %in% o$gbd_estimate_years) %>%
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
    # Tidy up...
    select(disease, country, year, age_bin, value = val) %>%
    arrange(disease, country, year, age_bin)
  
  warning("Improvements to GBD death extrapolations required...")
  
  # TEMP: Until we do a proper future projection of GBD estimates
  deaths_dt = 
    expand_grid(
      disease = unique(deaths_dt$disease),
      country = all_countries(), 
      age_bin = unique(deaths_dt$age_bin), 
      year    = o$years) %>%
    left_join(y  = deaths_dt, 
              by = names(.)) %>%
    arrange(disease, country, age_bin) %>%
    group_by(disease, country, age_bin) %>%
    fill(value, .direction = "down") %>%
    ungroup() %>%
    filter(!is.na(value)) %>%
    as.data.table()
  
  # Age bins in data before and after transformation
  age_bins = sort(unique(deaths_dt$age_bin))
  
  # Construct age datatable to expand age bins to single years
  age_dt = data.table(age = c(-1, o$ages)) %>%
    mutate(age_bin = ifelse(age %in% age_bins, age, NA)) %>%
    fill(age_bin, .direction = "down") %>%
    group_by(age_bin) %>%
    add_count(age_bin) %>%
    ungroup() %>%
    as.data.table()
  
  # Expand to all ages
  deaths_dt %>%
    full_join(y  = age_dt, 
              by = "age_bin", 
              relationship = "many-to-many") %>%
    mutate(deaths_disease = value / n) %>%
    select(disease, country, year, age, deaths_disease) %>%
    arrange(disease, country, year, age) %>%
    save_table("gbd_estimates")
  
  # Plot GBD death estimates by age
  plot_gbd_estimates()
}

# ---------------------------------------------------------
# Prepare GBD covariates for extrapolating to non-modelled countries
# ---------------------------------------------------------
prepare_gbd_covariates = function() {
  
  # TODO: We need to project these estimates to avoid losing impact
  #       estimates in geo-imputation model
  
  message(" - GBD covariates")
  
  # Prepare GBD 2019 HAQI for use as a covariate
  haqi_dt = fread(paste0(o$pth$input, "gbd19_haqi.csv")) %>%
    # Countries of interest...
    rename(country_name = location_name) %>%
    inner_join(y  = table("country"),
               by = "country_name") %>%
    # Tidy up...
    select(country, year = year_id, haqi = val) %>%
    mutate(haqi = haqi / 100) %>%
    arrange(country, year)
  
  # Function to load each SDI file
  load_sdi = function(name) {
    
    # Full file name
    file = paste0(o$pth$input, name, ".csv")
    
    # Load and convert to long form to allow binding
    sdi_dt = fread(file, header = TRUE) %>%
      rename(gbd_alt_name = Location) %>%
      pivot_longer(cols = -gbd_alt_name,
                   names_to  = "year",
                   values_to = "sdi") %>%
      as.data.table()
    
    return(sdi_dt)
  }
  
  # GBD country names
  gbd_name_dt = table("country") %>%
    mutate(gbd_alt_name = ifelse(
      test = is.na(gbd_alt_name), 
      yes  = country_name, 
      no   = gbd_alt_name)) %>%
    select(country, gbd_alt_name)
  
  # Prepare GBD 2019 SDI for use as a covariate
  sdi_dt = rbind(load_sdi("gbd19_sdi_1970"), 
                 load_sdi("gbd19_sdi_1990")) %>%
    # Countries of interest...
    inner_join(y  = gbd_name_dt,
               by = "gbd_alt_name") %>%
    select(country, year, sdi) %>%
    # Years of interest...
    mutate(year = as.integer(year)) %>%
    filter(year %in% o$years) %>%
    # Tidy up...
    arrange(country, year) %>%
    as.data.table()
  
  # Join metrics into single datatable
  gbd_covariates = 
    full_join(x  = sdi_dt, 
              y  = haqi_dt, 
              by = c("country", "year")) %>%
    arrange(country, year)
  
  warning("Improvements to GBD covariate extrapolations required...")
  
  # TEMP: Until we do a proper projection of GBD covariates
  #
  # NOTE: We only really need a forward projection here
  gbd_covariates = 
    expand_grid(country = all_countries(), 
                year    = o$years) %>%
    left_join(y  = gbd_covariates, 
              by = c("country", "year")) %>%
    group_by(country) %>%
    fill(sdi, haqi, .direction = "downup") %>%
    ungroup() %>%
    as.data.table()
  
  # Save in tables cache
  save_table(gbd_covariates, "gbd_covariates")
}

# ---------------------------------------------------------
# Prepare UNICEF stunting and maternal mortality estimates
# ---------------------------------------------------------
prepare_unicef = function() {
  
  message(" - UNICEF")
  
  # Read in WHO region data
  WHO_regions_dt = fread(paste0(o$pth$input, "WHO_country_codes.csv")) 
  
  # Read in UNICEF stunting data
  stunting_dt = fread(paste0(o$pth$input, "JME_Country_Estimates_May_2023.csv")) %>%
    filter(Indicator == "Stunting" &
             Estimate == "Point Estimate") %>%
    select(-c("Country and areas", "Note", "Indicator", "Measure", "Estimate")) %>% 
    pivot_longer(cols = -'ISO code',
                 names_to = "year",
                 values_to = "stunting") %>% 
    rename(country_code = 'ISO code') %>%
    full_join(WHO_regions_dt, by="country_code", relationship = "many-to-many") %>%
    select(-c(country, region_no, region_long)) %>%
    rename(country = country_code) %>%
    mutate(year = as.integer(year)) %>%
    as.data.table()
  
  
  maternal_mortality_dt = fread(paste0(o$pth$input, "unicef_maternal_mortality.csv"), header = TRUE) %>%
    select(-Country) %>%
    pivot_longer(cols = -'ISO Code',
                 names_to = "year",
                 values_to = "maternal_mortality") %>%
    rename(country_code = 'ISO Code') %>%
    full_join(WHO_regions_dt, by="country_code", relationship = "many-to-many") %>%
    select(-c(country, region_no, region_long)) %>%
    rename(country = country_code) %>%
    mutate(year = as.integer(year)) %>%
    filter(!is.na(year)) %>% 
    unique() %>% # remove duplicates for COD and TZN
    complete(country, year = 2000:2024, 
             fill = list(maternal_mortality = NA)) %>%
    as_tsibble(index = year,
               key = country) 
  
  # Interpolate missing values
  maternal_mortality_dt = maternal_mortality_dt %>%
    model(lm = TSLM(log(maternal_mortality) ~ trend())) %>%
    interpolate(maternal_mortality_dt)
  
  
   #browser()
  # Create table of UNICEF covariates
  unicef_dt = stunting_dt  %>%
    full_join(maternal_mortality_dt, by=c("country", "year")) %>%
    filter(!is.na(year))
  
  
  # Save in tables cache
  save_table(unicef_dt, "unicef_covariates")
}

# ---------------------------------------------------------
# Prepare World Bank literacy estimates
# ---------------------------------------------------------
prepare_world_bank = function() {
  
  message(" - World Bank")
 
  # Read in WHO region data
  WHO_regions_dt = fread(paste0(o$pth$input, "WHO_country_codes.csv")) 
  
  # Read in World Bank literacy data
  male_adult_literacy_dt = fread(paste0(o$pth$input, "literacy_male.csv"), header= TRUE) %>% 
    select(-c('Country Name', 'Indicator Name', 'Indicator Code')) %>% 
    pivot_longer(cols = -'Country Code',
                 names_to = "year",
                 values_to = "male_adult_literacy") %>% 
    rename(country_code = 'Country Code') %>%
    full_join(WHO_regions_dt, by="country_code", relationship = "many-to-many") %>%
    select(-c(country, region_no, region_long)) %>%
    rename(country = country_code) %>%
    filter(year != 'V68' & !is.na(year)) %>%
    mutate(year = as.integer(year)) %>%
    unique() %>% # remove duplicates for COD and TZN
    complete(country, year = 1974:2024, 
             fill = list(male_adult_literacy = NA)) %>%
    as_tsibble(index = year,
               key = country) 
  
  # Interpolate missing values
  male_adult_literacy_dt = male_adult_literacy_dt %>%
    # Omit countries with no data
    group_by(country) %>%
    filter(!all(is.na(male_adult_literacy))) %>% 
    ungroup() %>%
    # Fit time series linear regression model
    model(lm = TSLM(log(male_adult_literacy) ~ trend())) %>%
    interpolate(male_adult_literacy_dt)
  
  # Read in World Bank literacy data
  female_adult_literacy_dt = fread(paste0(o$pth$input, "literacy_female.csv"), header= TRUE) %>% 
    select(-c('Country Name', 'Indicator Name', 'Indicator Code')) %>% 
    pivot_longer(cols = -'Country Code',
                 names_to = "year",
                 values_to = "female_adult_literacy") %>% 
    rename(country_code = 'Country Code') %>%
    full_join(WHO_regions_dt, by="country_code", relationship = "many-to-many") %>%
    select(-c(country, region_no, region_long)) %>%
    rename(country = country_code) %>%
    filter(year != 'V68' & !is.na(year)) %>%
    mutate(year = as.integer(year)) %>%
    unique() %>% # remove duplicates for COD and TZN
    complete(country, year = 1974:2024, 
             fill = list(female_adult_literacy = NA)) %>%
    as_tsibble(index = year,
               key = country) 
  
  # Interpolate missing values
  female_adult_literacy_dt = female_adult_literacy_dt %>%
    # Omit countries with no data
    group_by(country) %>%
    filter(!all(is.na(female_adult_literacy))) %>% 
    ungroup() %>%
    # Fit time series linear regression model
    model(lm = TSLM(log(female_adult_literacy) ~ trend())) %>%
    interpolate(female_adult_literacy_dt)
  
  # Read in World Bank literacy data
  male_youth_literacy_dt = fread(paste0(o$pth$input, "literacy_male_15_24.csv"), header= TRUE) %>% 
    select(-c('Country Name', 'Indicator Name', 'Indicator Code')) %>% 
    pivot_longer(cols = -'Country Code',
                 names_to = "year",
                 values_to = "male_youth_literacy") %>% 
    rename(country_code = 'Country Code') %>%
    full_join(WHO_regions_dt, by="country_code", relationship = "many-to-many") %>%
    select(-c(country, region_no, region_long)) %>%
    rename(country = country_code) %>%
    filter(year != 'V68' & !is.na(year)) %>%
    mutate(year = as.integer(year)) %>%
    unique() %>% # remove duplicates for COD and TZN
    complete(country, year = 1974:2024, 
             fill = list(male_youth_literacy = NA)) %>%
    as_tsibble(index = year,
               key = country) 
  
  # Interpolate missing values
  male_youth_literacy_dt = male_youth_literacy_dt %>%
    # Omit countries with no data
    group_by(country) %>%
    filter(!all(is.na(male_youth_literacy))) %>% 
    ungroup() %>%
    # Fit time series linear regression model
    model(lm = TSLM(log(male_youth_literacy) ~ trend())) %>%
    interpolate(male_youth_literacy_dt)
  
  # Read in World Bank literacy data
  female_youth_literacy_dt = fread(paste0(o$pth$input, "literacy_female_15_24.csv"), header= TRUE) %>% 
    select(-c('Country Name', 'Indicator Name', 'Indicator Code')) %>% 
    pivot_longer(cols = -'Country Code',
                 names_to = "year",
                 values_to = "female_youth_literacy") %>% 
    rename(country_code = 'Country Code') %>%
    full_join(WHO_regions_dt, by="country_code", relationship = "many-to-many") %>%
    select(-c(country, region_no, region_long)) %>%
    rename(country = country_code) %>%
    filter(year != 'V68' & !is.na(year)) %>%
    mutate(year = as.integer(year)) %>%
    unique() %>% # remove duplicates for COD and TZN
    complete(country, year = 1974:2024, 
             fill = list(female_youth_literacy = NA)) %>%
    as_tsibble(index = year,
               key = country) 
  
  # Interpolate missing values
  female_youth_literacy_dt = female_youth_literacy_dt %>%
    # Omit countries with no data
    group_by(country) %>%
    filter(!all(is.na(female_youth_literacy))) %>% 
    ungroup() %>%
    # Fit time series linear regression model
    model(lm = TSLM(log(female_youth_literacy) ~ trend())) %>%
    interpolate(female_youth_literacy_dt)
  

  # Create table of World Bank covariates
  world_bank_dt = male_adult_literacy_dt  %>%
    full_join(male_youth_literacy_dt, by=c("country", "year")) %>%
    full_join(female_adult_literacy_dt, by=c("country", "year")) %>%
    full_join(female_youth_literacy_dt, by=c("country", "year"))
  
  # Save in tables cache
  save_table(world_bank_dt, "world_bank_covariates")
}

# ---------------------------------------------------------
# Prepare Gapminder covariates for extrapolating to non-modelled countries
# ---------------------------------------------------------
prepare_gapminder = function() {
  
  message(" - Gapminder covariates")
  
  # TODO: We can remove the region references, handled in table("countries")
  
  # Read in WHO region data
  WHO_regions_dt = fread(paste0(o$pth$input, "WHO_country_codes.csv")) 
  
  # Prepare Gapminder data for use as predictors
  # Gini coefficient
  gini_dt = fread(paste0(o$pth$input, "ddf--datapoints--gapminder_gini--by--geo--time.csv")) %>%
    rename(gini = gapminder_gini)

  # GPD per capita US$ inflation-adjusted
  gdp_dt = fread(paste0(o$pth$input, "ddf--datapoints--gdppercapita_us_inflation_adjusted--by--geo--time.csv")) %>%
    rename(gdp = gdppercapita_us_inflation_adjusted)
  
  # Doctors per 1000 population
  doctors_per_1000_dt = fread(paste0(o$pth$input, "ddf--datapoints--medical_doctors_per_1000_people--by--geo--time.csv")) %>%
    rename(doctors_per_1000 = medical_doctors_per_1000_people) %>%
    complete(geo, time = 1970:2024, 
             fill = list(doctors_per_1000 = NA)) %>%
    as_tsibble(index = time, 
               key = geo) 
  
  # Interpolate missing values
  doctors_per_1000_dt = doctors_per_1000_dt %>%
    model(lm = TSLM(log(doctors_per_1000) ~ trend())) %>%
    interpolate(doctors_per_1000_dt)
  
  # Population aged 0 to 14
  pop_0_to_14_dt = fread(paste0(o$pth$input, "ddf--datapoints--population_aged_0_14_years_both_sexes_percent--by--geo--time.csv")) %>%
    rename(pop_0to14 = population_aged_0_14_years_both_sexes_percent)
  
  # Population density
  pop_density_dt = fread(paste0(o$pth$input, "ddf--datapoints--population_density_per_square_km--by--geo--time.csv")) %>%
    rename(pop_density = population_density_per_square_km)
  
  # Urban population (%)
  urban_dt = fread(paste0(o$pth$input, "ddf--datapoints--urban_population_percent_of_total--by--geo--time.csv")) %>%
    rename(urban_percent = urban_population_percent_of_total)
  
  # HIV mortality
  hiv_mortality_dt = fread(paste0(o$pth$input, "ddf--datapoints--annual_hiv_deaths_number_all_ages--by--geo--time.csv")) %>%
    rename(hiv_mortality = annual_hiv_deaths_number_all_ages)
  
   # Health spending ($)
  health_spending_dt = fread(paste0(o$pth$input, "ddf--datapoints--total_health_spending_per_person_us--by--geo--time.csv")) %>%
    rename(health_spending = total_health_spending_per_person_us) %>%
    complete(geo, time = 1970:2024, 
             fill = list(health_spending = NA)) %>%
    as_tsibble(index = time, 
               key = geo) 
  
  # Interpolate missing values
  health_spending_dt = health_spending_dt %>%
    model(lm = TSLM(log(health_spending) ~ trend())) %>%
    interpolate(health_spending_dt)
  
  # Private health spending (%)
  private_health_spending_dt = fread(paste0(o$pth$input, "ddf--datapoints--private_share_of_total_health_spending_percent--by--geo--time.csv")) %>%
    rename(private_health = private_share_of_total_health_spending_percent) %>%
    complete(geo, time = 1970:2024, 
             fill = list(private_health = NA)) %>%
    as_tsibble(index = time, 
               key = geo) 
  
  # Interpolate missing values
  private_health_spending_dt = private_health_spending_dt %>%
    model(lm = TSLM(log(private_health) ~ trend())) %>%
    interpolate(private_health_spending_dt)
  
  # At least basic sanitation (%)
  sanitation_dt = fread(paste0(o$pth$input, "ddf--datapoints--at_least_basic_sanitation_overall_access_percent--by--geo--time.csv")) %>%
    rename(basic_sanitation = at_least_basic_sanitation_overall_access_percent) %>%
    complete(geo, time = 1970:2024, 
             fill = list(basic_sanitation = NA)) %>%
    as_tsibble(index = time, 
               key = geo) 
  
  # Interpolate missing values
  sanitation_dt = sanitation_dt %>%
    model(lm = TSLM(log(basic_sanitation) ~ trend())) %>%
    interpolate(sanitation_dt)
  
  
  # At least basic water source (%)
  water_dt = fread(paste0(o$pth$input, "ddf--datapoints--at_least_basic_water_source_overall_access_percent--by--geo--time.csv")) %>%
    rename(basic_water = at_least_basic_water_source_overall_access_percent) %>%
    complete(geo, time = 1970:2024, 
             fill = list(basic_water = NA)) %>%
    as_tsibble(index = time, 
               key = geo) 
  
  # Interpolate missing values
  water_dt = water_dt %>%
    model(lm = TSLM(log(basic_water) ~ trend())) %>%
    interpolate(water_dt)
  
  
  # Human development index (life expectancy, education, per-person income)
  hdi_dt = fread(paste0(o$pth$input, "ddf--datapoints--hdi_human_development_index--by--geo--time.csv")) %>%
    rename(HDI = hdi_human_development_index)
  
  # Attended births
  attended_births_dt = fread(paste0(o$pth$input, "ddf--datapoints--births_attended_by_skilled_health_staff_percent_of_total--by--geo--time.csv")) %>%
    rename(attended_births = births_attended_by_skilled_health_staff_percent_of_total) %>%
    complete(geo, time = 1990:2024, 
             fill = list(doctors_per_1000 = NA)) %>%
    as_tsibble(index = time, 
               key = geo) 
  
  # Interpolate missing values
  attended_births_dt = attended_births_dt %>%
    model(lm = TSLM(log(attended_births) ~ trend())) %>%
    interpolate(attended_births_dt)
  
  # Internet users
  internet_users_dt = fread(paste0(o$pth$input, "ddf--datapoints--internet_users--by--geo--time.csv")) 
  
  
  # Create table of Gapminder covariates
  gapminder_dt = gini_dt %>%
    full_join(gdp_dt, by=c("geo", "time")) %>%
    full_join(doctors_per_1000_dt, by=c("geo", "time")) %>%
    full_join(pop_0_to_14_dt, by=c("geo", "time")) %>%
    full_join(pop_density_dt, by=c("geo", "time")) %>%
    full_join(urban_dt, by=c("geo", "time")) %>%
    full_join(hiv_mortality_dt, by=c("geo", "time")) %>%
    full_join(health_spending_dt, by=c("geo", "time")) %>%
    full_join(private_health_spending_dt, by=c("geo", "time")) %>%
    full_join(internet_users_dt, by=c("geo", "time")) %>%
    full_join(attended_births_dt, by=c("geo", "time")) %>%
    full_join(water_dt, by=c("geo", "time")) %>%
    full_join(sanitation_dt, by=c("geo", "time")) %>%
    full_join(hdi_dt, by=c("geo", "time")) %>%
    mutate(country_code = toupper(geo)) %>%
    select(-geo) %>%
    relocate(country_code) %>%
    rename(year = time) %>%
    full_join(WHO_regions_dt, by="country_code", relationship = "many-to-many") %>%
    select(-country) %>%
    rename(country = country_code) %>%
    arrange(country, year) %>%
    filter(year >= 1964 & year <= 2024)  %>%# include ten years before EPI to allow for historical effect of covariates
    as.data.table()              
  

  # Check for Gapminder countries not linked to WHO regions
 # gapminder_dt %>% filter(is.na(region_short)) %>%
#    select(country, region_short) %>%
#    unique()
  
  gapminder_dt = gapminder_dt %>%
    filter(!is.na(region_short))
  
  # Save in tables cache
  save_table(gapminder_dt, "gapminder_covariates")
}

# ---------------------------------------------------------
# Prepare country income status classification over time
# ---------------------------------------------------------
prepare_income_status = function() {
  
  message(" - Income status")
  
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
  
  message(" - Demography data")
  
  # Details of data sets to load
  data_sets = list(
    pop = list(name = "pop",    var = "pop"),
    mx  = list(name = "deaths", var = "mxB"))
  
  # Iterate through data sets to load
  for (id in names(data_sets)) {
    
    # Index data set name and variable reference
    name = data_sets[[id]]$name
    var  = data_sets[[id]]$var
    
    # Load population size data
    if (id == "pop") {
      
      # Names of WPP2022 data files to load
      past   = paste0("popAge",     o$pop_bin, "dt")
      future = paste0("popprojAge", o$pop_bin, "dt") 
      
      # Load pop data from WPP github package
      data_list = data_package(past, future, package = "wpp2022")
      
      # Scale metrics by factor of 1k
      scaler_dt = expand_grid(
        country = all_countries(), 
        year    = o$years, 
        age     = 0 : 100, 
        scaler  = 1e3) %>%
        as.data.table()
    }
    
    # Load any other data type
    if (id != "pop") {
      
      # Name of WPP2022 data file - history and projection in one
      all_time = paste0(id, o$pop_bin, "dt") 
      
      # Load data from WPP github package
      data_list = data_package(all_time, package = "wpp2022")
      
      # We'll need to scale per 1k population
      scaler_dt = table("wpp_pop") %>%
        rename(scaler = pop)
    }
    
    # Combine past and future data
    data_dt = rbindlist(data_list, fill = TRUE) %>%
      # Select countries of interst...
      inner_join(y  = table("country"),  
                 by = "country_code") %>%
      select(country, year, age, value = !!var) %>%
      # Shift year by one (see github.com/PPgp/wpp2022 for details)...
      mutate(year = as.integer(year) + 1) %>%
      filter(year %in% o$years) %>%
      # Scale metrics...
      left_join(y = scaler_dt, 
                by = c("country", "year", "age")) %>%
      mutate(value = value * scaler) %>%
      select(-scaler) %>%
      # Tidy up...
      rename(!!name := value) %>%
      arrange(country, year, age)
    
    # Save in tables cache
    save_table(data_dt, paste1("wpp", name))
  }
}

# ---------------------------------------------------------
# Prepare age at birth by country and year
# ---------------------------------------------------------
prepare_birth_age = function() {
  
  # Construct path to data file
  #
  # SOURCE: https://w3.unece.org/PXWeb/en/Table?IndicatorCode=34
  data_file = paste0(o$pth$input, "age_at_birth.csv")
  
  # Load raw data
  data_dt = fread(data_file) %>%
    select(country = Alpha3Code, 
           year    = PeriodCode, 
           value   = Value) %>%
    # Remove any unknown countries...
    filter(country %in% all_countries(), 
           year    %in% o$years) %>%
    # Format values into numeric...
    mutate(value = str_remove(value, "\\.\\."), 
           value = as.numeric(value)) %>%
    filter(!is.na(value)) %>%
    # Take the national mean over time...
    group_by(country) %>%
    summarise(avg = mean(value)) %>%
    ungroup() %>%
    # Expand to all countries...
    right_join(y  = all_countries(as_dt = TRUE),
               by = "country") %>%
    # Impute missing countries with global mean...
    mutate(avg = ifelse(
      test = is.na(avg), 
      yes  = mean(avg, na.rm = TRUE), 
      no   = avg)) %>%
    # Convert to integer...
    mutate(avg = round(avg)) %>%
    arrange(country) %>% 
    as.data.table()
  
  # Construct age x country matrix
  birth_age_mat = matrix(
    data = 0, 
    nrow = length(o$ages),
    ncol = nrow(data_dt))
  
  # Range of viable ages around the mean
  range = -(o$birth_age_sd * 2) : (o$birth_age_sd * 2)
  
  # Iterate through countries
  for (i in seq_row(data_dt)) {
    
    # Average age at birth
    avg = data_dt[i, avg]
    
    # Distribution around this mean
    dist = dnorm(
      x    = avg + range, 
      mean = avg, 
      sd   = o$birth_age_sd)
    
    # Insert these values into matrix
    birth_age_mat[avg + range, i] = dist / sum(dist)
  }
  
  # COnvert matrix into long datatable
  birth_age_dt = birth_age_mat %>%
    as_named_dt(data_dt$country) %>%
    mutate(age = o$ages) %>%
    # Melt to tidy format...
    pivot_longer(cols = -age, 
                 names_to  = "country", 
                 values_to = "weight") %>%
    # Remove trivial values...
    filter(weight > 0) %>%
    select(country, age, weight) %>%
    arrange(country, age) %>%
    as.data.table()
  
  # Save to file for easier loading
  save_table(birth_age_dt, "birth_age")
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

