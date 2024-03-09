###########################################################
# NON-MODELLED
#
# Estimate impact (deaths and DALYs averted) for static model
# pathogens using Global Burden of Disease estimates.
#
###########################################################

# ---------------------------------------------------------
# Parent function for calculating impact for static model pathogens
# ---------------------------------------------------------
run_static = function() {
  
  # Only continue if specified by do_step
  if (!is.element(3, o$do_step)) return()
  
  message("* Estimating static model vaccine impact")
  
  # Pathogens diseases of interest (those statically modelled)
  diseases = unique(table("d_v_a")[source == "static", disease])
  
  # Associated full names
  names = table("disease_name") %>%
    filter(disease %in% diseases) %>%
    select(d = disease, 
           x = disease_name)
  
  # Iterate through these diseases
  for (disease in diseases) {
    
    message(" - ", names[d == disease, x])
    
    # Effective coverage considering waning immunity and boosters
    effective_coverage(disease)

    # Deaths averted considering effective coverage and GBD disease burden
    deaths_averted(disease)
    
    # Use deaths averted to calculate DALYs
    # dalys_averted(id)  # See dalys.R
  }
  
  # Compile all results
  compile_outputs(diseases)
  
  # ---- Plot results ----
  
  # Effective coverage with waning immunity for static model pathogens
  #plot_effective_coverage()
  
  # Deaths and DALYs averted for static model pathogens
  #plot_static()
}

# ---------------------------------------------------------
# Effective coverage considering waning immunity and boosters
# ---------------------------------------------------------
effective_coverage = function(disease) {
  
  # CALCULATION PROCESS: 
  # 1) Number of immune people over time - considering waning immunity
  # 2) Number of people covered with different vaccine / schedules
  # 3) Weight immunity by primary & booster dosing
  
  schedule_id = c(
    x  = "primary",
    BX = "booster", 
    PX = "pregnancy")
  
  # ---- Set up ----
  
  # Vaccines targeting this disease
  vaccine_dt = table("d_v_a") %>% 
    filter(disease == !!disease) %>%
    separate(
      col  = vaccine, 
      into = qc(type, schedule), 
      sep  = "_", 
      fill = "right", 
      remove = FALSE) %>%
    replace_na(list(schedule = "x")) %>%
    mutate(schedule = recode(schedule, !!!schedule_id)) %>%
    mutate(type = str_remove(vaccine, "([0-9]+|_.+)$")) %>% 
    select(type, vaccine, schedule) %>% 
    as.data.table()
  
  # Associated full names
  names = table("vaccine_name") %>%
    filter(vaccine %in% vaccine_dt$vaccine) %>%
    select(v = vaccine, 
           x = vaccine_name)
  
  # Vaccine immunity efficacy profiles
  profile_dt = table("vaccine_efficacy_profiles") %>%
    filter(vaccine %in% vaccine_dt$vaccine) %>%
    pivot_wider(id_cols = time, 
                names_from  = vaccine,
                values_from = profile) %>%
    as.data.table()
  
  # Vaccine coverage
  coverage_dt = table("coverage") %>%
    left_join(y  = table("d_v_a"), 
              by = "d_v_a_id") %>%
    filter(vaccine %in% vaccine_dt$vaccine) %>%
    select(vaccine, country, year, age, fvps)
  
  # ---- Waning immunity per vaccine / schedule ----
  
  # Initiate list for effective immunity
  immunity_list = list()
  
  # First deal with non-pregnancy vaccines
  vaccines_nonpx = vaccine_dt[schedule != "pregnancy", vaccine]
  
  # Iterate through these vaccines
  for (vaccine in vaccines_nonpx) {

    message("  > ", names[v == vaccine, x])

    # Immunity effiacy profile for this vaccine
    profile = profile_dt[[vaccine]]

    # For each entry, apply waning immunity over time
    immunity_list[[vaccine]] = coverage_dt %>%
      filter(vaccine == !!vaccine) %>%
      dtapply(waning_immunity, profile) %>%
      rbindlist() %>%
      mutate(vaccine = vaccine)
  }
  
  # Then deal with vaccines given during pregnancy
  vaccines_px = vaccine_dt[schedule == "pregnancy", vaccine]
  
  # # Iterate through these vaccines
  for (vaccine in vaccines_px) {
    
    message("  > ", names[v == vaccine, x])
    
    # First effect: on the pregnant woman...
    message("   ~ Effect on pregnant women")
    
    # Maternal doses assumed to be boosters
    booster_ref = vaccine_dt %>%
      filter(vaccine == !!vaccine) %>%
      select(type) %>%
      paste1("BX")
    
    # Mothers covered
    mother_coverage = coverage_dt %>%
      filter(vaccine == !!vaccine) %>%
      select(country, year, age, fvps)
    
    # Immunity effiacy profile for mothers
    mother_profile = profile_dt[[booster_ref]]
    
    # Apply this immunity
    mother_immunity = mother_coverage %>%
      dtapply(waning_immunity, mother_profile) %>%
      rbindlist()
    
    # Second effect: on the neonate(s)...
    message("   ~ Effect on neo-nates")
    
    # Neonates covered
    neonate_coverage = mother_coverage %>%
      lazy_dt() %>%
      group_by(country, year) %>%
      summarise(fvps = sum(fvps)) %>%
      ungroup() %>%
      mutate(age = -1, .before = fvps) %>%
      as.data.table()
    
    # Immunity among neonates as defined in vaccine_efficacy table
    neonate_profile = profile_dt[[vaccine]]
    
    # Apply this immunity
    neonate_immunity = neonate_coverage %>%
      dtapply(waning_immunity, neonate_profile) %>%
      rbindlist() %>%
      filter(effective > 0)
    
    # Combine effective coverage
    immunity_list[[vaccine]] = neonate_immunity %>%
      rbind(mother_immunity) %>%
      mutate(vaccine = vaccine) %>%
      as.data.table()
  }
  
  message("  > Concatenating waning immunity")
  
  # Concatenate results of each vaccine / schedule
  immunity_dt = expand_grid(
    vaccine = vaccine_dt$vaccine, 
    country = all_countries(),
    year    = o$years, 
    age     = c(-1, o$ages)) %>%
    # Append results...
    left_join(y  = rbindlist(immunity_list), 
              by = names(.)) %>%
    replace_na(list(
      covered   = 0,
      effective = 0)) %>%
    # Summarise where multiple sources of immunity...
    lazy_dt() %>%
    group_by(vaccine, country, year, age) %>%
    summarise(covered   = sum(covered),
              effective = sum(effective)) %>%
    ungroup() %>%
    # Append vaccine details...
    left_join(y  = vaccine_dt,
              by = "vaccine") %>%
    select(type, vaccine, schedule, country,
           year, age, covered, effective) %>%
    arrange(type, vaccine, country, year, age) %>%
    as.data.table()
  
  # And save the result in it's own right
  save_rds(immunity_dt, "static_d", "immunity", disease)
  
  # ---- Overall effective coverage ----
  
  message("  > Calculating effective coverage")
  
  # Assume all pregnant doses are boosters (in terms of the mother)
  immunity_dt %<>%
    # Assume doses for pregnant women are boosters...
    mutate(schedule = ifelse(
      test = schedule == "pregnancy" & age >= 0,
      yes  = "booster",
      no   = schedule)) %>%
    # In terms of neonates, this would be primary...
    mutate(schedule = ifelse(
      test = schedule == "pregnancy" & age == -1,
      yes  = "primary",
      no   = schedule)) %>%
    # Summarise over primary-boosters...
    lazy_dt() %>%
    group_by(type, schedule, country, year, age) %>%
    summarise(covered   = sum(covered),
              effective = sum(effective)) %>%
    ungroup() %>%
    as.data.table()

  # Inititate list for effective coverage
  effective_list = list()
  
  # Iterate through vacicne types
  for (vaccine_type in unique(vaccine_dt$type)) {
    
    # Convert effective FVPs to effective coverage
    effective_v_dt = immunity_dt %>%
      filter(type == vaccine_type) %>%
      # Weight between schedules...
      weight_booster() %>%
      # Convert to effective coverage...
      left_join(y  = table("wpp_pop"),
                by = c("country", "year", "age")) %>%
      mutate(coverage = pmin(effective / pop, o$max_coverage)) %>%
      replace_na(list(coverage = 0)) %>%
      # Append vaccine and disease details and tidy up...
      mutate(disease = disease, 
             type    = vaccine_type) %>%
      select(disease, type, country, year, 
             age, effective, pop, coverage) %>%
      arrange(disease, type, country, year, age) %>%
      as.data.table()
    
    # Store result for this vaccine type in list
    effective_list[[vaccine_type]] = effective_v_dt
    
    # And save the result in it's own right
    save_rds(effective_v_dt, "static_t", "effective_coverage", vaccine_type)
  }
  
  # Combine vaccine types (considered additive)
  effective_d_dt = rbindlist(effective_list) %>%
    group_by(disease, country, year, age) %>%
    summarise(effective = sum(effective), 
              pop       = mean(pop)) %>%
    ungroup() %>%
    mutate(coverage = pmin(effective / pop, o$max_coverage)) %>%
    replace_na(list(coverage = 0)) %>%
    select(-effective, -pop) %>%
    as.data.table()
  
  # Save this result to file
  save_rds(effective_d_dt, "static_d", "effective_coverage", disease)
}

# ---------------------------------------------------------
# People effectively covered over time considering waning immunity
# ---------------------------------------------------------
waning_immunity = function(data, profile) {
  
  # Extend ages to include neo-nates
  all_ages = c(-1, o$ages)
  
  # Indicies for years and ages
  year_idx = match(data$year, o$years)  : length(o$years)
  age_idx  = match(data$age,  all_ages) : length(o$ages)
  
  # Index upto only the smallest of these two vectors
  vec_idx = 1 : min(length(year_idx), length(age_idx))
  
  # Effective FVPs: initial efficacy and immunity decay
  effective_fvps = data$fvps * profile
  
  # These form the only non-trivial entries
  waning_immunity_dt = data.table(
    country   = data$country, 
    year      = o$years[year_idx[vec_idx]],
    age       = all_ages[age_idx[vec_idx]], 
    covered   = data$fvps, 
    effective = effective_fvps[vec_idx])
  
  return(waning_immunity_dt)
}

# ---------------------------------------------------------
# Weight between primary & booster for total effective FVPs
# ---------------------------------------------------------
weight_booster = function(immunity_dt) {
  
  # Number covered by schedule (primary and booster)
  covered_dt = immunity_dt %>%
    select(country, year, age, covered, schedule) %>%
    pivot_wider(names_from  = schedule, 
                values_from = covered) %>%
    # Append trivial column if no booster doses...
    bind_rows(data.table(booster = double())) %>%
    replace_na(list(booster = 0)) %>%
    as.data.table()
  
  # Proportion of primary cases that have booster (capped at 100%)
  weighting_dt = covered_dt %>%
    mutate(primary = pmax(primary, booster), 
           weight  = pmin(booster / primary, 1)) %>%
    filter(!is.nan(weight)) %>%
    select(country, year, age, weight) %>%
    as.data.table()
  
  # Effective number of FVP after booster considerations
  weighted_dt = immunity_dt %>%
    select(country, year, age, effective, schedule) %>%
    pivot_wider(names_from  = schedule, 
                values_from = effective) %>%
    # Append trivial column if no booster doses...
    bind_rows(data.table(booster = double())) %>%
    replace_na(list(booster = 0)) %>%
    # Append booster weighting details...
    left_join(y  = weighting_dt, 
              by = c("country", "year", "age")) %>%
    replace_na(list(weight = 0)) %>%
    # Weight efficacies based on number covered...
    mutate(effective = primary * (1 - weight) + booster) %>% 
    select(country, year, age, effective) %>%
    as.data.table()
  
  return(weighted_dt)
}

# ---------------------------------------------------------
# Deaths averted considering effective coverage and GBD disease burden
# ---------------------------------------------------------
deaths_averted = function(disease) {
  
  message("  > Calculating deaths averted")
  
  # ---- Deaths averted for this disease ----
  browser()
  # Load effective coverage for this disease from file
  effective_dt = read_rds("static_d", "effective_coverage", disease)
  
  # Load disease deaths, append coverage, and estimate deaths averted
  averted_disease = effective_dt %>%
    lazy_dt() %>%
    inner_join(y  = table("gbd_estimates"),
               by = c("disease", "country", "year", "age")) %>%
    # Estimate deaths without a vaccine and deaths averted...
    mutate(deaths_without = deaths_disease / (1 - coverage), 
           deaths_averted = deaths_without - deaths_disease) %>%
    # Summairse neo-nate and other infant deaths...
    mutate(age = pmax(age, 0)) %>%
    group_by(disease, country, year, age) %>%
    summarise(deaths_disease = sum(deaths_disease), 
              deaths_averted = sum(deaths_averted)) %>%
    ungroup() %>%
    as.data.table()
  
  # Save this result to file
  save_rds(averted_disease, "static_d", "deaths_averted", disease)
  
  # ---- Attribute impact to vaccine schedule ----
  
  # Relative weighting of effective coverage for each vaccine
  weight_dt = read_rds("static_d", "immunity", disease) %>%
    lazy_dt() %>%
    left_join(y  = table("d_v_a"),
              by = "vaccine") %>%
    group_by(country, year, age) %>%
    mutate(weight = effective / sum(effective)) %>%
    ungroup() %>%
    filter(!is.nan(weight)) %>%
    select(d_v_a_id, vaccine, country, year, age, weight) %>%
    as.data.table()
  
  # We'll split impact across schedules by previously calculated weighting
  averted_vaccine = averted_disease %>%
    lazy_dt() %>%
    # Attribute impact to each shedule...
    inner_join(y  = weight_dt, 
               by = c("country", "year", "age")) %>%
    # Summarise over age ...
    group_by(country, d_v_a_id, year) %>%
    summarise(impact = sum(deaths_averted * weight)) %>%
    ungroup() %>%
    # Tidy up...
    select(country, d_v_a_id, year, impact) %>%
    arrange(country, d_v_a_id, year) %>%
    as.data.table()
  
  # Number of 'FVP' by vaccine - yeah, confusing concept in the context of boosters
  fvps_dt = table("coverage") %>%
    lazy_dt() %>%
    left_join(y  = table("d_v_a"), 
              by = c("d_v_a_id")) %>%
    filter(disease == !!disease) %>%
    # Summarise over age ...
    group_by(d_v_a_id, country, year) %>%
    summarise(fvps = sum(fvps)) %>%
    ungroup() %>%
    as.data.table()
  
  # Full vaccine-specific results 
  result_list = averted_vaccine %>%
    inner_join(y  = fvps_dt, 
               by = c("d_v_a_id", "country", "year")) %>%
    select(d_v_a_id, country, year, fvps, impact) %>%
    split(.$d_v_a_id)
  
  # Function for saving vaccine-specific results
  save_result_fn = function(x, name)
    save_rds(x, "static_v", "deaths_averted", name)
  
  # Apply save function to each vaccine
  napply(result_list, save_result_fn)
}

# ---------------------------------------------------------
# Compile and save all static model outputs 
# ---------------------------------------------------------
compile_outputs = function(diseases) {
  
  message(" - Compiling results")
  
  # IDs targeting these disease
  d_v_a_id = table("d_v_a") %>% 
    filter(disease %in% diseases) %>%
    pull(d_v_a_id)
  
  # Vaccine types targeting these disease
  vaccine_types = table("d_v_a") %>% 
    filter(disease %in% diseases) %>%
    mutate(type = str_remove(vaccine, "([0-9]+|_.+)$")) %>% 
    pull(type) %>% 
    unique()
  
  # IDs for loading disease/vaccine results
  ids = list(
    disease = diseases,
    vaccine = d_v_a_id, 
    type    = vaccine_types)
  
  # All outputs to compile
  outputs = list(
    disease = qc(deaths_averted, effective_coverage), #, dalys_averted),
    vaccine = qc(deaths_averted), #, dalys_averted),
    type    = qc(effective_coverage))
  
  for (group in names(outputs)) {
    
    # Sub directory to load disease/vaccine results from
    load_dir = paste1("static", str_sub(group, 1, 1))
    
    # Function to read rds files
    load_fn = function(file) 
      read_rds(load_dir, file)
    
    # Loop through all outputs
    for (output in outputs[[group]]) {
      
      # Load all files related to this output
      output_dt = paste1(output, ids[[group]]) %>%
        lapply(load_fn) %>%
        rbindlist()
      
      # Save compiled results to file
      save_rds(output_dt, "static", output, group)
    }
  }
}

