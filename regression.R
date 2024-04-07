###########################################################
# REGRESSION
#
# Fit a user-defined series of models to cumulative impact 
# per FVP. This process has two use cases:
#  1) Impute impact for countries not modelled by VIMC
#  2) To infer drivers of impact (on complete set of estimates)
#
###########################################################

# ---------------------------------------------------------
# Parent function for regression modelling
# ---------------------------------------------------------
run_regression = function(case, metric) {

  # Only continue if specified by do_step
  if (case == "impute" & !is.element(4, o$do_step)) return()
  if (case == "infer"  & !is.element(7, o$do_step)) return()
  
  message("* Running regression: ", case, " ", metric)
  
  # TEMP: Use basic or full imputation method
  #
  # OPTIONS:
  #  basic_regression   - IA2030 method using GBD covariates
  #  perform_regression - Helen's time series regression, refactored code
  method = "perform_regression"
  
  # TEMP: Ignoring problematic cases for now
  ignore = NULL
  # c(12,  # Non-routine, small numbers
#   11,   # Yellow fever
#  16, 18, 22)  # DTP boosters causing problems

  # ---- Load data ----
  
  # Load response variable (impact per FVP)
  target = get_regression_data(case, metric)
  
  # Return out if no training data identified
  if (nrow(target) == 0)
    return()
  
  # ---- Perform regression ----
  
  # Which sources of public health impact are to be modelled
  if (case == "impute") use_sources = qc(vimc)
  if (case == "infer")  use_sources = qc(vimc, static, extern)
  
  # Call country imputation function
  predict_dt = table("d_v_a") %>%
    filter(source %in% use_sources) %>%
    # TEMP: Ignoring problematic cases for now...
    filter(!d_v_a_id %in% ignore) %>%  
    # Apply geographical imputation model...
    pull(d_v_a_id) %>%
    lapply(FUN = get(method), 
           target = target, 
           case   = case, 
           metric = metric) %>%
    rbindlist() %>%
    select(d_v_a_id, country, year, impact_impute)
  
  # ---- Use regression to impute missing countries ----
  
  # Apply imputations where needed
  # 
  # NOTE: A trivial process in the infer case
  impute_dt = target %>%
    select(d_v_a_id, country, year, fvps_cum, impact_cum) %>%
    # Append predictions...
    left_join(y  = predict_dt, 
              by = c("d_v_a_id", "country", "year")) %>%
    # Apply predictions where imputation is needed...
    mutate(impact = ifelse(
      test = is.na(impact_cum),
      yes  = impact_impute,
      no   = impact_cum)) %>%
    select(d_v_a_id, country, year, 
           fvps = fvps_cum, impact) %>%
    # Assume any missing values are zero impact...
    replace_na(list(impact = 0)) %>%
    # TEMP: Bound impact above by 1...
    mutate(impact = pmin(impact, 1))  # HCJ: This can be removed when fitting looks good
  
  # Save imputed results to file
  save_rds(impute_dt, case, case, metric, "result")
  
}

# ---------------------------------------------------------
# Define initial regression models to evaluate inclusion of lagged vaccination coverage
# ---------------------------------------------------------
define_coverage_models = function(case) {
  
  # List of available covariates
  covars = list(
    cov0  = "log(coverage)",
    cov1  = "log(coverage_minus_1)", 
    cov2  = "log(coverage_minus_2)", 
    cov3  = "log(coverage_minus_3)", 
    cov4  = "log(coverage_minus_4)") 
   
  # Define models 
  models = list(
    
    # Models for imputing missing countries
    impute = list(
      m1  = "cov0", 
      m2  = "cov0 + cov1", 
      m3  = "cov0 + cov1 + cov2", 
      m4  = "cov0 + cov1 + cov2 + cov3", 
      m5  = "cov0 + cov1 + cov2 + cov3 + cov4"), 
    
    # Models for inferring key drivers of impact 
    infer = list(
      m1  = "cov0", 
      m2  = "cov0 + cov1", 
      m3  = "cov0 + cov1 + cov2", 
      m4  = "cov0 + cov1 + cov2 + cov3", 
      m5  = "cov0 + cov1 + cov2 + cov3 + cov4")) 
  
  return(list(models[[case]], covars))
}

# ---------------------------------------------------------
# Define set of regression models to evaluate
# ---------------------------------------------------------
define_models = function(case) {
  
  # List of available covariates
  covars = list(
    cov0  = "log(coverage)",
    cov1  = "log(coverage_minus_1)", 
    cov2  = "log(coverage_minus_2)", 
    cov3  = "log(coverage_minus_3)", 
    cov4  = "log(coverage_minus_4)", 
    gini  = "gini", 
    gdp   = "gdp",
    lit_f = "lit_female", 
    lit_m = "lit_male",
    doc   = "doctors_per_1000",
    dens  = "pop_density",
    urban = "urban_percent",
    mat   = "maternal_mortality",
    stunt = "stunting",
    phs   = "private_health",
    san   = "basic_sanitation",
    wat   = "basic_water",
    ab    = "attended_births", 
    hs   = "health_spending", 
    hs1   = "health_spending_minus_1", 
    hs2   = "health_spending_minus_2")
  
 # if(phase == 1){
  # Define models (using shorthand covariate references)
#  models = list(
    
    # Models for imputing missing countries
 #   impute = list(
#      m1  = "cov0 + cov1 + cov2 + cov3 + cov4",
#      m2  = "cov0 + cov1 + cov2 + cov3 + cov4 + gini",
#      m3  = "cov0 + cov1 + cov2 + cov3 + cov4 + gdp",
#      m4  = "cov0 + cov1 + cov2 + cov3 + cov4 + lit_f",
#      m5  = "cov0 + cov1 + cov2 + cov3 + cov4 + lit_m",
#      m6  = "cov0 + cov1 + cov2 + cov3 + cov4 + doc",
#      m7  = "cov0 + cov1 + cov2 + cov3 + cov4 + dens",
#      m8  = "cov0 + cov1 + cov2 + cov3 + cov4 + urban",
#      m9  = "cov0 + cov1 + cov2 + cov3 + cov4 + stunt",
#      m10  = "cov0 + cov1 + cov2 + cov3 + cov4 + mat",
#      m11  = "cov0 + cov1 + cov2 + cov3 + cov4 + phs",
#      m12  = "cov0 + cov1 + cov2 + cov3 + cov4 + san",
#      m13  = "cov0 + cov1 + cov2 + cov3 + cov4 + wat",
#      m14  = "cov0 + cov1 + cov2 + cov3 + cov4 + ab",
#      m15  = "cov0 + cov1 + cov2 + cov3 + cov4 + hs"))
      
 #     return(list(models[[case]], covars))
#  }
      
#  if(phase == 2){
    # Define models (using shorthand covariate references)
   # models = list(
      # Models for imputing missing countries
    #  impute = list(
     # m16 = "cov0 + cov1 + cov2 + cov3 + cov4 + gini",
  #    m17  = "cov0 + cov1 + cov2 + cov3 + cov4 + gini + gdp",
  #    m18  = "cov0 + cov1 + cov2 + cov3 + cov4 + gini + lit_f",
  #    m19  = "cov0 + cov1 + cov2 + cov3 + cov4 + gini + lit_m",
  #    m20  = "cov0 + cov1 + cov2 + cov3 + cov4 + gini + doc",
    #  m21  = "cov0 + cov1 + cov2 + cov3 + cov4 + gini + dens",
  #    m22 = "cov0 + cov1 + cov2 + cov3 + cov4 + gini + urban",
     # m23 = "cov0 + cov1 + cov2 + cov3 + cov4 + gini + stunt",
    #  m24  = "cov0 + cov1 + cov2 + cov3 + cov4 + gini + mat",
  #    m25  = "cov0 + cov1 + cov2 + cov3 + cov4 + gini + phs",
  #    m26  = "cov0 + cov1 + cov2 + cov3 + cov4 + gini + san",
     # m27  = "cov0 + cov1 + cov2 + cov3 + cov4 + gini + wat",
  #    m28  = "cov0 + cov1 + cov2 + cov3 + cov4 + gini + ab",
    #  m29 = "cov0 + cov1 + cov2 + cov3 + cov4 + gini + hs",
      
  #    m30  = "cov0 + cov1 + cov2 + cov3 + cov4 + dens + gdp",
  #    m31  = "cov0 + cov1 + cov2 + cov3 + cov4 + dens + lit_f",
  #    m32  = "cov0 + cov1 + cov2 + cov3 + cov4 + dens + lit_m",
  #    m33  = "cov0 + cov1 + cov2 + cov3 + cov4 + dens + doc",
  #    m34  = "cov0 + cov1 + cov2 + cov3 + cov4 + dens + urban",
     # m35  = "cov0 + cov1 + cov2 + cov3 + cov4 + dens + stunt",
    #  m36  = "cov0 + cov1 + cov2 + cov3 + cov4 + dens + mat",
  #    m37  = "cov0 + cov1 + cov2 + cov3 + cov4 + dens + phs",
  #    m38  = "cov0 + cov1 + cov2 + cov3 + cov4 + dens + san",
     # m39  = "cov0 + cov1 + cov2 + cov3 + cov4 + dens + wat",
  #    m40  = "cov0 + cov1 + cov2 + cov3 + cov4 + dens + ab",
      #m41  = "cov0 + cov1 + cov2 + cov3 + cov4 + dens + hs",
  #    
  #    m42  = "cov0 + cov1 + cov2 + cov3 + cov4 + stunt + gdp",
 #     m43  = "cov0 + cov1 + cov2 + cov3 + cov4 + stunt + lit_f",
#      m44  = "cov0 + cov1 + cov2 + cov3 + cov4 + stunt + lit_m",
  #    m45  = "cov0 + cov1 + cov2 + cov3 + cov4 + stunt + doc",
  #    m46  = "cov0 + cov1 + cov2 + cov3 + cov4 + stunt + urban",
      #m47  = "cov0 + cov1 + cov2 + cov3 + cov4 + stunt + mat",
  #    m48  = "cov0 + cov1 + cov2 + cov3 + cov4 + stunt + phs",
  #    m49  = "cov0 + cov1 + cov2 + cov3 + cov4 + stunt + san",
      #m50  = "cov0 + cov1 + cov2 + cov3 + cov4 + stunt + wat",
  #    m51  = "cov0 + cov1 + cov2 + cov3 + cov4 + stunt + ab",
      #m52  = "cov0 + cov1 + cov2 + cov3 + cov4 + stunt + hs",

  #    m53  = "cov0 + cov1 + cov2 + cov3 + cov4 + mat + gdp",
  #    m54  = "cov0 + cov1 + cov2 + cov3 + cov4 + mat + lit_f",
  #    m55  = "cov0 + cov1 + cov2 + cov3 + cov4 + mat + lit_m",
  #    m56  = "cov0 + cov1 + cov2 + cov3 + cov4 + mat + doc",
  #    m57  = "cov0 + cov1 + cov2 + cov3 + cov4 + mat + urban",
  #    m58  = "cov0 + cov1 + cov2 + cov3 + cov4 + mat + phs",
  #    m59  = "cov0 + cov1 + cov2 + cov3 + cov4 + mat + san",
      #m60  = "cov0 + cov1 + cov2 + cov3 + cov4 + mat + wat",
  #    m61  = "cov0 + cov1 + cov2 + cov3 + cov4 + mat + ab",
      #m62  = "cov0 + cov1 + cov2 + cov3 + cov4 + mat + hs",
      
  #    m63  = "cov0 + cov1 + cov2 + cov3 + cov4 + wat + gdp",
  #    m64  = "cov0 + cov1 + cov2 + cov3 + cov4 + wat + lit_f",
   #   m65  = "cov0 + cov1 + cov2 + cov3 + cov4 + wat + lit_m",
  #    m66  = "cov0 + cov1 + cov2 + cov3 + cov4 + wat + doc",
  #    m67  = "cov0 + cov1 + cov2 + cov3 + cov4 + wat + urban",
  #    m68  = "cov0 + cov1 + cov2 + cov3 + cov4 + wat + phs",
  #    m69  = "cov0 + cov1 + cov2 + cov3 + cov4 + wat + san",
  #    m70  = "cov0 + cov1 + cov2 + cov3 + cov4 + wat + ab",
      #m71  = "cov0 + cov1 + cov2 + cov3 + cov4 + wat + hs",
  #   
    # Models for inferring key drivers of impact 
    #infer = list(
    #  x1  = "cov0", 
    #  x8  = "cov0 + cov1 + cov2 + cov3 + cov4 + pop14 + gini + ab", 
    #  x13 = "cov0 + hdi + pop14 + gini"))
  
 # return(list(models[[case]], covars))}
  
  #if(phase == 3){
    # Define models (using shorthand covariate references)
   # models = list(
      # Models for imputing missing countries
      #impute = list(
       # m301 = "cov0 + cov1 + cov2 + cov3 + cov4 + gini",
      #  m302 = "cov0 + cov1 + cov2 + cov3 + gini",
      #  m303 = "cov0 + cov1 + cov2 + gini",
      #  m304 = "cov0 + cov1 + gini",
      #  m305 = "cov0 + gini",
        
       # m306  = "cov0 + cov1 + cov2 + cov3 + cov4 + dens",
      #  m307  = "cov0 + cov1 + cov2 + cov3 + dens",
      #  m308  = "cov0 + cov1 + cov2 + dens",
      #  m309  = "cov0 + cov1 + dens",
      #  m310  = "cov0 + dens",
        
       # m311  = "cov0 + cov1 + cov2 + cov3 + cov4 + stunt",
      #  m312  = "cov0 + cov1 + cov2 + cov3 + stunt",
      #  m313  = "cov0 + cov1 + cov2 + stunt",
      #  m314  = "cov0 + cov1 + stunt",
      #  m315  = "cov0 + stunt",
       
       # m316  = "cov0 + cov1 + cov2 + cov3 + cov4 + mat",
      #  m317  = "cov0 + cov1 + cov2 + cov3 + mat",
      #  m318  = "cov0 + cov1 + cov2 + mat",
      #  m319  = "cov0 + cov1 + mat",
      #  m320  = "cov0 + mat",
      # 
       # m321  = "cov0 + cov1 + cov2 + cov3 + cov4 + wat",
      #  m322  = "cov0 + cov1 + cov2 + cov3 + wat",
      #  m323  = "cov0 + cov1 + cov2 + wat",
      #  m324  = "cov0 + cov1 + wat",
      #  m325  = "cov0 + wat"))
       
    #return(list(models[[case]], covars))
    
    #if(phase == 4){
    # Define models (using shorthand covariate references)
     models = list(
    # Models for imputing missing countries
    impute = list(
      m1  = "cov0", 
      m2  = "cov0 + cov1", 
      m3  = "cov0 + cov1 + cov2", 
      m4  = "cov0 + cov1 + cov2 + cov3", 
      m5  = "cov0 + cov1 + cov2 + cov3 + cov4", 
      
    m401 = "cov0 + cov1 + cov2 + cov3 + mat",
    m402 = "cov0 + cov1 + cov2 + cov3 + mat + gini",
    m403 = "cov0 + cov1 + cov2 + cov3 + mat + gini + stunt",
    m404 = "cov0 + cov1 + cov2 + cov3 + mat + gini + stunt + phs",
    m405 = "cov0 + cov1 + cov2 + cov3 + mat + gini + stunt + dens",
    m406 = "cov0 + cov1 + cov2 + cov3 + mat + gini + stunt + wat"))
  
return(list(models[[case]], covars))
  
}

# ---------------------------------------------------------
# Load/calculate target variable (impact per FVP)
# ---------------------------------------------------------
get_regression_data = function(case, metric) {
  
  # Population size of each country over time
  pop_dt = table("wpp_pop") %>%
    lazy_dt() %>%
    group_by(country, year) %>%
    summarise(pop = sum(pop)) %>%
    ungroup() %>%
    as.data.table()
  
  # Impact estimates in imputation case: VIMC pathogens, VIMC countries
  if (case == "impute") {
    outcomes_dt = table("vimc_estimates") %>%
      rename(impact = !!paste1(metric, "averted"))
  }
  
  # Impact estimates in imputation case: all modelled results
  if (case == "infer")
    outcomes_dt = read_rds("history", metric, "averted")
  
  # TODO: Probably best to remove previously imputed estimates in the infer case
  
  # Convert estimates to cumulative form
  impact_dt = outcomes_dt %>%
    lazy_dt() %>%
    # Sum impact over age...
    group_by(d_v_a_id, country, year) %>%
    summarise(impact_abs = sum(impact)) %>%
    ungroup() %>%
    mutate(impact_abs = pmax(impact_abs, 0)) %>%
    # Scale results to per capita...
    left_join(y  = pop_dt, 
              by = c("country", "year")) %>%
    mutate(impact_rel = impact_abs / pop) %>%
    select(-pop) %>%
    # Cumulative sum impact...
    arrange(d_v_a_id, country, year) %>%
    group_by(d_v_a_id, country) %>%
    mutate(impact_cum = cumsum(impact_rel)) %>%
    ungroup() %>%
    as.data.table()
  
  # Extract FVPs
  fvps_dt = table("coverage") %>%
    lazy_dt() %>%
    # Subset pathogens...
    filter(d_v_a_id %in% unique(impact_dt$d_v_a_id)) %>%
    # Summarise over age...
    group_by(d_v_a_id, country, year) %>%
    summarise(fvps_abs = sum(fvps)) %>%
    ungroup() %>%
    # Scale results to per capita...
    left_join(y  = pop_dt, 
              by = c("country", "year")) %>%
    mutate(fvps_rel = fvps_abs / pop) %>%
    select(-pop) %>%
    # Cumulative sum FVPs...
    arrange(d_v_a_id, country, year) %>%
    group_by(d_v_a_id, country) %>%
    mutate(fvps_cum = cumsum(fvps_rel)) %>%
    ungroup() %>%
    as.data.table()
  
  # Combine into single datatable
  target_dt = fvps_dt %>%
    left_join(y  = impact_dt, 
              by = c("d_v_a_id", "country", "year")) %>%
    # Impact per FVP...
    mutate(target = impact_cum / fvps_cum)
  
  # Save this datatable to file for plotting purposes
  save_rds(target_dt, case, case, metric, "target")
  
  # Throw a warning if no target data identified
  if (nrow(target_dt) == 0)
    warning("No training data identified")
  
  return(target_dt)
}

# ---------------------------------------------------------
# Perform regression (IA2030 style approach)
# ---------------------------------------------------------
basic_regression = function(d_v_a_id, target, case, metric) {
  
  # Details of this d_v_a
  d_v_a_name = data.table(d_v_a_id = d_v_a_id) %>%
    format_d_v_a_name() %>%
    pull(d_v_a_name)
  
  # Display progress message to user
  message(" > ", d_v_a_name)
  
  # ---- Append covariates ----
  
  # All-cause death rate by country, year, and age
  infant_mortality_dt = table("wpp_pop") %>%
    filter(age == 0) %>%  # TODO: Instead filter by table("wiise_vaccine") age group
    inner_join(y  = table("wpp_deaths"),
               by = c("country", "year", "age")) %>%
    # Calculate infant mortality rate...
    mutate(imr = deaths / pop) %>%
    select(country, year, imr)
  
  # Compile datatable of covariates
  covariates_dt = table("regression_covariates") %>%
    filter(metric %in% c("haqi", "sdi")) %>%
    pivot_wider(names_from = metric) %>%
    full_join(y  = infant_mortality_dt, 
              by = c("country", "year")) %>%
    arrange(country, year) %>%
    group_by(country) %>%
    fill(haqi, .direction = "up") %>%
    ungroup() %>%
    as.data.table()
  
  # Append covariates to target
  target_dt = target %>%
    filter(d_v_a_id == !!d_v_a_id) %>%
    # Append GBD indices and infant mortality...
    left_join(y  = covariates_dt,
              by = c("country", "year")) %>%
    # Calculate n years of estimates...
    mutate(n_years = 1) %>%
    group_by(country) %>%
    mutate(n_years = cumsum(n_years)) %>%
    ungroup() %>%
    as.data.table()
  
  # TODO: We can either include or exclude zero here - arguably with zero is better...
  
  # Data used to fit statistical model
  data_dt = target_dt %>%
    filter(!is.na(target)) %>%
    # filter(target > 0) %>%
    select(target, n_years, sdi, haqi, imr) %>%
    # Remove target outliers for better normalisation...
    mutate(lower = mean(target) - 3 * sd(target), 
           upper = mean(target) + 3 * sd(target), 
           outlier = target < lower | target > upper) %>%
    filter(outlier == FALSE) %>%
    select(-outlier, -lower, -upper)
  
  # Sanity check that we have no NAs here
  if (any(is.na(data_dt)))
    stop("NA values identified in predictors")
  
  # Values to predict for (including data used for fitting)
  pred_dt = target_dt %>%
    select(all_names(data_dt))
  
  # ---- Check for trivial case ----
  
  # Return out if no data available
  if (nrow(data_dt[target > 0]) < 10) {
    
    message(" !! Insufficient data for imputation !!")
    
    # Store trivial outcomes
    fit = list(data = data_dt, result = NULL)
    
    # Save to file
    save_rds(fit, "impute", "impute", metric, d_v_a_id)
    
    return()
  }
  
  # ---- Normalise predictors and response ----
  
  # Function to normalise ready for fitting
  transform_fn = function(x, a, b)
    y = t((x - a) / (b - a)) %>% as.data.table()
  
  # Function to back transform to original scale
  retransform_fn = function(y, a, b)
    x = y * (b["target"] - a["target"]) + a["target"]
  
  # Matrices of points to fit with and points to predict for
  data_mat = t(as.matrix(data_dt))
  pred_mat = t(as.matrix(pred_dt))
  
  # Min and max in data used for fitting
  a = rowMins(data_mat)
  b = rowMaxs(data_mat)
  
  # Use these min ana max values to normalise
  norm_data_dt = transform_fn(data_mat, a, b)
  norm_pred_dt = transform_fn(pred_mat, a, b)
  
  # ---- Fit a model to predict impact per FVP ----
  
  # Fit a GLM for impact per FVP using all covariates
  fit_model = glm(
    formula = target ~ n_years + sdi + haqi + imr, 
    data    = norm_data_dt)
  
  # Use fitted model to predict 
  result_dt = target_dt %>%
    select(country, d_v_a_id, year, fvps_cum, impact_cum) %>%
    # Predict impact per FVP...
    cbind(norm_pred_dt) %>%
    mutate(predict = predict(fit_model, .), 
           predict = pmax(predict, 0)) %>%
    # Remove predictors...
    select(country, d_v_a_id, year, fvps_cum, impact_cum, 
           target, predict) %>%
    # Back-transform target and prediction...
    mutate(target  = retransform_fn(target,  a, b), 
           predict = retransform_fn(predict, a, b)) %>%
    # Multiply through to obtain cumulative impact over time...
    mutate(impact_impute = fvps_cum * predict, 
           .after = impact_cum)
  
  # Sanity check that all predicted values are legitimate
  if (any(is.na(result_dt$predict)))
    stop("NA values identified in predicted impact")
  
  # Store the fitted model, the data used, and the result
  fit = list(
    model   = fit_model, 
    data    = norm_data_dt, 
    result  = result_dt)
  
  # Save to file
  save_rds(fit, case, case, metric, d_v_a_id)
  
  return(result_dt)
}

# ---------------------------------------------------------
# Perform regression
# ---------------------------------------------------------
perform_regression = function(d_v_a_id, target, case, metric) {

 # Extract name of this d-v-a
  d_v_a_name = table("d_v_a") %>%
    filter(d_v_a_id == !!d_v_a_id) %>%
    pull(d_v_a_name)
  
  # Display progress message to user
  message(" > ", d_v_a_name)

  # Load set of models to evaluate
  list[models, covars] = define_models(case)
  
  # Append all required covariates - see separate function
  target_ts = append_covariates(d_v_a_id, models, covars, target)
  
  # Income groupings
  # Load income status dictionary
  income_dict = table("income_dict")
  
  # Load income status of each country
  income_dt = table("income_status") %>%
    filter(year == max(year)-5) %>% # Income level 5 years ago
    left_join(y  = income_dict, 
              by = "income") %>%
    select(country, income = income_name)
  
  
  # ---- Evaluate all user-defined models ----
  
  message("  - Evaluating models")
  
  # Subset training data (which we have impact estimates for)
  data_ts = target_ts %>%
    # Remove zeros to allow for log transformation...
    filter(target > 0) %>%
    # Remove country if insufficient data points for fitting...
    group_by(country) %>%
    filter(n() >= o$min_data_requirement) %>% 
    ungroup()
 
  # Evaluate all models in parallel
  if (o$parallel$impute)
    model_list = mclapply(
      X   = names(models),
      FUN = evaluate_model, 
      models = models, 
      covars = covars, 
      data   = data_ts,
      mc.cores = o$n_cores,
      mc.preschedule = FALSE)
  
  # Evaluate all models consecutively
  if (!o$parallel$impute)
    model_list = lapply(
      X   = names(models),
      FUN = evaluate_model, 
      models = models, 
      covars = covars, 
      data   = data_ts)
  
  # ---- Model selection ----
   message("  - Model selection")
  
  # For each country, select the model with the best AICc
  model_choice = model_list %>%
    lapply(report) %>%
    rbindlist() %>%
    # Remove null models...
    filter(!is.infinite(AICc)) %>% 
    # Retain only best fit model (if equal, keep the first)...
    group_by(country) %>%
    slice_min(AICc, with_ties = FALSE) %>%
    unique() %>%
    # Reappend best model...
    left_join(y  = rbindlist(model_list), 
              by = c("country", "model_id")) %>%
    # Reduce down to keep only model and AICc... 
    select(country, model_id, tslm, AICc) %>%
   # mutate(model_id = as.factor(model_id)) %>%
    mutate(d_v_a_id = d_v_a_id, 
           .before = 1) %>%
    # Convert to mable class...
    as_mable(key   = country, 
             model = tslm) %>%
    suppressWarnings() %>%
    append_region_name() %>%
    left_join(y = income_dt,
              by = "country")
  
    # Extract parameters of best fitting model for each country
  model_fit = tidy(model_choice) %>%
    select(d_v_a_id, country, model_id, term, 
           estimate, std.error, p.value) %>%
    as.data.table()
  
  # ---- Predictions ----
  
  message("  - Model predictions")
 
  # Impute case: predict missing data using regional best models 
  if (case == "impute") {
    
    # Evaluate this model - see separate function
    predict_dt = evaluate_predictions(
      model_choice = model_choice, 
      model_list = model_list, 
      target     = target_ts, 
      case       = case,
      income_dt     = income_dt)
  }
  
  # Infer case: just a case of evaluating on the training data
  if (case == "infer") {
  
    # Evaluate models on the training data
    predict_dt = augment(model_choice) %>%
      select(country, year, prediction = .fitted) %>%
      as.data.table()
  }
  
  # ---- Format output ----
  
  # Apply predictions to impute missing impact estimates
  result_dt = target %>%
    filter(d_v_a_id == !!d_v_a_id) %>%
    # Append predictions...
    left_join(y  = predict_dt, 
              by = c("country", "year")) %>%
    select(d_v_a_id, country, year, fvps_cum, 
           impact_cum, target, prediction) %>%
    # Multiply through to obtain cumulative impact over time...
    mutate(impact_impute = fvps_cum * prediction, 
           .after = impact_cum)
  
  # Also format predictors for use in plotting
  data_dt = data_ts %>%
    mutate(d_v_a_id = d_v_a_id, 
           .before = 1) %>%
    as.data.table()
  
  # Store the data used, fitted model, and result
  fit = list(
    model  = model_choice,  # NOTE: Only for non-imputed
    report = model_fit,
    result = result_dt, 
    data   = data_dt)
  
  # Save to file
  save_rds(fit, case, case, metric, d_v_a_id)
  
  return(result_dt)
}

# ---------------------------------------------------------
# Append all required covariates
# ---------------------------------------------------------
append_covariates = function(d_v_a_id, models, covars, target) {
  
  # TODO: Are GBD covariates still used/needed?
  
  # ---- Identify covariates from model specification ----
  
  # Shorthand covariates used in specified models
  covars_used = unlist(models) %>%
    paste(collapse = " + ") %>%
    str_split_1(pattern = " \\+ ") %>%
    unique()
  
  # Associated names of covariate columns 
  covars_retain = covars[covars_used] %>%
    unlist(covars) %>%
    str_remove("^.*\\(+") %>%
    str_remove("\\)+$")
  
  # ---- Define covariates to be lagged ----
  
  # Details of covariates we wish to lag
  lag_dt = covars_retain %>%
    str_split("_minus_", simplify = TRUE) %>%
    as_named_dt(c("covar", "idx")) %>%
    filter(idx > 0)
  
  # Extract all to-be-lagged covariates
  covars_lag  = unique(lag_dt$covar)
  n_lag_years = max(as.numeric(lag_dt$idx))
  
  # Small function to apply lag to given covariate
  covar_lag_fn = function(dt) {
    
    # Iterate through covariates and years to lag
    for (i in covars_lag) {
      for (j in seq_len(n_lag_years)) {
        
        # Incrementally offset by one year
        dt[[paste1(i, "minus", j)]] = lag(dt[[i]], j)
      }
    }
    
    return(dt)
  }
  
  # ---- Format coverage (a key predictor) ----
  
  # Summarise vaccination coverage by country and year
  coverage_dt = table("coverage") %>%
    lazy_dt() %>%
    filter(d_v_a_id == !!d_v_a_id) %>%
    # Summarise over age groups...
    group_by(country, year) %>%
    summarise(fvps   = sum(fvps),
              cohort = sum(cohort)) %>%
    ungroup() %>%
    # Recalculate for whole population...
    mutate(coverage = fvps / cohort) %>%
    select(country, year, coverage) %>%
    as.data.table()
  
  # ---- Append all other covariates ----
  
  # Spread covariates to wide format
  covariates_dt = table("regression_covariates") %>%
    pivot_wider(names_from = metric) %>%
    as.data.table()
  
  # Create time-series tibble with all covariates
  target_ts = target %>%
    # Data for this d-v-a...
    filter(d_v_a_id == !!d_v_a_id) %>%
    select(-d_v_a_id) %>%
    # Append vaccination coverage...
    full_join(y = coverage_dt,  
              by = c("country", "year")) %>%
    # Append all possible covariates...
    full_join(y  = covariates_dt,
              by = c("country", "year")) %>%
    arrange(country, year) %>%
    # Lag any necessary covariates...
    group_by(country) %>%
    covar_lag_fn() %>%
    ungroup() %>%
    # Only retain covariates defined in models...
    select(country, year, target, 
           all_of(covars_retain)) %>%
    # Convert to time-series tibble...
    as_tsibble(index = year, 
               key   = country) 
  
  return(target_ts)
}

# ---------------------------------------------------------
# Evaluate given user-specified model
# ---------------------------------------------------------
evaluate_model = function(id, models, covars, data) {
  
  # Interpret covariate references
  model_str = interpret_covars(models[[id]], covars)
  
  # Construct full model string to be evaluated
  model_fn   = paste0("TSLM(log(target) ~ ", model_str, ")")
  model_eval = paste0("model(data, tslm = ", model_fn, ")")
  
  # Evaluate model and append model reference
  model_mab = eval_str(model_eval) %>%
    mutate(model_id = id, 
           .before  = tslm) %>%
    suppressWarnings()
  
  return(model_mab)
}

# ---------------------------------------------------------
# Evaluate chosen model for all settings
# ---------------------------------------------------------
evaluate_predictions = function(model_choice, model_list, target, case, income_dt) {
   # Full set of models available - we'll subset for this modal ID
  list[models, covars] = define_models(case)

  # Function to find mode
  mode <- function(x) {
    ux = unique(x[!is.na(x)])
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  # Find most commonly chosen model by region and/or income group
  group_choice_dt = model_choice %>%
                    select(d_v_a_id, country, region, income, model_id, tslm) %>%
                  # Find preferred model by income level
                   group_by(region, income) %>%
                   summarise(mode = mode(model_id)) 
                   
  
  # ---- Summarise predictor coefficients from training countries ---------
  coefficient_dt = model_choice %>%
    tidy() %>%
     lazy_dt() %>%
    # Median coefficient by region and income group to avoid outliers...
    group_by(region, income, model_id, term) %>%
    summarise(estimate = median(estimate, na.rm = TRUE)) %>%
    ungroup() %>%
    # Spread to wide format...
    pivot_wider(
      names_from  = term,
      names_glue  = "{term}_coefficient",
      values_from = estimate) %>%
    as.data.table()
  
  # Allocate chosen model and summarised coefficients to imputed countries
  # Choose model for imputed countries
  impute_choice_dt = model_choice %>%
    full_join(y = table("country"), 
              by =  "country") %>%
    as_tibble() %>%
    select(-c(income,tslm)) %>%
    append_region_name() %>%
    left_join(y = income_dt,
              by = "country") %>%
    fill(d_v_a_id, .direction = "down") %>%
    select(d_v_a_id, country, country_name, region, income, model_id) %>%
    filter(is.na(model_id)) %>%
    left_join(group_choice_dt,
              by = c("region", "income")) %>%
    arrange(region, income) %>%
    group_by(region) %>%
    mutate(model_id = mode) %>%
    # Use model selected for upper-middle income countries for high income countries
    fill(model_id, .direction = "up") %>%
    select(-mode) %>%
    #Choose predictor coefficients for imputed countries
    left_join(coefficient_dt,
              by = c("region", "income", "model_id")) %>%
    arrange(region, income, model_id) %>%
    group_by(region, model_id) %>%
    # Use coefficients for upper-middle income countries for high income countries
    fill(contains("coefficient"), .direction = "up") 
  
  # ---- Full summary of models and predictors for every country ----
  full_coefficient_dt = model_choice %>%
    tidy() %>%
    lazy_dt() %>%
    left_join(y = table("country"), 
     by =  "country") %>%
    mutate(region = region.x) %>%
     select(d_v_a_id, country, country_name, model_id, region, income, term, estimate) %>%
   # Spread to wide format...
    pivot_wider(
      names_from  = term,
      names_glue  = "{term}_coefficient",
      values_from = estimate) %>%
    as.data.table() %>%
    rbind(impute_choice_dt)
  
  # ---- Construct predictor function call ----
  
  # Small function to wrap a string in quotes
  quote = function(x, q = '"') 
    paste0(q, x, q)
  
  # Column names of predictors
  predict_covars = models %>% #[[id]] %>%
    interpret_covars(covars) %>%
    str_remove_all(" ") %>%
    str_split("\\+") %>%
    pluck(1)
  
  # Construct linear product of predictors and coefficients
  predict_str = predict_covars %>%
    paste1("coefficient") %>%
    quote("`") %>%
    paste(predict_covars, sep = " * ") %>%
    paste(collapse = " + ")
  
  # Construct complete function call to be evaluated
  predict_fn   = paste0("prediction = exp(", predict_str, ")")
  predict_eval = paste0("predictors %>% mutate(", predict_fn, ")")
  
  # Append coefficients to predictors
  predictors = target %>% 
    left_join(y = income_dt,
              by = "country") %>%
    left_join(y  = full_coefficient_dt, 
              by = "country") %>% 
    as.data.table()
  
  # Evaluate function call to predict all target values
  predict_dt = eval_str(predict_eval) %>%
    select(country, year, prediction)
  
  return(predict_dt)
}

# ---------------------------------------------------------
# Evaluate given user-specified model
# ---------------------------------------------------------
interpret_covars = function(model, covars) {
  
  # Interpret shorthand references
  for (covar in names(covars))
    model = str_replace_all(
      string  = model, 
      pattern = paste0("\\b", covar, "\\b"), 
      replacement = covars[[covar]])
  
  return(model)
}

# ---------------------------------------------------------
# Easily convert between year and period
# ---------------------------------------------------------
get_period = function() {
  
  # Indices of period change
  year_idx = seq(
    from = o$period_length, 
    to   = length(o$years), 
    by   = o$period_length) + 1
  
  # Format into full year-period datatable
  period_dt = tibble(year = o$years[year_idx]) %>%
    mutate(period = 1 : n()) %>%
    full_join(y  = tibble(year = o$years), 
              by = "year") %>%
    arrange(year) %>%
    fill(period, .direction = "updown") %>%
    as.data.table()
  
  return(period_dt)
}

