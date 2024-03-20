###########################################################
# SIA DATA EXPLORATION
#
# Sources of SIA data:
#   1) https://extranet.who.int/xmart4/WIISE/data/SIA_MAIN
#   2) https://www.who-immunization-repository.org/
#
###########################################################

# ---------------------------------------------------------
# Extract SIA coverage data
# ---------------------------------------------------------
coverage_sia = function(vimc_countries_dt) {
  
  message(" - Coverage data: SIA")
  
  # Throw an error if the necessary data is not available
  sia_file = paste0(o$pth$input, "sia_coverage.csv")
  if (!file.exists(sia_file))
    stop("This pipeline uses one non-publicly available file.\n",
         "Please request this file (sia_coverage.csv) from the developers.\n",
         "Contact: shattocka@who.int")
  
  # ---- Set up ----
  
  # Campaign activities (or 'all' for non-VIMC pathogens)
  d_v_a_dt = table("d_v_a") %>%
    filter(source != "extern") %>%
    bind_rows(table("d_v_a_extern")) %>%
    filter(activity %in% c("campaign", "all")) %>%
    select(d_v_a_id, vaccine)
  
  # Data dictionary for converting to v_a
  data_dict = table("vaccine_dict") %>%
    left_join(y  = d_v_a_dt, 
              by = "vaccine") %>%
    left_join(y  = table("regimen"), 
              by = "vaccine") %>%
    filter(!is.na(d_v_a_id)) %>%
    select(d_v_a_id, vaccine, intervention, schedule) %>%
    arrange(d_v_a_id, intervention)
  
  # ---- Load and wrangle data ----
  
  # Load and wrangle SIA data
  data_dt = fread(paste0(o$pth$input, "sia_coverage.csv")) %>%
    setnames(names(.), tolower(names(.))) %>% 
    # Select columns of interest...
    select(intervention = intervention_code,
           country      = iso3_code,          # Country ISO3 codes, 
           sia_type     = activity_type_code, # Catch up, national day, etc
           status       = activity_status,    # Done, ongoing, planned, etc
           age_group    = activity_age_group,
           cohort       = target, 
           doses, coverage, matches("year$|month$|day$")) %>%
    # Format data types...
    mutate_if(is.character, tolower) %>%
    mutate(across(.cols = c(cohort, coverage, doses), 
                  .fns  = as.numeric)) %>%
    # Attempt to impute doses and then filter out the unworkable...
    mutate(doses = ifelse(is.na(doses), cohort * coverage, doses)) %>%
    filter(doses > 0) %>%
    # Remove any activities not of interest...
    filter(!status   %in% c("cancelled", "postponed"), 
           !sia_type %in% c("routine")) %>%
    select(-status, -sia_type) %>%
    # Remove any unknown countries...
    mutate(country = toupper(country)) %>%
    filter(country %in% all_countries()) %>%
    # Remove any unknown interventions...
    filter(intervention %in% unique(data_dict$intervention)) %>%
    arrange(intervention, country) %>%
    # Deal with dates...
    format_sia_dates() %>%
    impute_sia_dates() %>%
    # Parse age groups...
    parse_age_groups()
  
  # intervention_dt = data_dt %>%
  #   left_join(y  = data_dict,
  #             by = "intervention",
  #             relationship = "many-to-many") %>%
  #   mutate(raw = doses / (schedule * pop), 
  #          coverage = pmin(raw, o$max_coverage)) %>%
  #   select(intervention, vaccine, coverage)
  # 
  # g = ggplot(intervention_dt) +
  #   aes(x = coverage,
  #       colour = intervention,
  #       fill   = intervention) +
  #   geom_histogram(
  #     binwidth = 0.05, 
  #     alpha    = 0.2) +
  #   facet_wrap(
  #     facets = vars(vaccine), 
  #     scales = "free_y")
  
  # Interpret 'interventions'
  sia_dt = data_dt %>%
    # Convert to d-v-a...
    left_join(y  = data_dict, 
              by = "intervention", 
              relationship = "many-to-many") %>%
    filter(!is.na(d_v_a_id)) %>%
    # Remove entires already covered by VIMC...
    left_join(y  = vimc_countries_dt, 
              by = c("d_v_a_id", "country", "year")) %>%
    filter(is.na(source)) %>%
    select(-source) %>%
    # Calculate FVPs...
    mutate(sheduled_doses = doses / schedule) %>%
    calculate_fvps() %>%
    # Tidy up...
    arrange(d_v_a_id, country, year, age) %>%
    mutate(source = "sia") %>%
    as.data.table()
  
  return(sia_dt)
}

# ---------------------------------------------------------
# SIA database has numerous date columns - combine into useable format
# ---------------------------------------------------------
format_sia_dates = function(sia_dt) {
  
  # TODO... do any of the entries with missing start dates have end dates??
  # If so, we could impute the start date as we do for missing end dates.
  
  # Define date columns in raw data set
  date_cols = c("plan", "postponed", "done")
  date_type = c("start", "end")
  
  # Function to create date strings from multiple columns
  date_fn = function(col) {
    
    # Combine all columns to create y-m-d string
    date_str = paste(
      sia_dt[[paste0(col, "_year")]], 
      sia_dt[[paste0(col, "_month")]], 
      sia_dt[[paste0(col, "_day")]], 
      sep = ".")
    
    # Trivialise any dates containing NA
    date_str[grepl("NA", date_str)] = ""
    
    return(date_str)
  }
  
  # Loop through date columns to create a single variable
  for (i in date_cols) {
    for (j in date_type) {
      
      # Construct new column name
      col = paste(i, j, sep = "_")
      
      # Apend this newly compiled column
      sia_dt[[col]] = date_fn(col)
    }
  }
  
  # First and last dates we are interested in
  data_from = paste(min(o$years),     1, 1, sep = ".")
  data_to   = paste(max(o$years + 1), 1, 1, sep = ".")
  
  # Create single start and end columns
  date_dt = sia_dt %>%
    mutate(start_date = paste0(plan_start, postponed_start, done_start), 
           end_date   = paste0(plan_end,   postponed_end,   done_end)) %>%
    # Convert to date format...
    mutate(start_date = format_date(start_date), 
           end_date   = format_date(end_date)) %>%
    # Remove the numerous now-redundant date columns...
    select(-matches("year$|month$|day$"), 
           -ends_with("start"), 
           -ends_with("end")) %>%
    # Remove ineligible dates...
    filter(!is.na(start_date),
           start_date >= format_date(data_from), 
           start_date <  format_date(data_to)) %>%
    arrange(intervention, country, start_date)
  
  return(date_dt)
}

# ---------------------------------------------------------
# Impute missing end dates, distribute over time, and sum over years
# ---------------------------------------------------------
impute_sia_dates = function(sia_dt) {
  
  average_fn = "mean"  # OPTION: "mean" or "median"
  
  # Impute missing end dates
  impute_dt = sia_dt %>%
    # Calculate average duration...
    mutate(days    = as.numeric(end_date - start_date), 
           average = get(average_fn)(days, na.rm = TRUE)) %>%
    # Fill in any missing end dates with duration average...
    mutate(fix_date = start_date + average, 
           end_date = if_else(is.na(end_date), fix_date, end_date)) %>%
    select(-days, -average, -fix_date)
  
  # Melt monthly dates to tidy format...
  
  # Single datatable column of all possible months
  all_months_dt = 
    seq(from = floor_date(min(impute_dt$start_date), "month"), 
        to   = floor_date(max(impute_dt$end_date),   "month"), 
        by   = "month") %>%
    as.character() %>%
    as_named_dt("month")
  
  # All months to distibute doses across (for which campaign has been 'run')
  #
  # NOTE: dtplyr would be handy here, but isn't yet implemented for rowwise() operations
  run_months_dt = impute_dt %>%
    mutate(start_date = floor_date(start_date, "month"),  # Beginning of month
           end_date   = floor_date(end_date,   "month"),  # Beginning of month
           end_date   = pmax(start_date, end_date)) %>%   # In case of end_date < start_date
    rowwise() %>%
    mutate(run_months = seq(start_date, end_date, by = "month") %>%  # All months to distibute across
             paste(collapse = " & ") %>%
             as.character(), 
           n_months = str_count(run_months, "&") + 1) %>%  # Number of months to distibute across
    ungroup() %>%
    as.data.table()
  
  # Expand for all possible and distrubte doses across months
  #
  # NOTES: 
  #  - We'll use this 'all months' dt for pretty plotting
  #  - Whilst this works, there could well be a more efficient way
  sia_month_dt = run_months_dt %>%
    expand_grid(all_months_dt) %>%                     # Full factorial for all possible months
    mutate(value = str_detect(run_months, month)) %>%  # Successful matches
    mutate(month = format_date(month), 
           doses = (doses / n_months) * value) %>%     # Divide total doses across the months
    select(-start_date, -end_date, -run_months, -n_months) %>%
    arrange(intervention, country, month) %>%
    as.data.table()
  
  # Remove these trivial dose entries and sum over year
  sia_year_dt = sia_month_dt %>%
    lazy_dt() %>%
    filter(doses > 0) %>%
    mutate(year = year(month)) %>%
    group_by(intervention, country, year, age_group) %>%
    summarise(doses = sum(doses)) %>%
    ungroup() %>%
    as.data.table()
  
  # Sanity check that we haven't changed number of doses
  # if (abs(sum(sia_year_dt$doses) - check_doses) > 1e-3)
  #   stop("We seem to have gained/lost doses here")
  
  return(sia_year_dt)
}

# ---------------------------------------------------------
# Parse age groups into age ranges
# ---------------------------------------------------------
parse_age_groups = function(sia_dt) {
  
  # All unique age group strings to parse
  age_groups = sort(unique(sia_dt$age_group))
  
  # Initialise working datatable
  group_dt = data.table(
    group = age_groups, 
    age   = age_groups)
  
  # ---- Parse age group stings ----
  
  # Regular expression to remove
  exp_rm = c(",.+", "\n.+", "\\+", "\\=", "[a-z, ]")
  
  # Regular expression to substitute
  exp_sub = c(
    "&"   = "and",
    "^<"  = "0&", 
    "^=<" = "0&",
    "^>"  = "100&", 
    "^>=" = "100&",
    "-<"  = "&", 
    "->"  = "&", 
    "-"   = "&", 
    " y"  = "#", 
    " m"  = "@") 
  
  # Several named special cases
  exp_txt = c(
    "all ages"    = "0-95",
    "adolescents" = "10-19",
    "adults"      = "18-60",
    "children"    = "1-12",
    "elderly"     = "60-95",
    "school"      = "5-16",
    "travellers"  = "18-60",
    "women"       = "18-60")
  
  # Parse key characters
  for (exp in names(exp_sub))
    group_dt[grepl(exp, age), age := gsub(exp, exp_sub[[exp]], age)]
  
  # Remove certain characters
  for (exp in exp_rm)
    group_dt[grepl(exp, age), age := gsub(exp, "", age)]
  
  # Parse strings that we don't yet have info for
  for (exp in names(exp_txt))
    group_dt[grepl(exp, group) & age == "", 
           age := gsub("-", "&", exp_txt[[exp]])]
  
  # ---- Apply parsed ages to data ----
  
  # Construct age group - age range dictionary
  age_dict = group_dt %>%
    mutate(n   = str_count(age, "&"), 
           age = ifelse(n > 1, "", age)) %>%
    # Split at denominator if it exists...
    separate_wider_delim(
      cols  = age, 
      delim = "&", 
      names = c("a1", "a2"), 
      too_few = "align_start", 
      cols_remove = FALSE) %>%
    # Extract units - years (#) or months (@)...
    mutate(u1 = str_extract(a1, "#|@"), 
           u2 = str_extract(a2, "#|@"), 
           u1 = ifelse(is.na(u1), u2, u1), 
           u2 = ifelse(is.na(u2), u1, u2)) %>%
    # Assume units are years if still trivial...
    mutate(u1 = ifelse(is.na(u1), "#", u1), 
           u2 = ifelse(is.na(u2), "#", u2)) %>%
    # Convert ages to numeric...
    mutate(a1 = str_remove(a1, "(#|@).*"), 
           a2 = str_remove(a2, "(#|@).*"), 
           a1 = suppressWarnings(as.numeric(a1)), 
           a2 = suppressWarnings(as.numeric(a2))) %>%
    # Convert all ages to year format...
    mutate(y1 = ifelse(u1 == "@", a1 / 12, a1), 
           y2 = ifelse(u2 == "@", a2 / 12, a2)) %>%
    # Deal with missing and infeasible values...
    mutate(y1 = ifelse(y1 > 1e3, NA, y1), 
           y1 = pmin(y1, max(o$ages)), 
           y2 = ifelse(is.na(y2), y1, y2)) %>%
    replace_na(list(y1 = 0, y2 = 0)) %>%
    # Ensure correct order and round...
    mutate(age_min = round(pmin(y1, y2)), 
           age_max = round(pmax(y1, y2))) %>%
    # Expand out for single age bins...
    select(age_group = group, age_min, age_max) %>%
    expand_grid(age = o$ages) %>%
    filter(age >= age_min, 
           age <= age_max) %>%
    # Tidy up...
    select(age_group, age) %>%
    as.data.table()
  
  # Convert to long form and distribute across age bins
  age_dt = sia_dt %>%
    lazy_dt() %>%
    rename(total_doses = doses) %>%
    # Expand with parsed age bins...
    left_join(y  = age_dict, 
              by = "age_group", 
              relationship = "many-to-many") %>%
    # Append population size...
    left_join(y  = table("wpp_pop"), 
              by = c("country", "year", "age")) %>%
    # Distribute across age bins...
    group_by(intervention, country, year, age_group) %>%
    mutate(doses = total_doses * pop / sum(pop)) %>%
    ungroup() %>%
    # Tidy up...
    select(intervention, country, year, age, doses, pop) %>%
    as.data.table()
  
  # Sanity check that we haven't lost/gained doses
  dose_diff = sum(sia_dt$doses) - sum(age_dt$doses)
  if (abs(dose_diff) > 1e-6)
    stop("Age disaggregation failed")
    
  return(age_dt)
}

