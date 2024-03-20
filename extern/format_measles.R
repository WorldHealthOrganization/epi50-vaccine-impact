
# ---- Load raw results ----

vaccine_dt    = fread(paste0(o$pth$extern, "epi50_measles_vaccine.csv"))
no_vaccine_dt = fread(paste0(o$pth$extern, "epi50_measles_no_vaccine.csv"))

raw_dt = vaccine_dt %>%
  select(iso, year, 
         vaccine.deaths = mean_deaths, 
         vaccine.cases  = mean_cases) %>%
  left_join(y  = no_vaccine_dt, 
            by = c("iso", "year")) %>%
  rename(country = iso, 
         no_vaccine.deaths = mean_deaths, 
         no_vaccine.cases  = mean_cases) %>%
  select(country, year, 
         starts_with("vaccine"), 
         starts_with("no_vaccine")) %>%
  filter(!is.na(vaccine.deaths), 
         !is.na(no_vaccine.deaths)) %>%
  mutate(vaccine.deaths = pmin(no_vaccine.deaths, vaccine.deaths), 
         vaccine.cases  = pmin(no_vaccine.cases,  vaccine.cases)) %>%
  pivot_longer(cols = -c(country, year)) %>%
  separate(col  = "name", 
           into = c("scenario", "metric"), 
           sep  = "\\.") %>%
  select(scenario, country, year, metric, value) %>%
  arrange(scenario, country, year, metric) %>%
  as.data.table()
  
# ---- Age structure ----

epi50_dynamice = read_rds("extern", "epi50_dynamice_results")

age_dt = epi50_dynamice %>%
  filter(metric == "deaths") %>%
  pivot_wider(names_from = scenario) %>%
  mutate(value = pmax(no_vaccine - vaccine, 0)) %>%
  group_by(country, year) %>%
  mutate(scaler = value / sum(value)) %>%
  ungroup() %>%
  select(country, year, age, scaler) %>%
  group_by(year, age) %>%
  mutate(impute = mean(scaler, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(scaler = ifelse(is.na(scaler), impute, scaler)) %>%
  select(-impute) %>%
  as.data.table()

# ---- Calculate DALYs ----

disability_weight = 0.002

life_exp_dt = table("wpp_life_exp") %>%
  select(-age)

dalys_dt = raw_dt %>%
  left_join(y  = age_dt, 
            by = c("country", "year"), 
            relationship = "many-to-many") %>%
  mutate(value = value * scaler) %>%
  select(all_names(epi50_dynamice)) %>%
  pivot_wider(names_from = metric) %>%
  left_join(y  = life_exp_dt, 
            by = c("country", "year")) %>%
  mutate(yll = deaths * pmax(life_exp - age, 0)) %>%
  mutate(yld = disability_weight * (cases - deaths)) %>%
  mutate(dalys = yll + yld) %>%
  select(scenario, country, year, age, dalys) %>%
  pivot_longer(cols = dalys,
               names_to = "metric") %>%
  select(all_names(epi50_dynamice)) %>%
  as.data.table()

# ---- Impute missing values ----

measles_dt = raw_dt %>%
  filter(metric == "deaths") %>%
  left_join(y  = age_dt, 
            by = c("country", "year"), 
            relationship = "many-to-many") %>%
  mutate(value = value * scaler) %>%
  select(all_names(epi50_dynamice)) %>%
  rbind(dalys_dt)

impute_dt = epi50_dynamice %>%
  rename(dynamice = value) %>%
  left_join(y  = measles_dt, 
            by = intersect(names(.), names(measles_dt))) %>%
  mutate(value = ifelse(
    test = is.na(value) & year > min(o$year), 
    yes  = dynamice, 
    no   = value))

init_dt = impute_dt %>%
  filter(!is.na(value), 
         year <= min(year) + 3) %>%
  group_by(scenario, country, age, metric) %>%
  summarise(value = mean(value)) %>%
  ungroup() %>%
  mutate(year = min(o$year)) %>%
  select(all_names(measles_dt)) %>%
  as.data.table()

# ---- Concatenate ----

epi50_measles = impute_dt %>%
  select(-dynamice) %>%
  filter(!is.na(value)) %>%
  rbind(init_dt) %>%
  arrange(scenario, country, year, age, metric)

save_rds(epi50_measles, "extern", "epi50_measles_results")

# ---- Plot outcomes ----

total_dt = epi50_measles %>%
  pivot_wider(names_from = scenario) %>% 
  mutate(averted = no_vaccine - vaccine) %>%
  group_by(metric, year) %>%
  summarise(averted = sum(averted)) %>%
  mutate(cum_averted = cumsum(averted)) %>%
  ungroup() %>%
  as.data.table()

g1 = ggplot(total_dt) + 
  aes(x = year, 
      y = averted) + 
  facet_wrap(
    facets = vars(metric), 
    scales = "free_y") + 
  geom_line()

g2 = ggplot(total_dt) + 
  aes(x = year, 
      y = cum_averted) + 
  facet_wrap(
    facets = vars(metric), 
    scales = "free_y") + 
  geom_line()

