###########################################################
# PLOTTING
#
# All plotting functionality in one place.
#
###########################################################

# **** Data visualisation ****

# ---------------------------------------------------------
# Plot methodology figure to be used in paper
# ---------------------------------------------------------
plot_scope = function() {
  
  message("  - Plotting country-disease scope")
  
  # Manually set tidy y axis limit
  y_max = 8  # In billions
  
  # Linear time interpolation for less pixilated figures
  smoothness = 10  # Higher value for smoother plot
  
  # Dictionary for full impact source descriptions
  impact_dict = c(
    extern = "Dynamic modelling (external to VIMC)",
    vimc   = "Dynamic modelling (contributing to VIMC)", 
    static = "Static modelling", 
    impute = "Geographic imputation model", 
    extrap = "Temporal extrapolation model")
  
  # Associated colours
  impact_colours = c("#EB7D5B", "#FED23F", "#B5D33D", "#6CA2EA", "#442288")
  
  # ---- Number of FVPs by pathogen ----
  
  # Number of FVPs by country
  fvps_dt = table("coverage") %>%
    filter(coverage > 0) %>%
    # Append disease details...
    left_join(y  = table("d_v_a"), 
              by = "d_v_a_id") %>%
    # Concept here is FVP, so remove birth dose & boosters...
    filter(!grepl("_(bd|bx|px)$", vaccine)) %>% 
    # Summarise over age...
    group_by(disease, country, year) %>%
    summarise(fvps = sum(fvps)) %>%
    ungroup() %>%
    as.data.table()
  
  # Total cumulative FVPs per disease
  total_dt = fvps_dt %>%
    # Total by disease...
    group_by(disease, year) %>%
    summarise(total = sum(fvps)) %>%
    ungroup() %>%
    # Cumulative over time...
    group_by(disease) %>%
    mutate(cum_total = cumsum(total)) %>%
    ungroup() %>%
    as.data.table()
  
  # ---- Source of impact estimates ----
  
  # Static model settings
  static_dt = table("gbd_estimates") %>%
    select(disease, country, year) %>%
    unique() %>%
    # filter(year %in% o$gbd_prediction_years) %>%
    mutate(class = "static")
  
  # External model settings
  extern_dt = table("extern_estimates") %>%
    left_join(y  = table("d_v_a"), 
              by = "d_v_a_id") %>%
    select(disease, country, year) %>%
    unique() %>%
    mutate(class = "extern")
  
  # VIMC settings
  vimc_dt = table("vimc_estimates") %>%
    left_join(y  = table("d_v_a"), 
              by = "d_v_a_id") %>%
    select(disease, country, year) %>%
    unique() %>%
    mutate(class = "vimc")
  
  # VIMC country imputation settings
  impute_dt = table("d_v_a") %>%
    filter(source == "vimc") %>%
    select(disease) %>%
    unique() %>%
    expand_grid(
      country = all_countries(), 
      year    = unique(vimc_dt$year)) %>%
    left_join(y  = vimc_dt, 
              by = c("disease", "country", "year")) %>%
    filter(is.na(class)) %>%
    mutate(class = "impute") %>%
    as.data.table()
  
  # ---- Construct plotting datatable ----
  
  # Combine all impact sources
  plot_dt = rbind(extern_dt, static_dt, vimc_dt, impute_dt) %>%
    # Append FVPs...
    right_join(y  = fvps_dt, 
               by = c("disease", "country", "year")) %>%
    # Anything not yet specified is time-exrapolated...
    mutate(class = ifelse(is.na(class), "extrap", class)) %>%
    # Summarise over all countries...
    group_by(disease, class, year) %>%
    summarise(fvps = sum(fvps)) %>%
    ungroup() %>%
    # Append cumulative total...
    left_join(y  = total_dt, 
              by = c("disease", "year")) %>%
    # Share of cumulative FVPs by class over time...
    mutate(share = fvps / total, 
           value = share * cum_total / 1e9) %>%
    # Expand for a more granular timescale...
    expand_grid(time = seq(
      from = 0 , 
      to   = 1 - 1 / smoothness, 
      by   = 1 / smoothness)) %>%
    mutate(time = year + time) %>%
    filter(time <= max(o$years)) %>%
    # Interpolate annual values for smoother plot...
    mutate(value = ifelse(time == year, value, NA)) %>%
    group_by(disease, class) %>%
    mutate(value = na_interpolation(value)) %>%
    ungroup() %>%
    # Use full disease names...
    left_join(y  = table("disease_name"), 
              by = "disease") %>%
    select(disease = disease_name, class, time, value) %>%
    # Use impact source descriptions...
    mutate(class = recode(class, !!!impact_dict), 
           class = factor(class, rev(impact_dict))) %>%
    # Set disease order...
    arrange(desc(value)) %>%
    mutate(disease = fct_inorder(disease)) %>%
    # Tidy up...
    arrange(disease, class, time) %>%
    as.data.table()
  
  # ---- Construct label datatable ----
  
  # Year range of analysis
  year1 = min(o$years)
  year2 = max(o$years)
  
  # Label description string
  label = paste0("Total (", year1, "-", year2, "):")
  
  # Construct labels: total FVPs over analysis timeframe
  label_dt = plot_dt %>%
    filter(time == year2) %>%
    # Total cumulative FVPs by disease...
    group_by(disease) %>%
    summarise(total = sum(value)) %>%
    ungroup() %>%
    # Construct labels...
    mutate(total = round(total, 2),
           label = paste(label, total, "billion")) %>%
    # Set coordinates...
    mutate(time  = year1 + 0.01 * (year2 - year1),
           value = y_max * 0.9) %>%
    as.data.table()
  
  # ---- Produce plot ----
  
  # Plot FVP over time per pathogen and impact source
  g = ggplot(plot_dt) +
    aes(x = time, 
        y = value) +
    geom_bar(
      mapping  = aes(fill = class), 
      stat     = "identity", 
      position = "stack",  
      width    = 1 / smoothness) +
    # Add total labels...
    geom_text(
      data    = label_dt, 
      mapping = aes(label = label), 
      size    = 3.5,
      hjust   = 0, 
      vjust   = 1) + 
    # Some intricate faceting...
    facet_rep_wrap(
      facets = vars(disease), 
      ncol   = 1, 
      labeller = label_wrap_gen(width = 20), 
      strip.position = "right", 
      repeat.tick.labels = FALSE) + 
    # Set colours and legend title...
    scale_fill_manual(
      values = FALSE, 
      name   = "Source of impact predictions") +
    guides(fill = guide_legend(reverse = TRUE)) +
    # Prettify x axis...
    scale_x_continuous(
      limits = c(year1 - 1 / smoothness, 
                 year2 + 1 / smoothness), 
      expand = expansion(mult = c(0, 0)), 
      breaks = seq(year1, year2, by = 5)) +  
    # Prettify y axis...
    scale_y_continuous(
      name   = paste("Cumulative number of people receiving",
                     "final primary dose (in billions)"), 
      limits = c(0, y_max),
      breaks = seq(2, y_max, by = 2),
      labels = comma,
      expand = expansion(mult = c(0, 0)))
  
  browser()
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(strip.text    = element_text(size = 14),
          strip.text.y  = element_text(angle = 0, hjust = 0),
          strip.background = element_blank(), 
          axis.title.x  = element_blank(),
          axis.title.y  = element_text(
            size = 18, margin = margin(l = 10, r = 20)),
          axis.text     = element_text(size = 8),
          axis.ticks    = element_blank(), 
          axis.line     = element_line(linewidth = 0.25),
          panel.spacing = unit(-0.6, "lines"), 
          panel.grid.major.y = element_line(linewidth = 0.25),
          legend.title  = element_text(size = 14),
          legend.text   = element_text(size = 11),
          legend.key    = element_blank(),
          legend.position = "right", 
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines"))
  
  # Save to file
  save_fig(g, "Analysis scope", dir = "data_visualisation")
  save_fig(g, "Figure S1", dir = "manuscript")
}

# ---------------------------------------------------------
# Plot total number of FVP over time for each d_v_a
# ---------------------------------------------------------
plot_total_fvps = function() {
  
  message("  - Plotting total number of FVP")
  
  # Flag for whether to plot FVPs cumulatively over time
  cumulative = TRUE
  
  # Number of FVPs by source of data
  source_dt = table("coverage_source") %>%
    # Summarise over countries and age...
    group_by(d_v_a_id, source, year) %>%
    summarise(fvps = sum(fvps) / 1e9) %>%
    ungroup() %>%
    # Cumulative FVPs...
    group_by(d_v_a_id, source) %>%
    mutate(fvps_cum = cumsum(fvps)) %>%
    ungroup() %>%
    # Tidy up...
    format_d_v_a_name() %>%
    filter(!is.na(d_v_a_name)) %>%
    arrange(d_v_a_name, source, year) %>%
    as.data.table()
  
  # Total FVPs (sum of all sources)
  #
  # NOTE: Not necessarily equal to sum of all sources as
  #       SIA are assumed to be only partially targeted
  total_dt = table("coverage") %>%
    # Summarise over countries and age...
    group_by(d_v_a_id, year) %>%
    summarise(fvps = sum(fvps) / 1e9) %>%
    ungroup() %>%
    # Cumulative FVPs...
    group_by(d_v_a_id) %>%
    mutate(fvps_cum = cumsum(fvps)) %>%
    ungroup() %>%
    # Tidy up...
    format_d_v_a_name() %>%
    arrange(d_v_a_id, year) %>%
    as.data.table()
  
  # Metric to use for y axis
  y = ifelse(cumulative, "fvps_cum", "fvps")
  
  # Plot FVPs over time for each d_v_a
  g = ggplot(source_dt) + 
    aes(x = year, y = !!sym(y)) + 
    geom_line(
      mapping   = aes(colour = source), 
      linewidth = 1.5) + 
    geom_line(
      data      = total_dt, 
      linetype  = "dashed",
      colour    = "black", 
      linewidth = 1.5) + 
    # Facet with strip text wrapping...
    facet_wrap(
      facets   = vars(d_v_a_name), 
      labeller = label_wrap_gen(width = 24), 
      scales   = "free_y") + 
    # Prettify x axis...
    scale_x_continuous(
      limits = c(min(o$years), max(o$years)), 
      expand = expansion(mult = c(0, 0)), 
      breaks = seq(min(o$years), max(o$years), by = 10)) +  
    # Prettify y axis...
    scale_y_continuous(
      name   = "Total receiving full schedule (in billions)", 
      labels = comma,
      expand = expansion(mult = c(0, NA)))
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.title.x  = element_blank(),
          axis.title.y  = element_text(
            size = 20, margin = margin(l = 10, r = 20)),
          axis.text     = element_text(size = 10),
          axis.text.x   = element_text(hjust = 1, angle = 50), 
          axis.line     = element_blank(),
          strip.text    = element_text(size = 12),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(1, "lines"),
          panel.grid.major.y = element_line(linewidth = 0.5),
          legend.title  = element_blank(),
          legend.text   = element_text(size = 14),
          legend.key    = element_blank(),
          legend.position = "right", 
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines"))
  
  # Save to file
  save_dir = "data_visualisation"
  save_fig(g, "FVPs by source", dir = save_dir)
}

# ---------------------------------------------------------
# Plot smoothed FVP for static model pathogens
# ---------------------------------------------------------
plot_smooth_fvps = function() {
  
  message("  - Plotting smoothed FVPs (static model pathogens)")
  
  # ---- Plot 1: data vs smoothing ----
  
  # Smoothing data - group by country and age
  smooth_dt = table("smoothed_fvps") %>%
    mutate(country_age = paste1(country, age), 
           fvps_smooth = fvps_smooth / 1e6, 
           fvps        = fvps / 1e6) %>%
    format_d_v_a_name() %>%
    select(d_v_a_name, country, country_age, 
           year, fvps, fvps_smooth)
  
  # Plot the data with associated smoothing
  g1 = ggplot(smooth_dt) +
    aes(x = year, 
        colour = country, 
        group  = country_age) +
    # Plot data...
    geom_point(
      mapping = aes(y = fvps),
      show.legend = FALSE) +
    # Plot smoothing...
    geom_line(
      mapping = aes(y = fvps_smooth),
      show.legend = FALSE) +
    # Simple faceting...
    facet_wrap(~d_v_a_name) + 
    # Prettify x axis...
    scale_x_continuous(
      expand = expansion(mult = c(0, 0)), 
      breaks = pretty_breaks()) +  
    # Prettify y axis...
    scale_y_continuous(
      name   = "Fully vaccinated people (in millions)", 
      labels = comma,
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05)))
  
  # Prettify theme
  g1 = g1 + theme_classic() + 
    theme(axis.title.x  = element_blank(),
          axis.title.y  = element_text(
            size = 20, margin = margin(l = 10, r = 20)),
          axis.text     = element_text(size = 10),
          # axis.text.x   = element_text(hjust = 1, angle = 50), 
          axis.line     = element_blank(),
          strip.text    = element_text(size = 14),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(1, "lines"))
  
  # ---- Plot 2: smoothing error ----
  
  # Total error for each vaccine
  diagnostic_dt = smooth_dt %>%
    # Summarise for each vaccine...
    group_by(d_v_a_name) %>%
    summarise(fvps = sum(fvps),
              fvps_smooth = sum(fvps_smooth)) %>%
    ungroup() %>%
    # Take the absolute difference...
    mutate(diff = fvps - fvps_smooth,
           abs  = abs(diff)) %>%
    # Calculate total error...
    mutate(err = round(100 * abs / fvps, 2), 
           err = paste0("Error: ", err, "%")) %>%
    # Append error to vaccine description...
    mutate(d_v_a_name = paste0(d_v_a_name, "\n", err)) %>%
    # Set plotting order by abs diff...
    arrange(abs) %>%
    mutate(d_v_a_name = fct_inorder(d_v_a_name)) %>%
    as.data.table()
  
  # Plot total smoothing errors
  g2 = ggplot(diagnostic_dt) +
    aes(x = d_v_a_name, 
        y = diff, 
        fill = abs) +
    geom_col(show.legend = FALSE) +
    coord_flip() + 
    # Prettify x axis (noting coord_flip)...
    scale_y_continuous(
      name   = "Total difference in FVPs after smoothing (in millions)", 
      labels = comma)
  
  # Prettify theme
  g2 = g2 + theme_classic() + 
    theme(axis.text     = element_text(size = 12),
          axis.title.x  = element_text(
            size = 20, margin = margin(b = 10, t = 20)),
          axis.title.y  = element_blank(),
          axis.line     = element_blank(),
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA))
  
  # ---- Save diagnostic plots ----
  
  # Directory to save to
  save_dir = "data_visualisation"
  
  # Save figures to file
  save_fig(g1, "Data smoothing outcomes",   dir = save_dir)
  save_fig(g2, "Data smoothing difference", dir = save_dir)
}

# ---------------------------------------------------------
# Plot coverage density by disease
# ---------------------------------------------------------
plot_coverage = function() {
  
  # All coverage values by disease, vaccine, and age class
  all_coverage_dt = table("coverage") %>%
    # Classify by age group...
    mutate(age_group = ifelse(age == 0, "infant", "other")) %>%
    # Classify by disease...
    left_join(y  = table("d_v_a"), 
              by = "d_v_a_id") %>%
    # Tidy up...
    select(disease, vaccine, age_group, coverage) %>%
    arrange(disease, vaccine, age_group)
  
  # Details for file destination
  save_name = "Coverage value density"
  save_dir  = "data_visualisation"
  
  # Produce plot for both disease and vaccine
  for (d_v in c("disease", "vaccine")) {
    
    # Select appropriate column: disease or vaccine
    plot_dt = all_coverage_dt %>%
      select(d_v = !!d_v, age_group, coverage)
    
    # Plot coverage value density
    g = ggplot(plot_dt) + 
      aes(x = coverage, 
          y = after_stat(scaled), 
          colour = age_group,
          fill   = age_group) + 
      geom_density(alpha = 0.2) + 
      facet_wrap(~d_v) + 
      # Prettify x axis...
      scale_x_continuous(
        name   = "Coverage",
        labels = percent,
        limits = c(0, 1), 
        expand = expansion(mult = c(0, 0)), 
        breaks = pretty_breaks()) +  
      # Prettify y axis...
      scale_y_continuous(
        name   = "Density", 
        limits = c(0, 1), 
        expand = expansion(mult = c(0, 0.1)), 
        breaks = pretty_breaks())
    
    # Prettify theme
    g = g + theme_classic() + 
      theme(axis.text     = element_text(size = 10),
            axis.title.x  = element_text(
              size = 20, margin = margin(b = 10, t = 20)),
            axis.title.y  = element_text(
              size = 20, margin = margin(l = 10, r = 20)),
            axis.line     = element_blank(),
            strip.text    = element_text(size = 14),
            strip.background = element_blank(), 
            panel.border  = element_rect(
              linewidth = 0.5, fill = NA),
            panel.spacing = unit(1, "lines"),
            legend.title  = element_blank(),
            legend.text   = element_text(size = 14),
            legend.key    = element_blank(),
            legend.position = "right", 
            legend.key.height = unit(2, "lines"),
            legend.key.width  = unit(2, "lines"))
    
    # Save to file
    save_fig(g, save_name, d_v, dir = save_dir)
  }
}

# ---------------------------------------------------------
# Plot age targets as defined by WIISE and VIMC coverage data
# ---------------------------------------------------------
plot_coverage_age_density = function() {
  
  message("  - Plotting coverage data density by age")
  
  # Construct plotting datatable
  plot_dt = table("coverage_source") %>%
    mutate(trans_age = pmax(age, 1), .after = age) %>%
    format_d_v_a_name() %>%
    filter(!is.na(d_v_a_name))
  
  # Plot age density of coverage data by source
  g = ggplot(plot_dt) +
    aes(x = trans_age, 
        y = after_stat(scaled), 
        colour = source,
        fill   = source) +
    geom_density(alpha = 0.2) +
    # Facet with strip text wrapping...
    facet_wrap(
      facets   = vars(d_v_a_name), 
      labeller = label_wrap_gen(width = 24)) + 
    # Prettify x axis...
    scale_x_continuous(
      name   = "Age (log2 scale)",
      trans  = "log2", 
      limits = c(1, 2^7), 
      expand = c(0, 0), 
      breaks = 2^(0:7)) +  
    # Prettify y axis...
    scale_y_continuous(
      name   = "Density", 
      limits = c(0, 1), 
      expand = expansion(mult = c(0, 0.1)), 
      breaks = pretty_breaks())
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.text     = element_text(size = 10),
          axis.title.x  = element_text(
            size = 20, margin = margin(b = 10, t = 20)),
          axis.title.y  = element_text(
            size = 20, margin = margin(l = 10, r = 20)),
          axis.line     = element_blank(),
          strip.text    = element_text(size = 12),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(1, "lines"),
          legend.title  = element_blank(),
          legend.text   = element_text(size = 12),
          legend.key    = element_blank(),
          legend.position = "right", 
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines"))
  
  # Save to file
  save_fig(g, "Coverage density by age", dir = "data_visualisation")
}

# ---------------------------------------------------------
# Plot countries with missing coverage data
# ---------------------------------------------------------
plot_missing_data = function() {
  
  message("  - Plotting countries with missing coverage data")
  
  # Country population (most recent year)
  pop_dt = table("wpp_pop") %>%
    filter(year == max(o$years)) %>%
    group_by(country) %>%
    summarise(pop = sum(pop)) %>%
    ungroup() %>%
    mutate(pop_norm = sqrt(pop / max(pop))) %>%
    as.data.table()
  
  # Identify countries with no data
  plot_dt = table("coverage") %>%
    left_join(y  = table("d_v_a"), 
              by = "d_v_a_id") %>%
    select(disease, country) %>%
    unique() %>%
    mutate(value = 0) %>%
    # Wrangle to sparse, long format...
    pivot_wider(id_cols     = "country", 
                names_from  = "disease") %>%
    pivot_longer(cols = -country, 
                 names_to = "disease") %>%
    # Append pop size for missing countries...
    left_join(y  = pop_dt, 
              by = "country") %>%
    mutate(value = ifelse(is.na(value), pop_norm, NA)) %>%
    # Append region details...
    left_join(y  = table("country"), 
              by = "country") %>%
    select(disease, region, country = country_name, value) %>%
    arrange(disease, region, country) %>%
    # Set plotting order...
    mutate(region  = fct_inorder(region),
           country = fct_inorder(country)) %>% 
    as.data.table()
  
  # Extract population limits for the data
  limits  = c(0, max(plot_dt$value, na.rm = TRUE))
  colours = colour_scheme("pals::brewer.reds", n = 10)
  
  # Plot missing data
  g = ggplot(plot_dt) +
    aes(x = disease, 
        y = fct_rev(country), 
        fill = value) + 
    # Plot missing countries in colour...
    geom_tile() + 
    # Facet by region for readability...
    facet_wrap(
      facets = vars(region), 
      scales = "free") + 
    # Prettify y axis...
    scale_y_discrete(
      name = "Countries with missing data") +
    # Set continuous colour bar...
    scale_fill_gradientn(
      na.value = "white",
      colours  = colours,
      limits   = limits, 
      breaks   = limits, 
      labels   = c("Smaller population", 
                   "Larger population"))
  
  # Prettify theme
  g = g + theme_classic() +
    theme(axis.title.x  = element_blank(),
          axis.title.y  = element_text(size = 28),
          axis.text.x   = element_text(
            size = 8, hjust = 1, angle = 50),
          axis.text.y   = element_text(size = 8),
          axis.line     = element_blank(),
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing.x = unit(0.5, "lines"), 
          panel.spacing.y = unit(0.5, "lines"), 
          strip.text    = element_text(size = 12),
          strip.background = element_blank(), 
          legend.title  = element_blank(),
          legend.text   = element_text(size = 10),
          legend.position = "bottom", 
          legend.key.width = unit(4, "cm"))
  
  # Save figure to file
  save_fig(g, "Missing data by country", 
           dir = "data_visualisation")
}

# **** Static models ****

# ---------------------------------------------------------
# Plot Global Burden of Disease burden estimates by age
# ---------------------------------------------------------
plot_gbd_estimates = function() {
  
  message("  - Plotting GBD burden estimates by age")
  
  # Flag for plotting recent extrapolation
  plot_extrap = FALSE
  
  # Define age groups and associated upper bounds
  age_groups = c(
    "Neonates"      = -1,
    "Other infants" = 0,
    "1-4 years"     = 5,
    "5-14"          = 15,
    "15-49"         = 50,
    "50-69"         = 70,
    "70+ years"     = max(o$ages))
  
  # Map each age bin to respective age group
  age_group_dt = data.table(age = c(-1, o$ages)) %>%
    mutate(group_idx = match(age, age_groups),
           group = names(age_groups[group_idx])) %>%
    fill(group, .direction = "up") %>%
    select(age, age_group = group)
  
  # Plot up to this year
  plot_to = ifelse(plot_extrap, max(o$years), max(o$gbd_estimate_years))
  
  # Load GBD predictions and categorise into age groups
  gbd_dt = table("gbd_estimates") %>%
    filter(year <= plot_to) %>%
    append_d_v_t_name() %>%
    # Append age group details...
    left_join(y  = age_group_dt,
              by = "age") %>%
    # Summarise for these broad age groups...
    lazy_dt() %>%
    group_by(disease, year, age_group) %>%
    summarise(deaths = sum(deaths_disease), 
              dalys  = sum(dalys_disease)) %>%
    ungroup() %>%
    # Melt to long format...
    pivot_longer(cols = c(deaths, dalys), 
                 names_to = "metric") %>%
    replace_na(list(value = 0)) %>%
    append_metric_name() %>%
    # Set age group factors for meaningful plotting order...
    mutate(age_group = factor(age_group, names(age_groups))) %>%
    select(disease, metric, year, age_group, value) %>%
    arrange(disease, metric, year, age_group) %>%
    as.data.table()
  
  # Plot deaths over time by age group
  g = ggplot(gbd_dt) +
    aes(x = year, 
        y = value, 
        fill = age_group) +
    geom_bar(stat = "identity") +
    # Facet by disease...
    facet_grid2(
      cols   = vars(disease), 
      rows   = vars(metric), 
      scales = "free_y", 
      independent = "y") + 
    # Set colour scheme...
    scale_fill_manual(
      name   = "Age group",
      values = colour_scheme(
        map = "brewer::paired", 
        n   = length(age_groups))) +
    # Prettify y axis...
    scale_y_continuous(
      labels = comma,
      expand = expansion(mult = c(0, 0.05)), 
      breaks = pretty_breaks())
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.title    = element_blank(),
          axis.text     = element_text(size = 10),
          axis.line     = element_blank(),
          strip.text    = element_text(size = 14),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(1, "lines"),
          legend.title  = element_blank(),
          legend.text   = element_text(size = 14),
          legend.key    = element_blank(),
          legend.position = "bottom", 
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines"))
  
  # Save to file
  save_fig(g, "GBD disease burden by age group", 
           dir = "static_models")
}

# ---------------------------------------------------------
# Plot proportion of GBD burden we have coverage data for
# ---------------------------------------------------------
plot_gbd_missing = function() {
  
  message("  - Plotting GBD burden by vaccine coverage status")
  
  # Dictionary for status groups
  status_dict = c(
    trivial     = "Countries with NO coverage data",
    non_trivial = "Countries with coverage data")
  
  # Countries which have non-trivial status
  coverage_dt = table("coverage") %>%
    left_join(y  = table("d_v_a"), 
              by = "d_v_a_id") %>%
    filter(source == "static") %>%
    select(disease, country) %>%
    unique() %>%
    mutate(status = "non_trivial")
  
  # GBD disease burden by status group
  status_dt = table("gbd_estimates") %>%
    group_by(disease, country, year) %>%
    summarise(deaths = sum(deaths_disease)) %>%
    ungroup() %>%
    left_join(y  = coverage_dt, 
              by = c("disease", "country")) %>%
    replace_na(list(status = "trivial")) %>%
    left_join(y  = table("disease_name"), 
              by = "disease") %>%
    group_by(disease_name, year, status) %>%
    summarise(deaths = sum(deaths)) %>%
    ungroup() %>%
    mutate(status = recode(status, !!!status_dict), 
           status = factor(status, status_dict)) %>%
    as.data.table()
  
  # Plot burden over time by vaccine coverage status
  g = ggplot(status_dt) +
    aes(x = year, 
        y = deaths, 
        fill = status) +
    geom_bar(stat = "identity") +
    # Facet by disease...
    facet_wrap(
      facets = vars(disease_name), 
      scales = "free_y") + 
    # Set colour scheme...
    scale_fill_manual(
      name   = "Age group",
      values = c("red", "gray")) +
    # Prettify y axis...
    scale_y_continuous(
      name   = "GBD-estimated number of deaths", 
      labels = comma,
      expand = expansion(mult = c(0, 0.05)), 
      breaks = pretty_breaks())
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.title.x  = element_blank(),
          axis.title.y  = element_text(
            size = 20, margin = margin(l = 10, r = 20)),
          axis.text     = element_text(size = 10),
          axis.line     = element_blank(),
          strip.text    = element_text(size = 14),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(1, "lines"),
          legend.title  = element_blank(),
          legend.text   = element_text(size = 14),
          legend.key    = element_blank(),
          legend.position = "bottom", 
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines"))
  
  # Save to file
  save_fig(g, "GBD deaths by vaccine coverage status", 
           dir = "static_models")
}

# ---------------------------------------------------------
# Plot static model pathogen vaccine efficacy with immunity decay
# ---------------------------------------------------------
plot_vaccine_efficacy = function() {
  
  message("  - Plotting vaccine efficacy profiles")
  
  schedule_dict = c(
    x  = "Primary schedule",
    bx = "Booster schedule", 
    px = "Pregnancy schedule")
  
  # Function to extract vaccine efficacy data
  extract_data_fn = function(dt) {
    
    # Extract and format into datatable
    data = dt$data %>%
      unlist() %>%
      as_named_dt("value") %>%
      mutate(name    = c("y", "x"), 
             vaccine = dt$vaccine) %>%
      pivot_wider(id_cols = vaccine) %>%
      as.data.table()
    
    return(data)
  }
  
  # Function to group similar vaccines but split by dose
  schedule_fn = function(dt) {
    
    # Append descriptive columns
    shedule_dt = dt %>%
      # Primary schedule or booster dose...
      separate(col  = vaccine, 
               into = c("type", "schedule"), 
               sep  = "_", 
               fill = "right", 
               remove = FALSE) %>%
      replace_na(list(schedule = "x")) %>%
      mutate(schedule = recode(schedule, !!!schedule_dict), 
             schedule = factor(schedule, schedule_dict)) %>%
      # Append vaccine type name
      mutate(type = str_remove(type, "[0-9]")) %>%
      left_join(y  = table("type_name"), 
                by = "type") %>%
      select(-type) %>%
      rename(type = type_name) %>%
      mutate(type = fct_inorder(type)) %>%
      as.data.table()
    
    return(shedule_dt)
  }
  
  # Extract all vaccine efficacy data points
  data_dt = table("vaccine_efficacy") %>%
    dtapply(extract_data_fn) %>%
    rbindlist() %>%
    inner_join(y  = table("d_v_a"), 
               by = "vaccine") %>%
    schedule_fn() %>%
    select(type, schedule, time = x, value = y)
  
  # Load vaccine efficacy profiles
  profile_dt = table("vaccine_efficacy_profiles") %>%
    schedule_fn() %>%
    select(type, schedule, time, value = profile)
  
  # Plot vaccine efficacy with waning immunity (if any)
  g = ggplot(profile_dt) + 
    aes(x = time, 
        y = value) + 
    # Plot immunity profile...
    geom_line(
      colour    = "black",
      linewidth = 1.5) + 
    # Plot coloured data points...
    geom_point(
      data    = data_dt, 
      mapping = aes(colour = type),
      size    = 2.5) + 
    # Facet by vaccine and schedule...
    facet_grid(
      rows     = vars(type), 
      cols     = vars(schedule),
      labeller = label_wrap_gen(width = 20)) + 
    # Prettify x axis...
    scale_x_continuous(
      name   = "Years after completion of full schedule",
      expand = expansion(mult = c(0.05, 0.05)), 
      breaks = pretty_breaks()) +
    # Prettify y axis...
    scale_y_continuous(
      name   = "Vaccine efficacy (death reduction)", 
      labels = percent,
      limits = c(0, 1), 
      expand = expansion(mult = c(0, 0.05)), 
      breaks = pretty_breaks()) + 
    # Prettify legend (needed for y spacing to take effect)...
    guides(fill = guide_legend(byrow = TRUE))
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.text     = element_text(size = 9),
          axis.title.x  = element_text(
            size = 16, margin = margin(b = 10, t = 20)),
          axis.title.y  = element_text(
            size = 16, margin = margin(l = 10, r = 20)),
          axis.line     = element_blank(),
          strip.text.x  = element_text(size = 13),
          strip.text.y  = element_text(size = 11),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(0.5, "lines"),
          panel.grid.major.y = element_line(linewidth = 0.25),
          legend.position = "none")
  
  # Save figure to file
  save_fig(g, "Vaccine efficacy profiles", dir = "static_models")
}

# ---------------------------------------------------------
# Plot effective coverage with waning immunity for static model pathogens
# ---------------------------------------------------------
plot_effective_coverage = function() {
  
  message("  - Plotting effective coverage by year and age")
  
  # Plot only up to a certain age
  age_max = 50
  
  # Manually define appropriate number of colours 
  colours = list(
    disease = "blues", 
    type    = "greens")
  
  # Number of breaks in continuous colour bar
  col_breaks = 40
  
  # Repeat for disease and vaccine type
  for (by in c("disease", "type")) {
    
    # Load previously calculated total coverage file
    effective_dt = read_rds("static", "effective_coverage", by)
    
    # Population weight over all countries
    plot_dt = effective_dt %>%
      append_d_v_t_name() %>%
      select(country, by = !!by, year, age, coverage) %>%
      filter(age >= 0,  # For simplicity, do not plot neonate immunity  
             age <= age_max) %>%
      # Append population size...
      left_join(y  = table("wpp_pop"),
                by = c("country", "year", "age")) %>%
      mutate(n = pop * coverage) %>%
      # Population weighted coverage...
      lazy_dt() %>%
      group_by(by, year, age) %>%
      summarise(effective_coverage = sum(n / sum(pop))) %>%
      ungroup() %>%
      as.data.table()
    
    # Construct continuous colour vector
    col_map = paste0("pals::brewer.", colours[by])
    col_vec = colour_scheme(col_map, n = col_breaks - 1)
    
    # Plot each pathogen by year ana age
    g = ggplot(plot_dt) + 
      aes(x = year, y = age, fill = effective_coverage) + 
      geom_tile() +
      facet_wrap(~by) + 
      # Set continuous colour bar...
      scale_fill_gradientn(
        colours = c("#FFFFFF", col_vec), 
        limits  = c(0, 1), 
        breaks  = pretty_breaks(), 
        label   = percent,
        guide   = guide_colourbar(
          title = "Effective coverage")) + 
      # Prettify x axis...
      scale_x_continuous(
        expand = c(0, 0), 
        breaks = seq(
          from = min(o$years), 
          to   = max(o$years), 
          by   = 5)) +
      # Prettify y axis...
      scale_y_continuous(
        name   = "Vaccine efficacy (death reduction)", 
        expand = c(0, 0), 
        breaks = pretty_breaks())
    
    # Prettify theme
    g = g + theme_classic() + 
      theme(axis.title.x  = element_blank(),
            axis.title.y  = element_text(
              size = 20, margin = margin(l = 10, r = 20)),
            axis.text     = element_text(size = 10),
            axis.text.x   = element_text(hjust = 1, angle = 50), 
            axis.line     = element_blank(),
            strip.text    = element_text(size = 14),
            strip.background = element_blank(), 
            panel.border  = element_rect(
              linewidth = 0.5, fill = NA),
            panel.spacing = unit(1, "lines"),
            legend.position    = "top", 
            legend.title  = element_text(
              size = 14, margin = margin(r = 20)),
            legend.text   = element_text(size = 10),
            legend.key.height = unit(1.5, "lines"),
            legend.key.width  = unit(6,   "lines"))
    
    # Save figure to file
    save_name = "Effective coverage"
    save_fig(g, save_name, by, dir = "static_models")
  }
}

# ---------------------------------------------------------
# Plot deaths and DALYs averted for static model pathogens
# ---------------------------------------------------------
plot_static = function() {
  
  message("  - Plotting static model impact results")
  
  # Flag for plotting recent extrapolation
  plot_extrap = FALSE
  
  # Disease burden / burden averted dictionary
  metric_dict = c(
    burden  = "estimated disease-specific burden (GBD 2019)", 
    averted = "estimated burden averted from static model")
  
  # Associated colours
  metric_colours = c("darkred", "navyblue")
  
  # ---- Plot by disease ----
  
  # Plot up to this year
  plot_to = ifelse(plot_extrap, max(o$years), max(o$gbd_prediction_years))
  
  # Repeat for each metric
  for (metric in o$metrics) {
    
    # Load previously calculated total coverage file
    averted_dt = read_rds("static", metric, "averted_disease")
    
    # Summarise results over country and age
    disease_dt = averted_dt %>%
      # lazy_dt() %>%
      filter(year <= plot_to) %>%
      pivot_longer(cols = c(burden, averted), 
                   names_to = "metric") %>%
      # Summarise over countries...
      group_by(disease, year, metric) %>%
      summarise(value = sum(value) / 1e6) %>%
      ungroup() %>%
      # Recode deaths disease/averted...
      left_join(y  = table("disease_name"), 
                by = "disease") %>%
      mutate(metric = recode(metric, !!!metric_dict), 
             metric = factor(metric, metric_dict)) %>%
      # Tidy up...
      select(metric, disease = disease_name, year, value) %>%
      arrange(metric, disease, year) %>%
      as.data.table()
    
    # Plot deaths and deaths averted by disease
    g = ggplot(disease_dt) + 
      aes(x = year, 
          y = value, 
          colour = metric) + 
      geom_line(linewidth = 2) + 
      facet_wrap(~disease) + 
      # Set colours and prettify legend...
      scale_colour_manual(values = metric_colours) + 
      guides(color = guide_legend(
        byrow = TRUE, ncol = 1)) +
      # Prettify x axis...
      scale_x_continuous(
        # limits = c(min(o$years), max(o$years)), 
        expand = expansion(mult = c(0, 0)), 
        breaks = pretty_breaks()) +  
      # Prettify y axis...
      scale_y_continuous(
        name   = "Number of people (in millions)", 
        labels = comma,
        expand = expansion(mult = c(0, 0.05)))
    
    # Set a figure title explaining metric
    g = g + ggtitle(
      label = table("metric_dict") %>%
        filter(metric == !!metric) %>%
        pull(metric_impact))
    
    # Prettify theme
    g = g + theme_classic() + 
      theme(plot.title    = element_text(size = 20, hjust = 0.5),
            axis.title.x  = element_blank(),
            axis.title.y  = element_text(
              size = 20, margin = margin(l = 10, r = 20)),
            axis.text     = element_text(size = 10),
            axis.text.x   = element_text(hjust = 1, angle = 50), 
            axis.line     = element_blank(),
            strip.text    = element_text(size = 14),
            strip.background = element_blank(), 
            panel.border  = element_rect(
              linewidth = 0.5, fill = NA),
            panel.spacing = unit(1, "lines"),
            panel.grid.major.y = element_line(linewidth = 0.5),
            legend.title  = element_blank(),
            legend.text   = element_text(size = 14),
            legend.key    = element_blank(),
            legend.position = "bottom", 
            legend.key.height = unit(2, "lines"),
            legend.key.width  = unit(3, "lines"))
    
    # Save figure to file
    save_fig(g, "Burden averted by disease", metric,
             dir = "static_models")
  }
  
  # ---- Plot by vaccine ----
  
  # Repeat for each metric
  for (metric in o$metrics) {
    
    # Load previously calculated total coverage file
    averted_dt = read_rds("static", metric, "averted_vaccine")
    
    # Full vaccine names
    name_dt = table("d_v_a") %>%
      filter(source == "static") %>%
      left_join(y  = table("vaccine_name"), 
                by = "vaccine") %>%
      select(d_v_a_id, name = vaccine_name)
    
    # Summarise results over country
    vaccine_dt = averted_dt %>%
      lazy_dt() %>%
      filter(year <= plot_to) %>%
      # Summarise over countries...
      group_by(d_v_a_id, year) %>%
      summarise(averted = sum(impact)) %>%
      ungroup() %>%
      # Convert from d_v_a to v_a...
      left_join(y  = name_dt, 
                by = "d_v_a_id") %>%
      select(vaccine = name, year, averted) %>%
      # Suitable vaccine order...
      mutate(vaccine = factor(vaccine, name_dt$name)) %>%
      arrange(vaccine, year) %>%
      as.data.table()
    
    # Plot deaths and deaths averted by vaccine
    g = ggplot(vaccine_dt) + 
      aes(x = year, 
          y = averted, 
          colour = vaccine) + 
      geom_line(linewidth   = 2, 
                show.legend = FALSE) + 
      # Facet by vaccine...
      facet_wrap(
        facets = vars(vaccine),
        scales = "free_y") + 
      # Prettify x axis...
      scale_x_continuous(
        # limits = c(min(o$years), max(o$years)), 
        expand = expansion(mult = c(0, 0)), 
        breaks = pretty_breaks()) +  
      # Prettify y axis...
      scale_y_continuous(
        name   = "Burden averted by vaccine", 
        labels = comma,
        limits = c(0, NA), 
        expand = expansion(mult = c(0, 0.05)))
    
    # Prettify theme
    g = g + theme_classic() + 
      theme(axis.title.x  = element_blank(),
            axis.title.y  = element_text(
              size = 18, margin = margin(l = 10, r = 20)),
            axis.text     = element_text(size = 9),
            axis.text.x   = element_text(hjust = 1, angle = 50), 
            axis.line     = element_blank(),
            strip.text    = element_text(size = 12),
            strip.background = element_blank(), 
            panel.border  = element_rect(
              linewidth = 0.5, fill = NA),
            panel.spacing = unit(1, "lines"),
            panel.grid.major.y = element_line(linewidth = 0.5))
    
    # Save figure to file
    save_fig(g, "Burden averted by vaccine", metric,
             dir = "static_models")
  }
}

# **** Regression (impute and infer) ****

# ---------------------------------------------------------
# Plot (any) correlation between covariates and imputation target
# ---------------------------------------------------------
plot_covariates = function(metric) {
  
  message("  - Plotting covariate-target relationships")
 # browser()
  # ---- Load data used for fitting ----
  
  # Function to load imputation data
  load_data_fn = function(id) {
    
    # Load imputation result for this d_v_a
    impute = read_rds("impute", "impute", metric, id)
    
    # If result is trivial, return trivial data
    if (is.null(impute$result)) {
      data = NULL 
      
    } else {  # Otherwise extract data
      data = impute$data %>%
        mutate(d_v_a_id = id, .before = 1)
    }
    
    return(data)
  }
  
  ignore = c(12)
  
  # Load imputation data for all d-v-a
  data_dt = table("d_v_a") %>%
    filter(source == "vimc") %>%
    filter(!d_v_a_id %in% ignore) %>% 
    pull(d_v_a_id) %>%
    lapply(load_data_fn) %>%
    rbindlist(fill= TRUE)
  
  # ---- Produce plot ----
  
  # Construct tidy plotting datatable
  plot_dt = data_dt %>%
    pivot_longer(cols = -c(d_v_a_id, target, country, year), 
                 names_to = "covariate") %>%
    arrange(d_v_a_id, covariate, target) %>%
    format_d_v_a_name() %>%
    select(d_v_a_name, target, covariate, value) %>%
    as.data.table()
 # browser()
  
  # Plot covariates vs imputation target
  g = ggplot(plot_dt) +
    aes(x = (target*1e5), 
        y = value, 
        colour = covariate) +
    geom_point(
      alpha = 0.2, 
      shape = 16, 
      show.legend = FALSE) +
    # Facet in 2 dims with label wrapping...
    facet_grid(
      rows     = vars(covariate), 
      cols     = vars(d_v_a_name),
      labeller = label_wrap_gen(width = 20), 
      scales   = "free") + 
    # Set colour scheme...
    scale_colour_manual(
      values = colour_scheme(
        map = "pals::kovesi.rainbow", 
        n   = n_unique(plot_dt$covariate))) +
    # Prettify x axis...
    scale_x_continuous(
      name   = "Normalised predictor", 
      limits = c(0, 1), 
      expand = c(0, 0), 
      breaks = pretty_breaks()) +  
    # Prettify y axis...
    scale_y_continuous(
      name   = "Normalised response (impact per FVP)", 
      limits = c(0, 1), 
      expand = c(0, 0), 
      breaks = pretty_breaks())
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.text     = element_text(size = 6),
          axis.text.x   = element_text(hjust = 1, angle = 50),
          axis.title.x  = element_text(
            size = 16, margin = margin(l = 4, r = 8)),
          axis.title.y  = element_text(
            size = 16, margin = margin(b = 4, t = 8)),
          axis.line     = element_blank(),
          strip.text    = element_text(size = 8),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(0.4, "lines"))
  
  # Save figure to file
  save_fig(g, "Covariate relationships", dir = "imputation")
}

# ---------------------------------------------------------
# Plot truth vs predicted for imputation training data
# ---------------------------------------------------------
plot_impute_quality = function(metric) {
 # browser()
  message("  - Plotting imputation quality of fit")
  
  # ---- Load results from fitting ----
  
  # Function to load imputation results
  load_results_fn = function(id)
    result = read_rds("impute", "impute", metric, id)$result
  
 # ignore = 12
  
  # Load imputation results for all d-v-a
  results_dt = table("d_v_a") %>%
    filter(source == "vimc") %>%
   # filter(!d_v_a_id %in% ignore) %>% 
    pull(d_v_a_id) %>%
    lapply(load_results_fn) %>%
    rbindlist()
  
  # ---- Construct plotting datatables ----
  
  # Prepare datatable for plotting
  plot_dt = results_dt %>%
    filter(!is.na(target)) %>%
    filter(!is.na(prediction)) %>%
    filter(target != 0) %>%
    select(-country) %>%
    as.data.table()
  
  # Maximum value in each facet (target or prediction)
  blank_dt = plot_dt %>%
    mutate(max_value = pmax(target, prediction)) %>%
    group_by(d_v_a_id) %>%
    summarise(max_value = max(max_value)) %>%
    ungroup() %>%
    expand_grid(type = c("target", "prediction")) %>%
    pivot_wider(names_from  = type, 
                values_from = max_value) %>%
    as.data.table()
  
  # ---- Produce plot ----
  
  # Single plot with multiple facets
  g = ggplot(plot_dt) +
    aes(x = target, 
        y = prediction, 
        color = d_v_a_id) +
    # Plot truth vs predicted...
    geom_point(alpha = 0.35, 
               shape = 16, 
               show.legend = FALSE) +
    # For square axes...
    geom_blank(data = blank_dt) +
    # To see quality of predict vs target...
    geom_abline(colour = "black") + 
    # Simple faceting with wrap labelling...
    facet_wrap(
      facets   = vars(d_v_a_id), 
      labeller = label_wrap_gen(width = 30), 
      scales   = "free") + 
    # Prettify x axis...
    scale_x_continuous(
      name   = "Imputation target", 
      labels = scientific,
      limits = c(0, NA), 
      expand = c(0, 0), 
      breaks = pretty_breaks()) +  
    # Prettify y axis...
    scale_y_continuous(
      name   = "Imputation prediction", 
      labels = scientific,
      limits = c(0, NA),
      expand = c(0, 0), 
      breaks = pretty_breaks())
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.text     = element_text(size = 8),
          axis.text.x   = element_text(hjust = 1, angle = 50),
          axis.title.x  = element_text(
            size = 18, margin = margin(l = 10, r = 20)),
          axis.title.y  = element_text(
            size = 18, margin = margin(b = 10, t = 20)),
          axis.line     = element_blank(),
          strip.text    = element_text(size = 12),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(0.5, "lines"))
  
  # Save figure to file
  save_fig(g, "Imputation fit", metric, dir = "imputation")
}

# --------------------------------------------------------
# Plot predictive performance for each country
# --------------------------------------------------------
plot_impute_perform = function(metric) {
 # browser()
  message("  - Plotting predictive performance by country")
  
  # ---- Load models from fitting ----
  # Function to load best model for each country and show results
  load_results_fn = function(id){
    model = read_rds("impute", "impute", metric, id)$model 
  }
  
 # ignore = 12
  
  # Load imputation results for all d-v-a
  diseases_dt = table("d_v_a") %>%
    filter(source == "vimc") %>%
   # filter(!d_v_a_id %in% ignore) %>% 
    pull(d_v_a_id)
  
  for(id in diseases_dt){
    model = load_results_fn(id)
    
    plot_dt = augment(model) %>%
      rename(prediction = .fitted)
    
    # Maximum value in each facet (target or prediction)
    blank_dt = plot_dt %>%
      mutate(max_value = pmax(target, prediction)) %>%
      group_by(country) %>%
      summarise(max_value = max(max_value)) %>%
      ungroup() %>%
      expand_grid(type = c("target", "prediction")) %>%
      pivot_wider(names_from  = type, 
                  values_from = max_value) %>%
      as.data.table()
    
    
    # ---- Produce plot ----
    
    # Single plot with multiple facets
    g = ggplot(plot_dt) +
      aes(x = target, 
          y = prediction) +
      # Plot truth vs predicted...
      geom_point(alpha = 0.35, 
                 shape = 16) +
      # For square axes...
      geom_blank(data = blank_dt) +
      # x=y
      geom_abline(colour = "black", 
                  intercept = 0,
                  slope = 1) + 
      # Simple faceting with wrap labelling...
      facet_wrap(
        facets   = ~country, 
        labeller = label_wrap_gen(width = 30), 
        ncol = 21,
        scales = "free") + 
      # Prettify x axis...
      scale_x_continuous(
        name   = "Imputation target", 
        labels = NULL,#scientific,
        limits = c(0, NA), 
        expand = c(0, 0), 
        breaks = pretty_breaks()) +  
      # Prettify y axis...
      scale_y_continuous(
        name   = "Imputation prediction", 
        labels = NULL,#scientific,
        limits = c(0, NA),
        expand = c(0, 0), 
        breaks = pretty_breaks()) +
      # Title 
      labs(title = paste("Predictive performance for", plot_dt$d_v_a_id))
    
    # Prettify theme
    g = g + theme_classic() + 
      theme(axis.text     = element_text(size = 8),
            axis.text.x   = element_text(hjust = 1, angle = 50),
            axis.title.x  = element_text(
              size = 16, margin = margin(l = 10, r = 20)),
            axis.title.y  = element_text(
              size = 16, margin = margin(b = 10, t = 20)),
            axis.line     = element_blank(),
            strip.text    = element_text(size = 12),
            strip.background = element_blank(), 
            panel.border  = element_rect(
              linewidth = 0.5, fill = NA),
            panel.spacing = unit(0.5, "lines"))
    
    
    # Details for file destination
    save_name = "Predictive performance by country"
    save_dir  = "imputation"
    
    # Save figure to file
    save_fig(g, save_name, id, metric, dir = save_dir)
  }
  
  return()
}

#-------------------------------------------------
# Plot model choice per region
#-------------------------------------------------
plot_model_choice = function(metric) {
  
  message("  - Plotting model choice by region")
  
  # ---- Load results from imputation ----
   # Function to load imputation results
  load_results_fn = function(id)
    result = read_rds("impute", "impute", metric, id)$model
  
  #ignore = 12
  
  # Load imputation results for all d-v-a
  choice_dt = table("d_v_a") %>%
    filter(source == "vimc") %>%
   # filter(!d_v_a_id %in% ignore) %>% 
    pull(d_v_a_id) %>%
    lapply(load_results_fn) %>%
    rbindlist() 
 
  # Group model choice by income 
   income_choice_dt = choice_dt %>%
                       group_by(d_v_a_id) %>%
                       count(income, model_id) %>%
                       group_by(model_id) %>%         
                       mutate(prop = prop.table(n)) 
   
  # Group model choice by region 
   region_choice_dt = choice_dt %>%
                       group_by(d_v_a_id) %>%
                       count(region, model_id) %>%
                       group_by(model_id) %>%         
                       mutate(prop = prop.table(n)) 
  
  # Details for file destination
  save_name = paste("Model_choice_by_region", metric, sep = "_")
  save_dir  = "data_visualisation"
  
  # Plot model choice by region
  g1 = ggplot(region_choice_dt,
              aes(x = model_id,
                  y = n,
                  fill = region)) + 
              geom_bar(stat = "identity") + 
  
  # Simple faceting with wrap labelling...
  facet_wrap(
    facets   = vars(d_v_a_id), 
    labeller = label_wrap_gen(width = 30), 
    scales   = "free") 
  
  # Prettify theme
  g1 = g1 + theme_classic() + 
    theme(axis.text     = element_text(size = 10),
          axis.title.x  = element_text(
            size = 20, margin = margin(b = 10, t = 20)),
          axis.title.y  = element_text(
            size = 20, margin = margin(l = 10, r = 20)),
          axis.line     = element_blank(),
          strip.text    = element_text(size = 14),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 1, fill = NA),
          panel.spacing = unit(1, "lines"),
          legend.title  = element_blank(),
          legend.text   = element_text(size = 14),
          legend.key    = element_blank(),
          legend.position = "right", 
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines"))
  
  
  # Plot model choice by income level
  g2 = ggplot(income_choice_dt,
              aes(x = model_id,
                  y = n,
                  fill = income)) + 
    geom_bar(stat = "identity") + 
    
    # Simple faceting with wrap labelling...
    facet_wrap(
      facets   = vars(d_v_a_id), 
      labeller = label_wrap_gen(width = 30), 
      scales   = "free") 
  
  # Prettify theme
  g2 = g2 + theme_classic() + 
    theme(axis.text     = element_text(size = 10),
          axis.title.x  = element_text(
            size = 20, margin = margin(b = 10, t = 20)),
          axis.title.y  = element_text(
            size = 20, margin = margin(l = 10, r = 20)),
          axis.line     = element_blank(),
          strip.text    = element_text(size = 14),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 1, fill = NA),
          panel.spacing = unit(1, "lines"),
          legend.title  = element_blank(),
          legend.text   = element_text(size = 14),
          legend.key    = element_blank(),
          legend.position = "right", 
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines"))
  
  # Save to file
  save_fig(g1, paste(save_name, "region", sep = "_"), dir = save_dir)
  save_fig(g2, paste(save_name, "income", sep = "_"), dir = save_dir)
  
  
}

# --------------------------------------------------------
# Plot fitted model for each country
# --------------------------------------------------------
plot_impute_fit = function(metric){
  message("  - Plotting model fit by country")
 
  # ---- Load models from fitting ----
  # Function to load best model for each country and show results
  load_results_fn = function(id){
    model = read_rds("impute", "impute", metric, id)$model 
  }
  
  # Load imputation results for all d-v-a
  diseases_dt = table("d_v_a") %>%
     filter(source == "vimc") %>%
     pull(d_v_a_id)
  
  for(id in diseases_dt){
    model = load_results_fn(id)
    
    plot_dt = augment(model) %>%
      rename(prediction = .fitted)
    
    # Maximum value in each facet (target or prediction)
    blank_dt = plot_dt %>%
      mutate(max_value = pmax(target, prediction)) %>%
      group_by(country) %>%
      summarise(max_value = max(max_value)) %>%
      ungroup() %>%
      expand_grid(type = c("target", "prediction")) %>%
      pivot_wider(names_from  = type, 
                  values_from = max_value) %>%
      as.data.table()
    
    
    # ---- Produce plot ----
    
    # Single plot with multiple facets
    g = ggplot(plot_dt) +
      aes(x = year) +
      # Plot fitting data
      geom_point(aes(y = target,
                     fill = "#6CA2EA")) +
      
      # Plot model output
      geom_line(aes(y = prediction,
                    colour = "black")) +
      
      # For square axes...
      geom_blank(data = blank_dt) +
      
      # Simple faceting with wrap labelling...
      facet_wrap(
        facets   = ~country, 
        labeller = label_wrap_gen(width = 30), 
        ncol = 21,
        scales = "free") + 
      # Prettify x axis...
      scale_x_continuous(
        name   = "Year", 
        labels = waiver(),
        #limits = c(1990, 2024), 
        expand = c(0, 0), 
        breaks = pretty_breaks()) +  
      # Prettify y axis...
      scale_y_continuous(
        name   = "Impact", 
        labels = NULL,#scientific,
        limits = c(0, NA),
        expand = c(0, 0), 
        breaks = pretty_breaks()) +
      # Title 
      labs(title = paste("Fitted model for", plot_dt$d_v_a_id))
    
    # Prettify theme
    g = g + theme_classic() + 
      theme(axis.text     = element_text(size = 8),
            axis.text.x   = element_text(hjust = 1, angle = 50),
            axis.title.x  = element_text(
              size = 16, margin = margin(l = 10, r = 20)),
            axis.title.y  = element_text(
              size = 16, margin = margin(b = 10, t = 20)),
            axis.line     = element_blank(),
            strip.text    = element_text(size = 12),
            strip.background = element_blank(), 
            panel.border  = element_rect(
              linewidth = 0.5, fill = NA),
            panel.spacing = unit(0.5, "lines"))
    
    
    # Details for file destination
    save_name = "Model fit by country"
    save_dir  = "imputation"
    
    # Save figure to file
    save_fig(g, save_name, id, metric, dir = save_dir)
  }
  
  return()
}

# --------------------------------------------------------
# Plot fitted model for validation countries
# --------------------------------------------------------
plot_validation = function(metric){
  message("  - Plotting validation")
  browser()
  # ---- Load models from fitting ----
  # Function to load best model for each country and show results
  load_results_fn = function(id){
    model = read_rds("impute", "impute", metric, id)$model 
  }
  
   # Load imputation results for all d-v-a
  diseases_dt = table("d_v_a") %>%
    filter(source == "vimc") %>%
    pull(d_v_a_id)
  
  for(id in diseases_dt){
    model = load_results_fn(id)
    
    plot_dt = augment(model) %>%
      rename(prediction = .fitted)
    
    # Maximum value in each facet (target or prediction)
    blank_dt = plot_dt %>%
      mutate(max_value = pmax(target, prediction)) %>%
      group_by(country) %>%
      summarise(max_value = max(max_value)) %>%
      ungroup() %>%
      expand_grid(type = c("target", "prediction")) %>%
      pivot_wider(names_from  = type, 
                  values_from = max_value) %>%
      as.data.table()
    
    
    # ---- Produce plot ----
    
    # Single plot with multiple facets
    g = ggplot(plot_dt) +
      aes(x = year) +
      # Plot fitting data
      geom_point(aes(y = target,
                     fill = "#6CA2EA")) +
      
      # Plot model output
      geom_line(aes(y = prediction,
                    colour = "black")) +
      
      # For square axes...
      geom_blank(data = blank_dt) +
      
      # Simple faceting with wrap labelling...
      facet_wrap(
        facets   = ~country, 
        labeller = label_wrap_gen(width = 30), 
        ncol = 21,
        scales = "free") + 
      # Prettify x axis...
      scale_x_continuous(
        name   = "Year", 
        labels = waiver(),
       # limits = c(1990, 2024), 
        expand = c(0, 0), 
        breaks = pretty_breaks()) +  
      # Prettify y axis...
      scale_y_continuous(
        name   = "Impact", 
        labels = NULL,#scientific,
        limits = c(0, NA),
        expand = c(0, 0), 
        breaks = pretty_breaks()) +
      # Title 
      labs(title = paste("Fitted model for", plot_dt$d_v_a_id))
    
    # Prettify theme
    g = g + theme_classic() + 
      theme(axis.text     = element_text(size = 8),
            axis.text.x   = element_text(hjust = 1, angle = 50),
            axis.title.x  = element_text(
              size = 16, margin = margin(l = 10, r = 20)),
            axis.title.y  = element_text(
              size = 16, margin = margin(b = 10, t = 20)),
            axis.line     = element_blank(),
            strip.text    = element_text(size = 12),
            strip.background = element_blank(), 
            panel.border  = element_rect(
              linewidth = 0.5, fill = NA),
            panel.spacing = unit(0.5, "lines"))
    
    
    # Details for file destination
    save_name = "Model fit for validation countries"
    save_dir  = "imputation"
    
    # Save figure to file
    save_fig(g, save_name, id, dir = save_dir)
  }
  
  return()
}
# **** Impact functions ****

# ---------------------------------------------------------
# Exploratory plots of data used to fit impact functions
# ---------------------------------------------------------
plot_impact_data = function(metric) {
  
  message("  - Plotting impact function fitting data")
  
  # Load data used for impact function fitting
  data_dt = read_rds("impact", "impact", metric, "data") %>%
    format_d_v_a_name()
  
  # ---- Plot 1: impact per FVP over time ----
  
  # Impact per FVP over time
  g1 = ggplot(data_dt) +
    aes(x = year, 
        y = impact_fvp, 
        colour = country) +
    geom_line(show.legend = FALSE) +
    # Faceting with wrap labelling...
    facet_wrap(
      facets   = vars(d_v_a_name), 
      labeller = label_wrap_gen(width = 30), 
      scales   = "free_y") + 
    # Set colour scheme...
    scale_colour_manual(
      values = colour_scheme(
        map = "viridis::viridis", 
        n   = n_unique(all_countries()))) +
    # Prettify x axis...
    scale_x_continuous(
      expand = c(0, 0), 
      breaks = pretty_breaks()) +  
    # Prettify y axis...
    scale_y_continuous(
      name   = "Impact per fully vaccinated", 
      labels = scientific,
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05)), 
      breaks = pretty_breaks())
  
  # Prettify theme
  g1 = g1 + theme_classic() + 
    theme(axis.text     = element_text(size = 8),
          axis.text.x   = element_text(hjust = 1, angle = 50),
          axis.title.x  = element_text(
            size = 16, margin = margin(l = 10, r = 20)),
          axis.title.y  = element_text(
            size = 16, margin = margin(b = 10, t = 20)),
          axis.line     = element_blank(),
          strip.text    = element_text(size = 10),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(0.5, "lines"))
  
  # Save figure to file
  save_fig(g1, "Data", "impact ratio", metric, 
           dir = "impact_functions")
  
  # ---- Plot 2: cumulative impact vs cumulative FVP ----
  
  # Cumulative FVPs vs cumulative deaths averted
  g2 = ggplot(data_dt) + 
    aes(x = fvps, 
        y = impact, 
        colour = country) +
    geom_line(show.legend = FALSE) +
    # Faceting with wrap labelling...
    facet_wrap(
      facets   = vars(d_v_a_name), 
      labeller = label_wrap_gen(width = 30), 
      scales   = "free") + 
    # Set colour scheme...
    scale_colour_manual(
      values = colour_scheme(
        map = "viridis::viridis", 
        n   = n_unique(all_countries()))) +
    # Prettify x axis...
    scale_x_continuous(
      name   = "Cumulative FVPs per capita (including new birth cohorts)", 
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05)), 
      breaks = pretty_breaks()) +  
    # Prettify y axis...
    scale_y_continuous(
      name   = "Cumulative impact per capita (including new birth cohorts)", 
      labels = scientific,
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05)), 
      breaks = pretty_breaks())
  
  # Prettify theme
  g2 = g2 + theme_classic() + 
    theme(axis.text     = element_text(size = 8),
          axis.title.x  = element_text(
            size = 16, margin = margin(l = 10, r = 20)),
          axis.title.y  = element_text(
            size = 16, margin = margin(b = 10, t = 20)),
          axis.line     = element_blank(),
          strip.text    = element_text(size = 10),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(0.5, "lines"))
  
  # Save figure to file
  save_fig(g2, "Data", "cumulative FVP vs impact", metric, 
           dir = "impact_functions")
}

# ---------------------------------------------------------
# Plot impact ratios - either all or initial only
# ---------------------------------------------------------
plot_impact_fvps = function(metric, scope) {
  
  width = 0.3
  
  # Function for averaging (mean or median)
  avg_fn = get("mean")
  
  # ---- Load data based on scope ----
  
  # All time plot
  if (scope == "all_time") {
    
    message("  - Plotting all-time impact per FVP: ", metric)
    
    # Load initial ratio data
    impact_dt = read_rds("history", "impact", metric, "data")
    
    # Set a descriptive y-axis title
    y_lab = "Impact per fully vaccinated person (log10 scale)"
  }
  
  # Initial year plot
  if (scope == "initial") {
    
    message("  - Plotting initial impact per FVP: ", metric)
    
    # Load initial ratio data
    impact_dt = read_rds("history", "initial_ratio", metric) %>%
      rename(impact_fvp = initial_ratio)
    
    # Set a descriptive y-axis title
    y_lab = "Initial impact per FVP used for back projection (log10 scale)"
  }
  
  # ---- Classify by income status ----
  
  # Load income status dictionary
  income_dict = table("income_dict")
  
  # Load income status of each country
  income_dt = table("income_status") %>%
    filter(year == max(year)) %>%
    left_join(y  = income_dict, 
              by = "income") %>%
    select(country, income = income_name)
  
  # Classify by income group
  data_dt = impact_dt %>%
    filter(impact_fvp > 0) %>%
    # Append d_v_a description...
    format_d_v_a_name() %>%
    # Append income status description...
    left_join(y  = income_dt, 
              by = "country") %>%
    mutate(income = factor(
      x      = income, 
      levels = rev(income_dict$income_name))) %>%
    select(d_v_a = d_v_a_name, income, impact_fvp) %>%
    unique()
  
  # ---- Plotting coordinates ----
  
  # Set x values for each d_v_a
  x_major = data_dt %>%
    group_by(d_v_a) %>%
    summarise(order = avg_fn(impact_fvp)) %>%
    ungroup() %>%
    arrange(desc(order)) %>%
    mutate(x_major = 1 : n()) %>%
    select(-order) %>%
    as.data.table()
  
  # Offset each income group
  x_minor = data_dt %>%
    select(income) %>%
    unique() %>%
    arrange(income) %>%
    mutate(x_minor = seq(
      from = -width / 2, 
      to   = width / 2, 
      length.out = n()))
  
  # Width of offset of points within income group
  x_spray = width / (nrow(income_dict) * 2)
  
  # Bring all x and y values together
  plot_dt = data_dt %>%
    left_join(y  = x_major, 
              by = "d_v_a") %>%
    left_join(y  = x_minor, 
              by = "income") %>%
    mutate(x = x_major + x_minor) %>%
    select(d_v_a, income, x, y = impact_fvp) %>%
    arrange(x)
  
  # Extract bounds (after transformation)
  lb = floor(min(log10(data_dt$impact_fvp)))
  ub = ceiling(max(log10(data_dt$impact_fvp)))
  
  # ---- Vaccine and income status averages ----
  
  # Average by vaccine and income group
  income_average_dt = plot_dt %>%
    group_by(income, x) %>%
    summarise(y = avg_fn(y)) %>%
    ungroup() %>%
    arrange(x) %>%
    as.data.table()
  
  # Average by vaccine (over all income groups)
  vaccine_average_dt = data_dt %>%
    group_by(d_v_a) %>%
    summarise(y = avg_fn(impact_fvp)) %>%
    ungroup() %>%
    left_join(y  = x_major,
              by = "d_v_a") %>%
    select(d_v_a, x = x_major, y) %>%
    arrange(x) %>%
    as.data.table()
  
  # ---- Produce primary plot ----
  
  # Plot impact per FVP
  g = ggplot(plot_dt) +
    aes(x = x, y = y) +
    # Plot all points by vaccine and income...
    geom_point(
      mapping  = aes(colour = income),
      alpha    = 0.1,
      shape    = 16,
      stroke   = 0,
      position = position_jitter(
        width = x_spray,
        seed  = 1)) + 
    # Average of each vaccine...
    geom_point(
      data   = vaccine_average_dt,
      colour = "black",
      shape = 95, # Horizontal lines
      size  = 20) +
    # Average of each income group...
    geom_point(
      data    = income_average_dt, 
      mapping = aes(fill = income), 
      color   = "black", 
      shape   = 23, 
      size    = 3) + 
    # Prettify x axis...
    scale_x_continuous(
      breaks = x_major$x_major, 
      labels = x_major$d_v_a) +
    # Prettify y axis...
    scale_y_continuous(
      name   = y_lab,
      trans  = "log10",
      labels = scientific, 
      limits = c(10 ^ lb, 10 ^ ub),
      expand = c(0, 0),
      breaks = 10 ^ rev(ub : lb))
  
  # ---- Prettify and save ----
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.text     = element_text(size = 10),
          axis.text.x   = element_text(hjust = 1, angle = 50),
          axis.title.x  = element_blank(),
          axis.title.y  = element_text(
            size = 16, margin = margin(l = 10, r = 20)),
          axis.line     = element_blank(),
          panel.grid.major.y = element_line(linewidth = 0.25),
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(1, "lines"))
  
  # Filename to save to depending on scope
  save_stem  = "Density of impact per FVP"
  save_scope = str_replace(scope, "_", " ")
  
  # Save these figures to file
  save_fig(g, save_stem, save_scope, metric, 
           dir = "impact_functions")
  
  # Save all time plot as key manuscript figure
  # if (scope == "all_time" && metric == "deaths")
  #   save_fig(g, "Figure 3", dir = "manuscript")
}

# ---------------------------------------------------------
# Plot function selection statistics
# ---------------------------------------------------------
plot_model_selection = function(metric) {
  
  message("  - Plotting impact model selection")
  
  # Define colour scheme
  colour_map = "brewer::set1"
  
  # Load stuff: best fit functions and associtaed coefficients
  choice_dt = read_rds("impact", "model_choice", metric) %>%
    format_d_v_a_name()
  
  # ---- Plot function count ----
  
  # Simple plotting function with a few features
  plot_selection = function(var, type = "count", stat = "number") {
    
    # Determine order - with 'focus' function first
    fn_dict = fn_set(dict = TRUE)
    
    # Number of times each model is optimal
    selection_dt = choice_dt %>% 
      rename(var = !!var) %>% 
      # Number and proportion of each fn...
      count(var, fn, name = "number") %>%
      group_by(var) %>%
      mutate(total = sum(number)) %>%
      ungroup() %>%
      mutate(proportion = number / total) %>%
      # Set appropriate plotting order...
      rename(val = !!stat) %>%
      select(var, fn, val) %>%
      pivot_wider(names_from  = fn, 
                  values_from = val, 
                  values_fill = 0) %>%
      arrange_at(intersect(
        x = names(fn_dict), 
        y = names(.))) %>%
      # Final formatting...
      pivot_longer(cols = -var, 
                   names_to  = "fn", 
                   values_to = "val") %>%
      mutate(fn  = recode(fn, !!!fn_dict), 
             fn  = factor(fn, fn_dict),
             var = fct_rev(fct_inorder(var))) %>%
      as.data.table()
    
    # Check figure type flag
    if (type == "count") {
      
      # Number of occurances
      g = ggplot(selection_dt) + 
        aes(x    = var, 
            y    = val, 
            fill = fct_rev(fn)) + 
        geom_col() + 
        coord_flip() + 
        # Set colour scheme...
        scale_fill_manual(
          values = colour_scheme(
            map = colour_map, 
            n   = n_unique(selection_dt$fn))) +
        # Prettify x axis (noting coord_flip)...
        scale_y_continuous(
          name   = paste(first_cap(stat), "of countries"), 
          expand = expansion(mult = c(0, 0.05)), 
          breaks = pretty_breaks()) + 
        # Prettify legend...
        guides(fill = guide_legend(
          reverse = TRUE))
      
      # Prettify theme
      g = g + theme_classic() + 
        theme(axis.text     = element_text(size = 12),
              axis.title.x  = element_text(
                size = 20, margin = margin(b = 10, t = 20)),
              axis.title.y  = element_blank(),
              axis.line     = element_blank(),
              panel.border  = element_rect(
                linewidth = 0.5, fill = NA),
              legend.title  = element_blank(),
              legend.text   = element_text(size = 12),
              legend.key    = element_blank(),
              legend.position = "right", 
              legend.key.height = unit(2, "lines"),
              legend.key.width  = unit(2, "lines"))
    }
    
    # Check figure type flag
    if (type == "density") {
      
      # Density of occurances
      g = ggplot(selection_dt) + 
        aes(x    = val, 
            fill = fn) +
        geom_bar() + 
        # Prettify x axis...
        scale_x_continuous(
          name   = "Number of model selections for a country", 
          breaks = pretty_breaks()) +  
        # Prettify y axis...
        scale_y_continuous(
          name   = paste(first_cap(stat), "of countries"), 
          expand = expansion(mult = c(0, 0.05)), 
          breaks = pretty_breaks())
      
      # Prettify theme
      g = g + theme_classic() + 
        theme(axis.text     = element_text(size = 12),
              axis.title.x  = element_text(
                size = 20, margin = margin(b = 10, t = 20)),
              axis.title.y  = element_text(
                size = 20, margin = margin(l = 10, r = 20)),
              axis.line     = element_blank(),
              panel.border  = element_rect(
                linewidth = 0.5, fill = NA),
              legend.title  = element_blank(),
              legend.text   = element_text(size = 12),
              legend.key    = element_blank(),
              legend.position = "right", 
              legend.key.height = unit(2, "lines"),
              legend.key.width  = unit(2, "lines"))
    }
    
    return(g)
  }
  
  # ---- A variety of plots ----
  
  # Figure sub-directory to save to 
  save_dir = "impact_functions"
  
  # Plot by disease-vaccine-activity
  g1 = plot_selection("d_v_a_name", stat = "number")
  g2 = plot_selection("d_v_a_name", stat = "proportion")
  
  # Filename stem
  save_name = "Selection by pathogen"
  
  # Save both figures
  save_fig(g1, save_name, metric, "number",     dir = save_dir)
  save_fig(g2, save_name, metric, "proportion", dir = save_dir)
  
  # Plot by country
  g3 = plot_selection("country", type = "density")
  
  # Save the last figure
  save_name = "Selection density by country"
  save_fig(g3, save_name, metric, dir = save_dir)
}

# ---------------------------------------------------------
# Plot impact function evaluation
# ---------------------------------------------------------
plot_model_fits = function(metric) {
  
  message("  - Plotting impact function fits")
  
  # Load data used for impact function fitting
  data_dt = read_rds("impact", "impact", metric, "data") %>%
    format_d_v_a_name()
  
  # Evaluate only as far as we have data
  max_data = data_dt %>%
    group_by(d_v_a_name, country) %>%
    slice_max(fvps, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(d_v_a_name, country, x_max = fvps) %>%
    # Increment up one so we plot slightly past the data...
    mutate(x_max = x_max + o$eval_x_scale / 100) %>%
    as.data.table()
  
  browser() # Need to load up mean coef datatable here
  
  # Evaluate selected impact function
  best_fit = evaluate_impact_function(
    data = NULL, 
    coef = NULL) %>%
    format_d_v_a_name()
  
  # Apply max_data so we only plot up to data (or just past)
  plot_dt = best_fit %>%
    left_join(y  = max_data,
              by = c("d_v_a_name", "country")) %>%
    filter(fvps < x_max) %>%
    format_d_v_a_name() %>%
    select(d_v_a_name, country, fvps, impact)
  
  # Plot function evaluation against the data
  g = ggplot(plot_dt) +
    aes(x = fvps, 
        y = impact, 
        colour = country) +
    # Plot data, then fit on top...
    geom_point(
      data = data_dt,
      size = 0.75,
      alpha = 0.5,
      show.legend = FALSE) +
    geom_line(show.legend = FALSE) +
    # Faceting with wrap labelling...
    facet_wrap(
      facets   = vars(d_v_a_name), 
      labeller = label_wrap_gen(width = 30), 
      scales   = "free") + 
    # Set colour scheme...
    scale_colour_manual(
      values = colour_scheme(
        map = "viridis::viridis", 
        n   = n_unique(all_countries()))) +
    # Prettify x axis...
    scale_x_continuous(
      name   = "Cumulative FVPs per capita (including new birth cohorts)", 
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05)), 
      breaks = pretty_breaks()) +  
    # Prettify y axis...
    scale_y_continuous(
      name   = "Cumulative impact per capita", 
      labels = scientific,
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05)), 
      breaks = pretty_breaks())
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.text     = element_text(size = 8),
          axis.text.x   = element_text(hjust = 1, angle = 50),
          axis.title.x  = element_text(
            size = 16, margin = margin(b = 10, t = 20)),
          axis.title.y  = element_text(
            size = 16, margin = margin(l = 10, r = 20)),
          axis.line     = element_blank(),
          strip.text    = element_text(size = 10),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(0.5, "lines"))
  
  save_fig(g, "Impact function evaluation", metric, 
           dir = "impact_functions")
}

# ---------------------------------------------------------
# Plot impact vs coverage by vaccine, income, and decade 
# ---------------------------------------------------------
plot_impact_coverage = function(metric) {
  
  message("  - Plotting impact against coverage")
  
  browser()
  
  # Function for averaging (mean or median)
  avg_fn = get("mean")
  
  # ---- Classify by income status ----
  
  # Load initial ratio data
  impact_dt = read_rds("impact", "data")
  
  # Load income status dictionary
  income_dict = table("income_dict")
  
  # Load income status of each country
  income_dt = table("income_status") %>%
    filter(year == max(year)) %>%
    left_join(y  = income_dict, 
              by = "income") %>%
    select(country, income = income_name)
  
  browser()
  
  # Classify by income group
  data_dt = impact_dt %>%
    filter(impact_fvp > 0) %>%
    # Append income status description...
    left_join(y  = income_dt, 
              by = "country") %>%
    mutate(income = factor(
      x      = income, 
      levels = rev(income_dict$income_name))) %>%
    select(d_v_a = d_v_a_name, income, impact_fvp) %>%
    unique()
  
  # ---- Plotting coordinates ----
  
  # Set x values for each d_v_a
  x_major = data_dt %>%
    group_by(d_v_a) %>%
    summarise(order = avg_fn(impact_fvp)) %>%
    ungroup() %>%
    arrange(desc(order)) %>%
    mutate(x_major = 1 : n()) %>%
    select(-order) %>%
    as.data.table()
  
  # Offset each income group
  x_minor = data_dt %>%
    select(income) %>%
    unique() %>%
    arrange(income) %>%
    mutate(x_minor = seq(
      from = -width / 2, 
      to   = width / 2, 
      length.out = n()))
  
  # Width of offset of points within income group
  x_spray = width / (nrow(income_dict) * 2)
  
  # Bring all x and y values together
  plot_dt = data_dt %>%
    left_join(y  = x_major, 
              by = "d_v_a") %>%
    left_join(y  = x_minor, 
              by = "income") %>%
    mutate(x = x_major + x_minor) %>%
    select(d_v_a, income, x, y = impact_fvp) %>%
    arrange(x)
  
  # Extract bounds (after transformation)
  lb = floor(min(log10(data_dt$impact_fvp)))
  ub = ceiling(max(log10(data_dt$impact_fvp)))
  
  # ---- Vaccine and income status averages ----
  
  income_average_dt = plot_dt %>%
    group_by(income, x) %>%
    summarise(y = avg_fn(y)) %>%
    ungroup() %>%
    arrange(x) %>%
    as.data.table()
  
  vaccine_average_dt = data_dt %>%
    group_by(d_v_a) %>%
    summarise(y = avg_fn(impact_fvp)) %>%
    ungroup() %>%
    left_join(y  = x_major,
              by = "d_v_a") %>%
    select(d_v_a, x = x_major, y) %>%
    arrange(x) %>%
    as.data.table()
  
  # ---- Produce primary plot ----
  
  # Plot impact per FVP
  g = ggplot(plot_dt) +
    aes(x = x, y = y) +
    # Plot all points by vaccine and income...
    geom_point(
      mapping  = aes(colour = income),
      alpha    = 0.1,
      shape    = 16,
      stroke   = 0,
      position = position_jitter(
        width = x_spray,
        seed  = 1)) + 
    # Average of each vaccine...
    geom_point(
      data   = vaccine_average_dt,
      colour = "black",
      shape = 95, # Horizontal lines
      size  = 20) +
    # Average of each income group...
    geom_point(
      data    = income_average_dt, 
      mapping = aes(fill = income), 
      color   = "black", 
      shape   = 23, 
      size    = 3) + 
    # Prettify x axis...
    scale_x_continuous(
      breaks = x_major$x_major, 
      labels = x_major$d_v_a) +
    # Prettify y axis...
    scale_y_continuous(
      name   = "Impact per fully vaccinated person (log10 scale)",
      trans  = "log10",
      labels = scientific, 
      limits = c(10 ^ lb, 10 ^ ub),
      expand = c(0, 0),
      breaks = 10 ^ rev(ub : lb))
  
  # ---- Prettify and save ----
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.text     = element_text(size = 10),
          axis.text.x   = element_text(hjust = 1, angle = 50),
          axis.title.x  = element_blank(),
          axis.title.y  = element_text(
            size = 16, margin = margin(l = 10, r = 20)),
          axis.line     = element_blank(),
          panel.grid.major.y = element_line(linewidth = 0.25),
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(1, "lines"))
  
  # Filename to save to depending on scope
  save_stem  = "Impact - coverage by income and decade"
  save_scope = str_replace(scope, "_", " ")
  
  # Save these figures to file
  save_fig(g, save_stem, save_scope, dir = "impact_functions")
  
  # Save all time plot as key manuscript figure
  # if (scope == "all_time")
  #   save_fig(g, "Figure 3", dir = "manuscript")
}

# **** Historical impact ****

# ---------------------------------------------------------
# Main results plot - historical impact over time
# ---------------------------------------------------------
plot_historical_impact = function(region = NULL) {
  
  # ---- Figure properties ----
  
  # Dictionary for temporal and cumulative subplots
  metrics = qc(deaths, yll, dalys)
  
  # Set scale for set of metrics
  scale = list(
    deaths = list(val = 1e6, str = "million"),
    yll    = list(val = 1e9, str = "billion"), 
    dalys  = list(val = 1e9, str = "billion"))
  
  # Diseases to combine into one colour (set to NULL to turn off)
  grouping = qc(dip, hepb, je, mena, pcv, rota, rubella, yf) 
  
  # Descriptive name for this grouping
  other = "Other pathogens"
  
  # Custom y axes limits for set of metrics
  y_lims = c(200, 11, 11)
  
  # Year range string
  range = paste(range(o$years), collapse = "-")
  
  # ---- Regional scope ----
  
  # Regional scope
  if (!is.null(region)) {
    scope = table("region_dict") %>%
      filter(region == !!region) %>%
      pull(region_name)
  }
  
  # Global scope
  if (is.null(region)) {
    region = all_regions()
    scope  = NULL
  }
  
  # Display what we're plotting to user
  msg = "  - Plotting historical impact"
  message(paste(c(msg, scope), collapse = ": "))
  
  # ---- Load results ----
  
  # Function to load results for a given metric
  load_fn = function(metric) {
    
    # Load result and append metric ID
    impact_dt = read_rds("history", "burden_averted", metric) %>%
      mutate(metric = !!metric)
    
    return(impact_dt)
  }
  
  # Load results and apply initial formatting
  impact_dt = lapply(metrics, load_fn) %>%
    rbindlist() %>%
    # Remove 'other' pathogens - these go last...
    left_join(y  = table("d_v_a"), 
              by = "d_v_a_id") %>%
    left_join(y  = table("disease_name"), 
              by = "disease") %>%
    # Tidy up...
    select(metric, disease, disease_name, 
           country, year, impact)
  
  # ---- Summarize results ----
  
  # Cumulative impact by metric and disease group
  results_dt = impact_dt %>%
    # Subset to regional scope
    left_join(y  = table("country"), 
              by = "country") %>%
    filter(region %in% !!region) %>%
    # Combine subset of diseases...
    mutate(group = ifelse(
      test = disease %in% grouping, 
      yes  = other, 
      no   = disease_name)) %>%
    # Cumulative results for each disease...
    group_by(metric, group, year) %>%
    summarise(impact = sum(impact)) %>%
    mutate(impact = cumsum(impact)) %>%
    ungroup() %>%
    rename(disease = group) %>%
    # All metrics bounded above by DALYs...
    group_by(disease, year) %>%
    mutate(impact = pmin(
      impact, impact[metric == "dalys"])) %>%
    ungroup() %>%
    as.data.table()

  # ---- Disease totals ----
  
  # Format scaling details into joinable datatable
  scale_dt = list2dt(scale) %>%
    mutate(metric = names(scale), .before = 1)
  
  # Construct labels: total FVPs over analysis timeframe
  label_dt = results_dt %>%
    # Total burden averted over all years...
    group_by(disease, metric) %>%
    summarise(total = max(impact)) %>%
    ungroup() %>%
    # Apply metric scaling...
    left_join(y  = scale_dt, 
              by = "metric") %>%
    mutate(total = paste(round(total / val, 1), str)) %>%
    select(-val, -str) %>%
    # Append metric and total strings...
    left_join(y  = table("metric_dict"), 
              by = "metric") %>%
    mutate(total = paste0(metric_short, ": ", total)) %>%
    select(-metric_short, -metric_long, -metric_impact) %>%
    # Format into single string per disease...
    pivot_wider(names_from  = metric,
                values_from = total) %>%
    mutate(name = paste0("**", disease, "**"), 
           .after = disease) %>%
    select(disease, name, all_of(metrics)) %>%
    # Collapse into single string with line breaks...
    unite(col = "label", names(.)[-1], 
          sep = "<br>", remove = FALSE) %>%
    select(disease, label) %>%
    as.data.table()
  
  # ---- Disease colours ----
  
  # Common plotting order of diseases based on global deaths averted
  diseases = impact_dt %>%
    filter(metric == metrics[1], 
           !disease %in% grouping) %>%
    # Order by decreasing global impact...
    group_by(disease_name) %>%
    summarise(order = sum(impact)) %>%
    ungroup() %>%
    arrange(-order) %>%
    # Covert to vector
    pull(disease_name) %>%
    c(other)
  
  # Set colours as named vector for consistent colouring
  colours = colours_who("category", length(diseases))
  
  # ---- Construct plotting datatables ----
  
  # Set order and append total labels to plotting data
  plot_dt = results_dt %>%
    # Apply scaler...
    left_join(y  = scale_dt,
              by = "metric") %>%
    mutate(value = impact / val) %>%
    # Full metric description...
    left_join(y  = table("metric_dict"), 
              by = "metric") %>%
    mutate(metric_id = factor(metric, metrics), 
           metric = paste0(
             metric_impact, "\n(cumulative ", 
             range, ", in ", str, "s)")) %>%
    # Append labels...
    left_join(y  = label_dt, 
              by = "disease") %>%
    # Plotting order...
    mutate(disease = factor(disease, diseases)) %>%
    arrange(metric_id, disease) %>%
    mutate(metric = fct_inorder(metric), 
           label  = fct_inorder(label)) %>%
    select(label, metric_id, metric, year, value)
  
  # Blank plot to get custom y axes limits
  blank_dt = plot_dt %>%
    select(metric) %>%
    unique() %>%
    mutate(year  = min(o$years), 
           value = y_lims)
  
  # We'll also plot metric totals in each facet
  text_dt = plot_dt %>%
    # Total impact by metric...
    filter(year %in% max(o$years)) %>%
    group_by(metric_id, metric) %>%
    summarise(value = round(sum(value), 1)) %>%
    ungroup() %>%
    # Construct total impact text...
    left_join(y  = scale_dt, 
              by = c("metric_id" = "metric")) %>%
    mutate(text = paste("Total:", value, str)) %>%
    select(metric, text) %>%
    # Set facet positioning...
    left_join(blank_dt, by = "metric") %>%
    mutate(value = 0.95 * value, 
           year  = 0.1 * diff(range(o$years)) + year) %>%
    select(metric, year, value, text)
  
  # ---- Produce plot ----
  
  # Stacked yearly bar plot
  g = ggplot(plot_dt) +
    aes(x = year, 
        y = value) + 
    # Plot bars...
    geom_col(
      mapping = aes(
        fill = label)) + 
    # Add total text...
    geom_text(
      data    = text_dt, 
      mapping = aes(
        label = text), 
      hjust   = "left", 
      size    = 4) + 
    # Facet by temporal-cumulative metric...
    facet_wrap(
      facets = vars(metric), 
      scales = "free_y", 
      nrow   = 1) + 
    # Set colours...
    scale_fill_manual(values = colours) + 
    # Prettify y axis...
    geom_blank(data = blank_dt) + 
    scale_y_continuous(
      labels = comma, 
      breaks = pretty_breaks(),
      expand = expansion(mult = c(0, 0))) +
    # Prettiy x axis...
    scale_x_continuous(
      limits = c(min(o$years) - 1, 
                 max(o$years) + 1), 
      expand = expansion(mult = c(0, 0)), 
      breaks = seq(
        from = min(o$years), 
        to   = max(o$years), 
        by   = 5)) +
    # Prettify legend (needed for y spacing to take effect)...
    guides(fill = guide_legend(
      byrow = TRUE))
  
  # Plot region name unless plotting global results
  if (!is.null(scope))
    g = g + ggtitle(scope)
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.title    = element_blank(),
          axis.text     = element_text(size = 9),
          axis.text.x   = element_text(hjust = 1, angle = 50), 
          axis.line     = element_blank(),
          plot.title    = element_text(
            margin = margin(t = 10, b = 20), 
            size   = 18,
            hjust  = 0.5), 
          strip.text    = element_text(size = 12),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(0.5, "lines"),
          panel.grid.major.y = element_line(linewidth = 0.25),
          legend.title  = element_blank(),
          legend.text   = element_markdown(size = 11),
          legend.key    = element_blank(),
          legend.position = "right", 
          legend.spacing.y  = unit(2, "lines"),
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines"))
  
  # Save these figures to file
  save_fig(g, "Historical impact", scope, 
           dir = "historical_impact")
  
  # Save global results as main manuscript figure
  if (is.null(scope))
  save_fig(g, "Figure 1", dir = "manuscript")
}

# ---------------------------------------------------------
# Main results plot - historical impact over time
# ---------------------------------------------------------
plot_temporal_impact = function(metric) {
  
  message("  - Plotting temporal impact: ", metric)
  
  # Results by disease and region
  plot_dt = read_rds("history", "all_samples", metric) %>%
    # Append full region and disease names...
    append_region_name() %>%
    left_join(y  = table("d_v_a"), 
              by = "d_v_a_id") %>%
    left_join(y  = table("disease_name"), 
              by = "disease") %>%
    select(disease = disease_name, region, year, sample, impact) %>%
    # Results for each disease and each region...
    group_by(disease, region, year, sample) %>%
    summarise(impact = sum(impact)) %>%
    ungroup() %>%
    # Remove only single data points...
    add_count(disease, region) %>%
    mutate(n = n / (o$uncertainty_samples + 1)) %>%
    filter(n > 1) %>%
    select(-n) %>%
    # Summarising uncertainty...
    summarise_uncertainty(cumulative = FALSE) %>%
    as.data.table()
  
  # Colour scheme
  colours = colours_who_region() %>%
    enframe() %>%
    left_join(y  = table("region_dict"), 
              by = c("name" = "region")) %>%
    select(region_name, value) %>%
    deframe()
  
  # Y axis label describing metric
  axis_label = table("metric_dict") %>%
    filter(metric == !!metric) %>%
    pull(metric_impact)
  
  # Stacked yearly bar plot
  g = ggplot(plot_dt) +
    aes(x = year, 
        y = impact) + 
    # Plot uncertainty bands...
    geom_ribbon(
      mapping = aes(
        ymin = lower,
        ymax = upper,
        fill = region),
      colour = NA,
      alpha  = 0.5) +
    # Plot best prediction line...
    geom_line(
      mapping = aes(colour = region), 
      linewidth = 1.2) + 
    # Facet by temporal-cumulative metric...
    facet_wrap(
      facets = vars(disease), 
      scales = "free_y") + 
    # Set colours...
    scale_colour_manual(
      values = colours) + 
    # Prettify y axis...
    scale_y_continuous(
      name   = axis_label,
      labels = comma, 
      breaks = pretty_breaks(),
      expand = expansion(mult = c(0, 0.05))) +
    # Prettiy x axis...
    scale_x_continuous(
      limits = c(min(o$years), 
                 max(o$years)), 
      expand = expansion(mult = c(0, 0)), 
      breaks = seq(
        from = min(o$years), 
        to   = max(o$years), 
        by   = 5))
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.title.x  = element_blank(),
          axis.title.y  = element_text(size = 20),
          axis.text     = element_text(size = 9),
          axis.text.x   = element_text(hjust = 1, angle = 50), 
          axis.line     = element_blank(),
          strip.text    = element_text(size = 12),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(0.5, "lines"),
          panel.grid.major.y = element_line(linewidth = 0.25),
          legend.title  = element_blank(),
          legend.text   = element_text(size = 12),
          legend.key    = element_blank(),
          legend.position = "bottom", 
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines"))
  
  # Save these figures to file
  save_fig(g, "Temporal impact by region", metric, 
           dir = "historical_impact")
}

# ---------------------------------------------------------
# Plot change in infant mortality rates over time
# ---------------------------------------------------------
plot_infant_mortality = function() {
  
  message("  - Plotting infant mortality rates")
  
  # ---- Figure properties ----
  
  age_bound = 0
  
  # Set colour scheme for vaccine coverage
  colour_map = "viridis::viridis" # pals::kovesi.rainbow "viridis::viridis"
  
  year_range = paste(range(o$years), collapse = "-")
  
  metric_dict = list(
    rate   = paste("Infant mortality rate", year_range), 
    deaths = paste("Cumulative infant deaths", year_range), 
    cov    = "Global vaccine coverage")
  
  metric_lab = list(
    rate   = "percent", 
    deaths = "comma") 
  
  # Dictionary for each scenario
  line_dict = list(
    vaccine    = "Vaccination as obserevd",
    no_vaccine = "Hypothetical scenario:<br>no historical vaccination",
    no_other   = "Hypothetical scenario:<br>no improvement in infant survival")
  
  # Dictionary for each scenario
  fill_dict = list( 
    vaccine = "Reduction in infant mortality<br>attributable to vaccination", 
    other   = "Reduction in infant mortality<br>attributable to other causes")
  
  fill_colours = c(colours_who("logo", 1), "grey30")
  
  # ---- Compile results ----
  
  # Load mortality rates in vaccine and no vaccine scenarios
  mortality = mortality_rates(
    age_bound = age_bound, 
    grouping  = "none")
  
  # Mortality rate over time in each scenario
  rate_dt = mortality$rate %>%
    mutate(metric = "rate") %>%
    select(metric, scenario, year, value) %>%
    as.data.table()
  
  # Cumulative number of deaths in each scenario
  deaths_dt = mortality$value %>%
    group_by(scenario) %>%
    mutate(value = cumsum(value)) %>%
    ungroup() %>%
    mutate(metric = "deaths") %>%
    select(metric, scenario, year, value) %>%
    as.data.table()
  
  # Vaccine coverage over time - global average
  coverage_dt = table("coverage_global") %>%
    filter(coverage > 0,
           !grepl("_BD$", vaccine)) %>%
    # Combine DTP3 predictions...
    mutate(vaccine = ifelse(
      test = vaccine %in% qc(dip, tet, aper),
      yes  = "DTP3",
      no   = vaccine)) %>%
    group_by(vaccine, year) %>%
    summarise(coverage = mean(coverage)) %>%
    ungroup() %>%
    # Append descriptive vaccine names...
    left_join(y  = table("vaccine_name"),
              by = "vaccine") %>%
    select(vaccine = vaccine_name, year, coverage) %>%
    replace_na(list(vaccine = "DTP third dose")) %>%
    # Metric subplot order...
    mutate(metric = metric_dict$cov,
           metric = factor(metric, metric_dict)) %>%
    # Tidy up...
    arrange(vaccine, year) %>%
    as.data.table()
  
  # ---- Impact attributable to vaccination ----
  
  # Percentage of mortality decline attributable to vaccination
  attributable_dt = rate_dt %>%
    mutate(value = value * 100) %>%
    pivot_wider(names_from = scenario) %>%
    mutate(relative = 100 * (no_vaccine - vaccine) / 
             (no_vaccine[1] - vaccine)) %>%
    filter(year > o$years[1]) %>%
    mutate(n_years = 1 : n(), 
           relative_mean = cumsum(relative) / n_years) %>%
    as.data.table()
  
  attributable = round(mean(attributable_dt$relative))
  
  # Append this very important outcome to legend labels
  fill_dict$vaccine %<>% paste0(": **", attributable, "%**")
  fill_dict$other   %<>% paste0(": **", 100 - attributable, "%**")
  
  # ---- Construct plotting datatables ----
  
  # Combine into single plotting datatable
  lines_dt = rate_dt %>%
    bind_rows(deaths_dt) %>%
    # Smaller linewidth for hypothetical scenarios...
    mutate(size = ifelse(scenario == "vaccine", 2, 1), 
           .after = scenario) %>%
    # Scenario order...
    rename(line = scenario) %>%
    mutate(line = recode(line, !!!line_dict), 
           line = factor(line, line_dict)) %>%
    # Metric subplot order...
    mutate(metric = recode(metric, !!!metric_dict),
           metric = factor(metric, metric_dict)) %>%
    arrange(metric, line, year)
  
  # Second plotting datatable to highlight vaccine contribution
  area_dt = rate_dt %>%
    rbind(deaths_dt) %>%
    pivot_wider(names_from = scenario) %>%
    # Fill impact attributable to vaccination and other factors...
    mutate(y0.vaccine = pmin(vaccine, no_vaccine),
           y1.vaccine = pmax(vaccine, no_vaccine), 
           y0.other   = pmin(no_vaccine, no_other),
           y1.other   = pmax(no_vaccine, no_other)) %>%
    select(-vaccine, -no_vaccine, -no_other) %>%
    # Format into plottable datatable...
    pivot_longer(cols = -c(metric, year)) %>%
    separate(col  = "name", 
             into = c("var", "fill"), 
             sep  = "\\.") %>%
    pivot_wider(names_from = var) %>%
    # Scenario order...
    mutate(fill = recode(fill, !!!fill_dict), 
           fill = factor(fill, fill_dict)) %>%
    # Metric subplot order...
    mutate(metric = recode(metric, !!!metric_dict),
           metric = factor(metric, metric_dict)) %>%
    # Tidy up...
    select(metric, fill, year, y0, y1) %>%
    arrange(metric, fill, year) %>%
    as.data.table()
  
  # ---- Produce plot ----
  
  # Function to prettify theme
  theme_fn = function(g) {
    
    # Prettify theme
    g = g + theme_classic() + 
      theme(axis.title    = element_blank(),
            axis.text     = element_text(size = 10),
            axis.text.x   = element_text(hjust = 1, angle = 50), 
            axis.line     = element_blank(),
            strip.text    = element_text(size = 16),
            strip.background = element_blank(), 
            panel.border  = element_rect(
              linewidth = 0.5, fill = NA),
            panel.spacing = unit(1, "lines"),
            legend.title  = element_blank(),
            legend.text   = element_markdown(size = 11),
            legend.position   = "right",
            legend.spacing.y  = unit(0.5, "lines"),
            legend.key.height = unit(1.1, "lines"),
            legend.key.width  = unit(1.8, "lines"))
    
    return(g)
  }
  
  # Function to create subplot
  plot_fn = function(plot_metric) {
    
    sub = metric_dict[[plot_metric]]
    
    # Plot infant mortality
    g = ggplot(area_dt[metric == sub]) +
      aes(x = year) +
      # Plot area...
      geom_ribbon(
        mapping = aes(
          ymin = y0,
          ymax = y1, 
          fill = fill),
        colour = NA, 
        alpha  = 0.5) +
      # Plot lines...
      geom_line(
        data    = lines_dt[metric == sub],
        mapping = aes(
          y = value,
          linetype  = line,
          linewidth = size), 
        colour = "black") +
      # Facet by metric...
      facet_wrap(
        facets = vars(metric), 
        scales = "free_y", 
        ncol   = 1) +
      # Set colour scheme...
      scale_fill_manual(
        values = fill_colours) +
      # Set linewidth values
      scale_linewidth(
        range = c(0.5, 1.2)) +
      # Prettify y axis...
      scale_y_continuous(
        labels = get(metric_lab[[plot_metric]]), 
        limits = c(0, NA),
        expand = expansion(mult = c(0, 0.05)), 
        breaks = pretty_breaks()) +
      # Prettiy x axis...
      scale_x_continuous(
        limits = c(min(o$years), max(o$years)), 
        expand = expansion(mult = c(0, 0)), 
        breaks = seq(min(o$years), max(o$years), by = 5)) +
      # Prettify legends...
      guides(
        linewidth = "none", 
        linetype = guide_legend(
          order     = 1, 
          keyheight = 1.8, 
          reverse   = TRUE),
        fill = guide_legend(
          order     = 2, 
          keyheight = 1.8, 
          reverse   = TRUE))
    
    g = theme_fn(g)
    
    return(g)
  }
  
  # Plot coverage over time
  g_cov = ggplot(coverage_dt) +
    aes(x = year, 
        y = coverage, 
        colour = vaccine) +
    # Plot lines...
    geom_line(linewidth = 1.2) +
    # Trivial facet...
    facet_wrap(
      facets = vars(metric)) +  
    # Set colour scheme...
    scale_color_manual(
      values = rev(colour_scheme(
        map = colour_map, 
        n   = n_unique(coverage_dt$vaccine)))) +
    # Set linewidth values
    scale_linewidth(
      range = c(0.5, 1.2)) +
    # Prettify y axis...
    scale_y_continuous(
      labels = percent,
      limits = c(0, 1),
      expand = c(0, 0),
      breaks = pretty_breaks()) + 
    # Prettiy x axis...
    scale_x_continuous(
      limits = c(min(o$years), max(o$years)), 
      expand = expansion(mult = c(0, 0)), 
      breaks = seq(min(o$years), max(o$years), by = 5))
  
  # Create subplots
  g1 = plot_fn("rate") 
  g2 = plot_fn("deaths") + theme(legend.position = "none")
  g3 = theme_fn(g_cov)
  
  # Patch subplots together
  g = g1 / g2 / g3

  # Save figure to file
  save_fig(g, "Infant mortality rates", dir = "historical_impact")
  
  # Also save as main manuscript figure
  save_fig(g, "Figure 2", dir = "manuscript")
}

# ---------------------------------------------------------
# Plot vaccine contribution to infant mortality decrease by region
# ---------------------------------------------------------
plot_mortality_region = function() {
  
  message("  - Plotting infant mortality decrease by region")
  
  # Load mortality rates in vaccine and no vaccine scenarios
  region_dt = mortality_rates(grouping = "region")$rate
  global_dt = mortality_rates(grouping = "none")$rate
  
  plot_dt = global_dt %>%
    mutate(group = "World") %>%
    rbind(region_dt) %>%
    pivot_wider(names_from = scenario) %>%
    group_by(group) %>%
    mutate(relative = (no_vaccine - vaccine) / 
             (no_vaccine[1] - vaccine), 
           n_years = 1 : n(), 
           contribution = cumsum(relative) / n_years) %>%
    ungroup() %>%
    mutate(contribution = pmin(contribution, 1)) %>%
    filter(year >= 1980) %>%
    select(group, year, contribution) %>%
    as.data.table()
  
  # Plot all metrics over time
  g = ggplot(plot_dt) +
    aes(x = year, 
        y = contribution, 
        colour = group) +
    geom_line() 
}

# ---------------------------------------------------------
# Plot regional differences in infant mortality changes
# ---------------------------------------------------------
plot_mortality_change = function() {
  
  message("  - Plotting regional changes in infant mortality")
  
  # Description of metric and year scope
  metric_str = "infant mortality rate"
  years_str  = paste0("(", paste(range(o$years), collapse = "-"), ")")
  
  # Dictionary for metric type
  type_dict = list(
    abs = paste("Absolute decrease in", metric_str, years_str),
    rel = paste("Relative decrease in", metric_str, years_str), 
    att = paste("Contribution of vaccination to decrease in", 
                metric_str, years_str))
  
  # Metric type colour scheme (+1 for global average)
  colours = c("#EBAA2D", "#DF721F", "#71C2A9", "#808080")
  
  # Mortality rates - by region and global average
  region_dt = mortality_rates(grouping = "region")$rate
  world_dt  = mortality_rates(grouping = "none")$rate %>%
    mutate(group = "World")
  
  abs_dt = rbind(region_dt, world_dt) %>%
    filter(scenario == "vaccine", 
           year %in% range(o$years)) %>%
    group_by(group) %>%
    summarise(start = value[1], 
              end   = value[2], 
              abs   = start - end) %>%
    ungroup() %>%
    arrange(-abs) %>%
    as.data.table()
  
  rel_dt = rbind(region_dt, world_dt) %>%
    filter(scenario == "vaccine", 
           year %in% range(o$years)) %>%
    group_by(group) %>%
    summarise(rel = 1 - min(value) / max(value)) %>%
    ungroup() %>%
    as.data.table()
  
  # Percentage of mortality decline attributable to vaccination
  att_dt = rbind(region_dt, world_dt) %>%
    pivot_wider(names_from = scenario) %>%
    select(group, year, vaccine, no_vaccine) %>%
    group_by(group) %>%
    mutate(att = (no_vaccine - vaccine) / 
             (no_vaccine[1] - vaccine)) %>%
    filter(year > o$years[1]) %>%
    summarise(att = mean(att)) %>%
    ungroup() %>%
    as.data.table()
  
  # Order regions by absolute decrease in mortality rates
  order = abs_dt %>%
    filter(group != "World") %>%
    arrange(desc(abs)) %>%
    pull(group) %>%
    c("World")
  
  # Labels bars for clarity
  label_dt = abs_dt %>%
    # Convert to percentage...
    mutate(across(
      .cols = c(abs, start, end), 
      .fns  = ~ . * 100)) %>%
    # Construct string explaining start and end values...
    mutate(label = sprintf(
      "%.1f%%\n(%.1f%% to %.1f%%)", 
      abs, start, end)) %>%
    # Only required for absolute bars...
    mutate(type = "abs") %>%
    select(group, type, label)
  
  # Construct plotting datatable
  plot_dt = abs_dt %>%
    select(-start, -end) %>%
    # Append all results...
    left_join(rel_dt, by = "group") %>%
    left_join(att_dt, by = "group") %>%
    pivot_longer(cols = -group, 
                 names_to = "type") %>%
    # Labels: simple and expansive versions...
    mutate(str = sprintf("%.0f%%", value * 100)) %>%
    left_join(y  = label_dt, 
              by = c("group", "type")) %>%
    mutate(label = ifelse(is.na(label), str, label)) %>%
    # Position according to metric type...
    group_by(type) %>%
    mutate(nudge   = max(value) * 0.05, 
           lab_pos = value + nudge) %>%
    ungroup() %>%
    select(-str, -nudge) %>%
    # Define colour groups (world is special case)...
    mutate(colour = match(type, names(type_dict)), 
           colour = ifelse(
             test = group == "World", 
             yes  = max(colour) + 1, 
             no   = colour), 
           colour = as.character(colour)) %>%
    # Set factors for desired ordering...
    mutate(type  = recode(type, !!!type_dict), 
           type  = factor(type, type_dict), 
           group = factor(group, order)) %>%
    as.data.table()
  
  # Plot regional results as bars
  g = ggplot(plot_dt) +
    aes(x = group, 
        y = value,  
        fill = colour) +
    # Plot bars...
    geom_bar(
      stat  = "identity", 
      width = 0.6) +
    # Add descriptive text...
    geom_text(
      mapping = aes(
        label = label, 
        y     = lab_pos), 
      vjust   = 0, 
      size    = 3.5) + 
    # Facet by type...
    facet_wrap(
      facets = vars(type), 
      ncol   = 1, 
      scales = "free") + 
    # Set colours...
    scale_fill_manual(values = colours) +
    # Prettify y axis...
    scale_y_continuous(
      labels = percent, 
      limits = c(0, NA), 
      expand = expansion(mult = c(0, 0.2)), 
      breaks = pretty_breaks())
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.title    = element_blank(),
          axis.text.x   = element_text(size = 14),
          axis.text.y   = element_text(size = 10),
          axis.line     = element_blank(),
          strip.text    = element_text(size = 20),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(2, "lines"),
          panel.grid.major.y = element_line(linewidth = 0.25),
          legend.position = "none")
  
  # Save figure to file
  save_fig(g, "Infant mortality change by region", 
           dir = "historical_impact")
  
  # Also save as main manuscript figure
  # save_fig(g, "Figure 3", dir = "manuscript")
}

# ---------------------------------------------------------
# Plot absolute and relative probability of death in 2024
# ---------------------------------------------------------
plot_prob_death_age = function() {
  
  message("  - Plotting probability of death by age")
  
  age_bins = 4
  
  # Dictionary for each vaccine case
  #
  # NOTE: Using reverse order here for prettier plot
  case_dict = list( 
    no_vaccine = "No historical vaccination",
    vaccine    = "Vaccination as obserevd")
  
  # Construct a somewhat elaborate y axis label
  y_label = bquote(
    ~.("Marginal probability of death in n")
    ^.("th") 
    ~.("year of life")
    ~.(paste0("(", max(o$years), ")")))
  
  # Deaths averted by vaccination
  averted_dt = read_rds("history", "burden_averted_deaths") %>%
    filter(year == max(o$years)) %>%
    # Apply age structure...
    left_join(y  = table("impact_age_multiplier"), 
              by = "d_v_a_id", 
              relationship = "many-to-many") %>%
    mutate(impact = impact * scaler) %>%
    # Summarise over region...
    append_region_name() %>%
    group_by(region, age) %>%
    summarise(averted = sum(impact)) %>%
    ungroup() %>%
    # Avoid zero for log scaling ...
    mutate(age = age + 1) %>%
    filter(age %in% 2 ^ (0 : age_bins)) %>%
    select(region, age, averted) %>%
    as.data.table()
  
  # Probability of death by age
  plot_dt = table("wpp_pop") %>%
    filter(year == max(o$years)) %>%
    left_join(y  = table("wpp_deaths"), 
              by = c("country", "year", "age")) %>%
    # Avoid zero for log scaling ...
    mutate(age = age + 1) %>%
    filter(age %in% 2 ^ (0 : age_bins)) %>%
    # Average prob of death by region and year...
    append_region_name() %>%
    group_by(region, age) %>%
    summarise(pop    = sum(pop), 
              deaths = sum(deaths)) %>%
    ungroup() %>%
    # Append deaths averted...
    left_join(y  = averted_dt, 
              by = c("region", "age")) %>%
    replace_na(list(averted = 0)) %>%
    mutate(vaccine    = deaths / pop, 
           no_vaccine = (deaths + averted) / pop) %>%
    # Melt to long format ready for plotting...
    select(region, age, vaccine, no_vaccine) %>%
    pivot_longer(cols = c(vaccine, no_vaccine), 
                 names_to = "case") %>%
    # Set regional order...
    arrange(desc(value)) %>%
    mutate(region = fct_inorder(region)) %>%
    # Set vaccine case order...
    mutate(case = recode(case, !!!case_dict), 
           case = factor(case, case_dict)) %>%
    arrange(region, case, age) %>%
    as.data.table()
  
  # Plot both scenarios as side-by-side bars
  g = ggplot(plot_dt) +
    aes(x = age, 
        y = value, 
        fill = case) +
    geom_col(position = "dodge") +
    # Facet by region...
    facet_wrap(
      facets = vars(region), 
      scales = "free_y") +
    # Set colour scheme...
    scale_fill_manual(
      values = colours_who("logo", n = 2)) + 
    # Prettify x axis...
    scale_x_continuous(
      name   = "Age (log2 scale)",
      trans  = "log2", 
      breaks = 2 ^ (0 : age_bins)) +
    # Prettify y axis...
    scale_y_continuous(
      name   = y_label,
      labels = percent, 
      limits = c(0, NA), 
      expand = expansion(mult = c(0, 0.05)), 
      breaks = pretty_breaks())
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.title   = element_text(size = 20),
          axis.title.x = element_text(
            margin = margin(t = 20, b = 10)),
          axis.title.y = element_text(
            margin = margin(l = 10, r = 20)),
          axis.text    = element_text(size = 12),
          axis.ticks   = element_blank(), 
          axis.line    = element_line(linewidth = 0.25),
          strip.text   = element_text(size = 18),
          strip.background = element_blank(), 
          panel.spacing = unit(1, "lines"), 
          panel.grid.major.y = element_line(linewidth = 0.25),
          legend.title = element_blank(),
          legend.text  = element_text(size = 16),
          legend.key   = element_blank(),
          legend.position = "bottom", 
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines"))
  
  # Save figure to file
  save_fig(g, "Probability of death by age", dir = "historical_impact")
}

# ---------------------------------------------------------
# Plot absolute and relative probability of death in 2024
# ---------------------------------------------------------
plot_survival_increase = function(log_age = FALSE) {
  
  message("  - Plotting increase in childhood survival")
  
  # Year to create snapshot for - makes sence to be final year of analysis
  snapshot_year = max(o$years)
  
  # Age bound (nice round value or a square)
  age_bound = ifelse(log_age, 64, 50)
  
  # Construct pplot title
  title = "Historical vaccination compared to hypothetical no vaccination"

  # Construct a somewhat elaborate y axis label
  y_label = paste0(
    "Marginal increase in survival probability in n<sup>th</sup> ", 
    "year of life (", max(o$years), " snapshot)")
  
  # Dictionary for absolute and marginal metrics
  metric_dict = c(
    relative = "Relative", 
    absolute = "Absolute")
  
  # predictiond child deaths averted by vaccination
  averted_dt = read_rds("history", "burden_averted_deaths") %>%
    # Results in year of interest...
    filter(year == snapshot_year) %>%
    # Apply age effect...
    left_join(y  = table("impact_age_multiplier"), 
              by = "d_v_a_id", 
              relationship = "many-to-many") %>%
    mutate(impact = impact * scaler) %>%
    # Summarise results over region...
    append_region_name() %>%
    group_by(region, age) %>%
    summarise(averted = sum(impact)) %>%
    ungroup() %>%
    as.data.table()
  
  # Calculate absolute and relative effect in year of interest
  survival_dt = table("wpp_pop") %>%
    left_join(y  = table("wpp_deaths"), 
              by = c("country", "year", "age")) %>%
    # Retain only years and ages of interest...
    filter(year == snapshot_year, 
           age  <= age_bound - log_age * 1) %>%
    # Append region...
    append_region_name() %>%
    # Average prob of death by region and year...
    group_by(region, age) %>%
    summarise(pop    = sum(pop), 
              deaths = sum(deaths)) %>%
    ungroup() %>%
    # Append deaths averted results...
    left_join(y  = averted_dt, 
              by = c("region", "age")) %>%
    mutate(age = age + log_age * 1) %>%
    # Cumulative results...
    group_by(region) %>%
    mutate(c_pop   = cumsum(pop), 
           c_death = cumsum(deaths), 
           c_avert = cumsum(averted)) %>%
    ungroup() %>%
    # Compute absolute and relative results...
    mutate(vaccine    = c_death / c_pop, 
           no_vaccine = (c_death + c_avert) / c_pop, 
           absolute   = no_vaccine - vaccine,
           relative   = 1 - vaccine / no_vaccine) %>%
    # Melt to long format...
    select(region, age, pop, absolute, relative) %>%
    pivot_longer(cols = c(absolute, relative), 
                 names_to = "metric") %>%
    # Tidy up...
    select(metric, region, age, pop, value) %>%
    arrange(metric, region, age) %>%
    as.data.table()
  
  # Combine with global values
  plot_dt = survival_dt %>%
    # Calculate global average...
    group_by(metric, age) %>%
    mutate(weight = pop / sum(pop)) %>%
    summarise(value = sum(value * weight)) %>%
    ungroup() %>%
    # Append regional results...
    mutate(region = "World", 
           size   = 2) %>%
    bind_rows(survival_dt) %>%
    select(-pop) %>%
    replace_na(list(size = 1)) %>%
    # Metric order...
    mutate(metric = recode(metric, !!!metric_dict), 
           metric = factor(metric, metric_dict)) %>%
    # Tidy up...
    arrange(size, region) %>%
    mutate(region = fct_inorder(region)) %>%
    as.data.table()
  
  # ---- Produce plot ----
  
  # WHO regional colours, with global average in black
  colours = unname(c(colours_who_region(), "#000000"))
  
  # Plot absolute and relative difference in snapshot year
  g = ggplot(plot_dt) +
    aes(x = age, 
        y = value, 
        colour = fct_rev(region), 
        linewidth = size) +
    geom_line() +
    # Facet by absolute/relative...
    facet_wrap(
      facets = vars(metric), 
      scales = "free_y", 
      ncol   = 1) + 
    # Set colour scheme
    scale_colour_manual(
      values = rev(colours)) +
    scale_linewidth(
      range = c(1.2, 3)) +
    # Add plot title...
    ggtitle(title) +
    # Prettify y axis...
    scale_y_continuous(
      name   = y_label, 
      labels = percent, 
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05)), 
      breaks = pretty_breaks()) +  
    # Prettify legend...
    guides(
      linetype  = "none", 
      linewidth = "none", 
      colour = guide_legend(
        reverse = TRUE, 
        override.aes = list(
          linewidth = 2)))
  
  # Transform x axis to log2 scale
  if (log_age == TRUE) {
    g = g + scale_x_continuous(
      name   = "Age (log2 scale)",
      trans  = "log2", 
      breaks = 2 ^ (1 : (sqrt(age_bound) - 1) - 1))
  }
  
  # Keep x axis on real number line
  if (log_age == FALSE) {
    g = g + scale_x_continuous(
      name   = "Age in years", 
      expand = expansion(add = 1),
      breaks = seq(0, age_bound, by = 5))
  }
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.title   = element_text(size = 16),
          axis.title.x = element_text(
            margin = margin(t = 20, b = 10)),
          axis.title.y = element_markdown(
            margin = margin(l = 10, r = 20)),
          axis.text = element_text(size = 11),
          axis.line  = element_blank(),
          strip.text = element_text(size = 14),
          strip.background = element_blank(), 
          plot.title    = element_text(
            margin = margin(t = 10, b = 20), 
            size   = 18,
            hjust  = 0.5), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(1, "lines"),
          legend.title  = element_blank(),
          legend.text   = element_text(size = 13),
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines"))
  
  # Save figure to file
  save_fig(g, "Increase in survival", dir = "historical_impact")

  # Also save as main manuscript figure
  save_fig(g, "Figure 3", dir = "manuscript")
}

# ---------------------------------------------------------
# Plot measles deaths in context of all cause deaths
# ---------------------------------------------------------
plot_measles_in_context = function() {
  
  message("  - Plotting measles in context")
  
  # Upper age bound for predictions
  age_bound = 5
  
  # Metrics to plot
  metric_dict = list(
    cum  = "Cumulative number of deaths in children under %i",
    abs  = "Annual number of deaths in children under %i",
    norm = "Normalised cause of death for children under %i",
    cov  = "Measles vaccine coverage (%i-%i)")
  
  # Cause of death
  cause_dict = list(
    all_cause = "All cause deaths",
    measles   = "Measles as cause of death", 
    other     = "Other cause of death")
  
  # Define vaccine scenarios
  case_dict = list( 
    vaccine    = "Vaccination as obserevd", 
    no_vaccine = "No historical vaccination")
  
  # Name of d-v-a to plot
  measles_id = "Measles"
  
  # ID of coverage values to plot
  coverage_ids = c(101, 102)
  
  # Colour schemes
  colours = list(
    area  = c("#B0B0B0", "#0189C5"), 
    line  = c("#000000", "#0189C5"), 
    doses = c("#C70A0A", "#901616", "#400808"))
  
  # ---- Deaths by cause ----
  
  # Death data - WPP and external measles models
  deaths = list(
    wpp    = table("wpp_deaths"), 
    extern = table("extern_deaths"))
  
  # Measles deaths with and without vaccine
  measles_deaths = deaths$extern %>%
    format_d_v_a_name() %>%
    # Pathogen and age group of interest...
    filter(d_v_a_name == measles_id, 
           age        <= age_bound) %>%
    pivot_longer(cols = c(vaccine, no_vaccine), 
                 names_to = "case") %>%
    # Summarise over age bins of interest...
    group_by(case, year) %>%
    summarise(measles = sum(value)) %>%
    ungroup() %>%
    as.data.table()
  
  # All cause deaths from WPP
  all_deaths = deaths$wpp %>%
    # Age group of interest...
    filter(age <= age_bound) %>%
    group_by(year) %>%
    summarise(deaths = sum(deaths)) %>%
    ungroup() %>%
    # Repeat for each vaccine scenario...
    expand_grid(case = names(case_dict)) %>%
    select(case, year, deaths) %>%
    # Append number of deaths averted from vaccination...
    left_join(y  = measles_deaths, 
              by = c("case", "year")) %>%
    group_by(year) %>%
    mutate(averted = measles - measles[case == "vaccine"]) %>%
    ungroup() %>%
    # Increment all cause deaths in scenario of no vaccination...
    mutate(all_cause = deaths + averted) %>%
    select(case, year, all_cause) %>%
    arrange(case, year) %>%
    as.data.table()
  
  # Combine causes of death for measles in context
  context_dt = all_deaths %>%
    left_join(y  = measles_deaths, 
              by = c("case", "year")) %>%
    mutate(other = all_cause - measles) %>% 
    # Melt cause of death to long format...
    pivot_longer(cols = names(cause_dict), 
                 names_to  = "cause", 
                 values_to = "abs") %>%
    # Calculate normalised cause of death...
    group_by(case, year) %>%
    mutate(norm = abs / abs[cause == "all_cause"]) %>%
    ungroup() %>%
    # Calculate cumulative values...
    group_by(case, cause) %>%
    mutate(cum = cumsum(abs)) %>%
    ungroup() %>%
    # Convert to tidy format...
    pivot_longer(cols = any_names(metric_dict), 
                 names_to = "metric") %>%
    # Set appropriate order...
    select(case, metric, cause, year, value) %>%
    arrange(case, metric, cause, year) %>%
    as.data.table()
  
  # ---- Construct plotting datatables ----
  
  # First update metric dictionary with values
  metric_dict = list(
    cum  = sprintf(metric_dict$cum,  age_bound),
    abs  = sprintf(metric_dict$abs,  age_bound),
    norm = sprintf(metric_dict$norm, age_bound),
    cov  = sprintf(metric_dict$cov,  min(o$years), max(o$years)))
  
  # Lines for measles and all cause
  line_dt = context_dt %>%
    filter(cause != "other") %>%
    # Don't plot normalised all cause as always one...
    filter(!(metric == "norm" & 
               cause == "all_cause")) %>%
    # Recode variables and set ordering...
    mutate(case   = recode(case,   !!!case_dict), 
           case   = factor(case,   case_dict), 
           metric = recode(metric, !!!metric_dict), 
           metric = factor(metric, metric_dict), 
           cause  = recode(cause,  !!!cause_dict),
           cause  = factor(cause,  cause_dict))
  
  # Area for measles and other cause in vaccine case only
  area_dt = context_dt %>%
    filter(cause != "all_cause", 
           case  == "vaccine") %>%
    select(-case) %>%
    # Recode variables and set ordering...
    mutate(metric = recode(metric, !!!metric_dict), 
           metric = factor(metric, metric_dict), 
           cause  = recode(cause,  !!!cause_dict),
           cause  = factor(cause,  cause_dict))
  
  # Gloabl measles coverage over time
  coverage_dt = table("coverage_global") %>%
    filter(d_v_a_id %in% coverage_ids) %>%
    # Append metric details...
    left_join(y  = table("vaccine_name"), 
              by = "vaccine") %>%
    mutate(metric = metric_dict$cov, 
           metric = factor(metric, metric_dict)) %>%
    # Tidy up...
    select(metric, vaccine_name, year, value = coverage) %>%
    as.data.table()
  
  # ---- Produce plot ----
  
  # Plot measles in context of all cause deaths
  g = ggplot(area_dt) +
    aes(x = year, 
        y = value) +
    # Plot cause of death...
    geom_area(
      mapping = aes(fill = fct_rev(cause)),
      alpha   = 0.5,
      show.legend = FALSE) +
    # Plot all cause deaths...
    geom_line(
      data    = line_dt, # [case == case_dict[[1]], ]
      mapping = aes(
        colour = cause, 
        linetype = case), 
      linewidth = 1.1) +
    # Set primary colour scheme...
    scale_colour_manual(
      name   = "Cause of death", 
      values = colours$line) + 
    scale_fill_manual(
      values = colours$area) +
    # Prettify primary legends...
    guides(colour = guide_legend(
      order   = 2, 
      reverse = TRUE)) + 
    # Plot number of doses...
    ggnewscale::new_scale_color() + 
    geom_line(
      data    = coverage_dt,
      mapping = aes(colour = vaccine_name)) +
    # Set doses colour scheme...
    scale_colour_manual(
      name   = "Vaccine", 
      values = colours$doses) + 
    # Facet by metric...
    facet_wrap(
      facets = vars(metric),
      scales = "free_y",
      ncol   = 1) +
    # Prettiy x axis...
    scale_x_continuous(
      limits = c(min(o$years), max(o$years)), 
      expand = expansion(mult = c(0, 0)), 
      breaks = seq(min(o$years), max(o$years), by = 5)) +
    # Prettify y axis...
    expand_limits(y = c(0, 1)) + 
    scale_y_continuous(
      labels = comma, 
      expand = expansion(mult = c(0, 0)), 
      breaks = pretty_breaks()) +
    # Prettify primary legends...
    guides(
      colour    = guide_legend(
        order   = 3), 
      linetype  = guide_legend(
        order   = 1, 
        title   = "Vaccine scenario"))
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.title    = element_blank(),
          axis.text     = element_text(size = 10),
          axis.text.x   = element_text(hjust = 1, angle = 50), 
          axis.line     = element_blank(),
          strip.text    = element_text(size = 16),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(1, "lines"),
          legend.title  = element_text(size = 14),
          legend.text   = element_text(size = 12),
          legend.margin = margin(t = 40, b = 40, l = 10),
          legend.position = "right", 
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines"))
  
  # Save figure to file
  save_fig(g, "Measles in context", dir = "historical_impact")
  save_fig(g, "Figure S2", dir = "manuscript")
}

# ---------------------------------------------------------
# Plot comparison of EPI50 outcomes vs VIMC outcomes
# ---------------------------------------------------------
plot_vimc_comparison = function() {
  
  message("  - Plotting comparison of EPI50 vs VIMC outcomes")
  
  # Dictionary for source of predictions
  type_dict = c(
    impact_epi50 = "Deaths averted: EPI50", 
    impact_vimc  = "Deaths averted: VIMC", 
    fvps_epi50   = "Fully vaccinated people: EPI50", 
    fvps_vimc    = "Fully vaccinated people: VIMC")
  
  # Dictionary for metrics of interest
  metric_dict = c(
    impact = "Deaths averted", 
    fvps   = "Fully vaccinated people")
  
  # ---- Load EPI50 and VIMC outcomes ----
  
  # Store all plotting values in list
  plot_list = list(
    
    # EPI50 impact 
    a = read_rds("history", "burden_averted_deaths") %>%
      lazy_dt() %>%
      group_by(d_v_a_id) %>%
      summarise(value = round(sum(impact))) %>%
      ungroup() %>%
      mutate(metric = "impact", 
             type   = "epi50") %>%
      as.data.table(), 
    
    # VIMC impact
    b = table("vimc_predictions") %>%
      lazy_dt() %>%
      group_by(d_v_a_id) %>%
      summarise(value = sum(deaths_averted)) %>%
      ungroup() %>%
      mutate(metric = "impact", 
             type   = "vimc") %>%
      as.data.table(), 
    
    # EPI50 FVPs
    c = table("coverage") %>%
      lazy_dt() %>%
      group_by(d_v_a_id) %>%
      summarise(value = round(sum(fvps))) %>%
      ungroup() %>%
      mutate(metric = "fvps", 
             type   = "epi50") %>%
      as.data.table(),
    
    # VIMC FVPs
    d = table("coverage_source") %>%
      filter(source == "vimc") %>%
      lazy_dt() %>%
      group_by(d_v_a_id) %>%
      summarise(value = round(sum(fvps))) %>%
      ungroup() %>%
      mutate(metric = "fvps", 
             type   = "vimc") %>%
      as.data.table())
  
  # ---- Concatenate all outcomes ---- 

  # Combine all plotting data
  plot_dt = rbindlist(plot_list) %>%
    left_join(y  = table("d_v_a"), 
              by = "d_v_a_id") %>%
    filter(source == "vimc") %>%
    # Append d_v_a names...
    select(d_v_a_name, metric, type, value) %>%
    mutate(type = paste1(metric, type)) %>%
    # Set plotting order...
    mutate(metric = recode(metric, !!!metric_dict), 
           metric = factor(metric, metric_dict), 
           type   = recode(type,   !!!type_dict), 
           type   = factor(type,   type_dict)) %>%
    as.data.table()
  
  # ---- Produce plot ----
  
  # Plot bars of EPI50 outcomes vs VIMC
  g = ggplot(plot_dt) + 
    aes(x = d_v_a_name, 
        y = value, 
        fill = type) + 
    geom_col(position = "dodge") + 
    # Facet by metric...
    facet_wrap(
      facets = vars(metric), 
      scales = "free_y", 
      ncol   = 1) + 
    # Set paired colour scheme...
    scale_fill_manual(
      values = rev(colour_scheme(
        map = "brewer::paired", 
        n   = length(type_dict)))) + 
    # Prettify y axis...
    scale_y_continuous(
      labels = comma)
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.title    = element_blank(),
          axis.text     = element_text(size = 11),
          axis.text.x   = element_text(hjust = 1, angle = 50), 
          axis.line     = element_blank(),
          strip.text    = element_text(size = 14),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(1, "lines"),
          panel.grid.major.y = element_line(linewidth = 0.25),
          legend.title  = element_blank(),
          legend.text   = element_text(size = 11),
          legend.key    = element_blank(),
          legend.position = "right", 
          legend.spacing.y  = unit(1, "lines"),
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines"))
  
  # Save theis figure to file
  save_fig(g, "Results comparison with VIMC", 
           dir = "historical_impact")
}

# **** Other (to be reviewed) ****

# ---------------------------------------------------------
# Helen's exploratory figures
# ---------------------------------------------------------
plot_exploratory = function() {
  
  browser() # To be fully integrated...
  
  # Plot predictions of regression predictors by region
  plot_model_fit_df = model_fit %>%
    select(country, region_short, term, prediction) %>%
    pivot_wider(names_from = term,
                values_from = prediction) %>%
    group_by(region_short) %>%
    as.data.table()
  
  # Plot Spearman rank correlation between coefficients of regression
  plot_model_fit_df %>% 
    select(`log(coverage)`,
           `log(coverage_minus_1)`,
           `log(coverage_minus_2)`,
           `log(coverage_minus_3)`,
           `log(coverage_minus_4)`,
           gini,
           HDI,
           `(Intercept)`,
           attended_births) %>%
    ggpairs()
  
  ggpairs(plot_model_fit_df, columns = c(3,5,6), ggplot2::aes(colour=region_short))
  
  # Plot density of coefficients of predictors by region
  ggpairs(plot_model_fit_df,
          columns = c(3,6:10),
          aes(colour = region_short,
              alpha = 0.5))
  
  # Plot data vs fitted for a single country
  plot_df = augment(model_1) %>%
    filter(country == "THA")
  
  ggplot(data = plot_df, aes(x = target, y = .fitted)) +
    geom_point() +
    labs(
      y = "Fitted (predicted values)",
      x = "Data (actual values)",
      title = paste("Vaccine impact of", d_v_a_name, "in", plot_df$country)
    ) +
    geom_abline(intercept = 0, slope = 1)
  
  # Plot data vs. fitted for all countries
  plot_df = augment(best_model)
  
  ggplot(data = plot_df, aes(x = target, y = .fitted)) +
    geom_point() +
    labs(
      y = "Fitted (predicted values)",
      x = "Data (actual values)",
      title = paste("Vaccine impact of", d_v_a_name)
    ) +
    geom_abline(intercept = 0, slope = 1) +
    facet_wrap(~country, ncol = 21)
  
  # Plot model fit for a single country
  plot_df = augment(best_model) %>%
    filter(country == "AGO")
  
  ggplot(data = plot_df, aes(x = year)) +
    geom_point(aes(y = target, colour = "Data")) +
    geom_line(aes(y = .fitted, colour = "Fitted")) +
    labs(y = NULL,
         title = paste("Vaccine impact of", d_v_a_name, "in", plot_df$country)
    ) +
    scale_colour_manual(values=c(Data="black",Fitted="#D55E00")) +
    guides(colour = guide_legend(title = NULL))
  
  # Plot model fit for all countries
  plot_df = augment(best_model)
  
  ggplot(data = plot_df, aes(x = year)) +
    geom_point(aes(y = target, colour = "Data")) +
    geom_line(aes(y = .fitted, colour = "Fitted")) +
    labs(y = NULL,
         title = paste("Vaccine impact of", d_v_a_name)
    ) +
    scale_colour_manual(values=c(Data="black",Fitted="#D55E00")) +
    guides(colour = guide_legend(title = NULL))  +
    facet_wrap(~country, ncol = 21)
  
  # Manually explore associations between predictor variables for different geographical regions and time points
  explore_dt =  data_dt %>% as.data.table() %>% # Transform to data table to remove country as categorical variable
    filter(#year > 2000 & year <= 2020 &
      region_short == "AFR" &
        target > 2e-20) %>%
    select(-country) %>%
    select(target, gini, health_spending, coverage, 
           coverage_minus_1, coverage_minus_2, coverage_minus_3, sdi)
  
  explore_dt %>% 
    ggpairs(upper = list(continuous = wrap("cor", method = "spearman"))) # Use Spearman rank correlation to account for outliers
  
  # Explore model selection by region
  ggplot(data = model_choice, aes(x = model_number)) +
    geom_histogram() +
    facet_wrap(~region_short)
  
  # Explore prob. density of coefficients of predictors
  plot_df = model_fit %>% filter(term == "pop_0to14")
  
  ggplot(data = plot_df, aes(x=prediction) ) +
    geom_density() +
    facet_wrap(~region_short, ncol=1)
}

# ---------------------------------------------------------
# Helen's tornado plot of predictor coefficients
# ---------------------------------------------------------
plot_tornado = function() {
  
  # TODO: Split by decade, facet by region OR d_v_a_id OR predictor
  
  # Explore density of coefficients of predictors
  plot_dt = impute_1$model_fit %>%
    select(-c(country, d_v_a_id.x, d_v_a_id.y, model_number, .model, AICc)) %>%
    filter(#!term == "HDI" &
      !region_short == "NA" &
        p.value <= 0.05) %>%
    mutate(model = region_short)
  
  g = dwplot(plot_dt) + 
    geom_vline(xintercept = 0, linetype = 2) +
    ggtitle(paste0("Predicting ", d_v_a_name, " vaccine impact"))
  
  g = g + theme_classic()
}

# **** Helper functions ****

# ---------------------------------------------------------
# Append metric (deaths and DALYs) descriptive names
# ---------------------------------------------------------
append_metric_name = function(dt) {
  
  # Append metric names
  name_dt = dt %>%
    left_join(y  = table("metric_dict"), 
              by = "metric") %>%
    select(-metric, -metric_short) %>%
    rename(metric = metric_long) %>%
    mutate(metric = fct_inorder(metric))
  
  return(name_dt)
}

# ---------------------------------------------------------
# Append full region names
# ---------------------------------------------------------
append_region_name = function(dt) {
  
  # Throw error if region already appended - too confusing
  if ("region" %in% names(dt))
    stop("Provide plotting datatable without 'region'")
  
  # Country - region name mapping
  region_dt = table("country") %>%
    left_join(y  = table("region_dict"), 
              by = "region") %>%
    select(region = region_name, country) %>%
    # Convert to factors...
    arrange(region, country) %>%
    mutate(region = fct_inorder(region))
  
  # Append region names
  name_dt = dt %>%
    left_join(y  = region_dt, 
              by = "country")
  
  return(name_dt)
}


# ---------------------------------------------------------
# Full descriptive names for disease, vaccine, or vaccine type
# ---------------------------------------------------------
append_d_v_t_name = function(dt) {
  
  # All columns to update
  d_v_t = intersect(
    x = names(dt), 
    y = qc(disease, vaccine, type))
  
  # Iterate through columns to update
  for (x in d_v_t) {
    
    # Column name of full description
    x_name = paste1(x, "name")
    
    # Plotting order
    x_ord = table(x_name)[[x_name]]
    
    # Append full name description
    dt %<>%
      lazy_dt() %>%
      left_join(y  = table(x_name), 
                by = x) %>%
      select(-all_of(x)) %>%
      # Rename and reorder columns...
      rename(.x := all_of(x_name)) %>%
      mutate(.x = factor(.x, x_ord)) %>%
      rename(!!x := .x) %>%
      select(all_names(dt)) %>%
      as.data.table()
  }
  
  return(dt)
}

# ---------------------------------------------------------
# Set natural order of d_v_a names - appending if necessary
# ---------------------------------------------------------
format_d_v_a_name = function(dt) {
  
  # Natural order of d_v_a_name (not necessarily alphabetical)
  d_v_a_dt = table("d_v_a") %>%
    select(d_v_a_id, d_v_a_name)
  
  # Check if d_v_a_name is already defined
  if (!"d_v_a_name" %in% names(dt)) {
    
    # Append name if it's not already in the data
    dt %<>%
      left_join(y  = d_v_a_dt, 
                by = "d_v_a_id")
  }
  
  # Set order using factors according to d_v_a_dt
  order_dt = dt %>%
    mutate(d_v_a_name = factor(
      x      = d_v_a_name, 
      levels = d_v_a_dt$d_v_a_name)) %>%
    arrange(d_v_a_name)
  
  return(order_dt)
}

# ---------------------------------------------------------
# Get vector of colours for a given variable
# ---------------------------------------------------------
get_palette = function(variable) {
  
  # Input argument must be specified in o$palette
  if (!variable %in% names(o$palette))
    stop("Input '", variable, "' must be one of: ", 
         paste(names(o$palette), collapse = ", "))
  
  # Number of colours to generate based on input argument
  n_colours = list(
    disease = n_unique(table("d_v_a")$disease), 
    region  = n_unique(table("country")$region), 
    income  = n_unique(table("income_status")$income))
  
  # Shorthand for palette and number of colours
  p = o$palette[[variable]]
  n = n_colours[[variable]]
  
  # Generate vector of colours
  colours = colour_scheme(p, n = n)
  
  return(colours)
}

# ---------------------------------------------------------
# Save a ggplot figure to file with default settings
# ---------------------------------------------------------
save_fig = function(g, ..., dir = NULL) {
  
  # Collapse inputs into vector of strings
  fig_name_parts = unlist(list(...))
  
  # Construct file name to concatenate with file path
  save_name = paste(fig_name_parts, collapse = " - ")
  
  # Construct path to save file to
  save_path = o$pth$figures
  if (!is.null(dir)) {
    ext_path  = paste(dir, collapse = file_sep())
    save_path = paste0(save_path, ext_path, file_sep())
  }
  
  # Create directory if it exists
  if (!dir.exists(save_path))
    dir.create(save_path)
  
  # Repeat the saving process for each image format in figure_format
  for (fig_format in o$figure_format) {
    full_path = paste0(save_path, save_name, ".", fig_format)
    
    # Save figure (size specified in options.R)
    ggsave(full_path, 
           plot   = g, 
           device = fig_format, 
           dpi    = o$save_resolution, 
           width  = o$save_width, 
           height = o$save_height, 
           units  = o$save_units)
  }
}

    
 