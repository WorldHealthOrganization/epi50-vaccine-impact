###########################################################
# RESULTS
#
# Call plotting functions as defined by o$plot_xxx flags (see
# options.R). All plotting functions themselves live in 
# plotting.R.
#
###########################################################

# ---------------------------------------------------------
# Call plots as defined by o$plot_xxx flags
# ---------------------------------------------------------
run_results = function() {
  
  # Only continue if specified by do_step
  if (!is.element(9, o$do_step)) return()
  
  message("* Producing results")
  
  # ---- Input data plots ----

  # Check plotting flag
  if (o$plot_inputs) {
    
    # Methodology pathogen-country-scope figure
    plot_scope()

    # Total number of FVP over time by source
    plot_total_fvps()

    # Plot coverage density by disease
    plot_coverage()

    # Coverage data density by age
    plot_coverage_age_density()

    # Missing coverage data by country
    plot_missing_data()

    # Global Burden of Disease death estimates by age
    plot_gbd_estimates()
    
    # Proportion of GBD burden we have coverage data for
    plot_gbd_missing()
  }

  # ---- Static model plots ----

  # Check plotting flag
  if (o$plot_static) {

    # Plot vaccine efficacy profiles for static model pathogens
    plot_vaccine_efficacy()

    # Effective coverage with waning immunity for static model pathogens
    plot_effective_coverage()

    # Deaths and DALYs averted for static model pathogens
    plot_static()
  }
  
  # ---- Imputation plots ----
  
  # Check plotting flag
  if (o$plot_imputation) {
    
    # Plot model choice by region
    plot_model_choice()
    
    # Plot predictor and response relationships
    plot_covariates()

    # Plot imputation quality of fit
    plot_impute_quality()
    
    # Plot train-predict countries
    plot_impute_countries()
  }
  
  # ---- Impact function plots ----
  
  # Check plotting flag
  if (o$plot_impact) {
    
    # Exploratory plots of data used to fit impact functions
    plot_impact_data()

    # Plot all-time impact per FVPs
    plot_impact_fvps(scope = "all_time")
    
    # Plot impact vs coverage by vaccine, income, and decade 
    # plot_impact_coverage()
    
    # Plot function selection statistics
    plot_model_selection()
    
    # Plot impact function evaluation
    plot_model_fits()
  }
  
  # ---- Historical results ----
  
  # Check plotting flag
  if (o$plot_history) {
    
    # Inital impact ratios used to back project
    # plot_impact_fvps(scope = "initial")

    # Main results plot - historical impact over time
   #plot_historical_impact()

    # Child mortality rates over time with and without vaccination
    #plot_child_mortality()

    # Regional differences in child mortality changes
    #plot_mortality_change()

    # Measles deaths in context of all cause deaths
    #plot_measles_in_context()

    #
    plot_prob_death_age()

    #
   # plot_survival_increase()
    
    # Plot comparison of EPI50 outcomes vs VIMC outcomes
    #plot_vimc_comparison()
  }
}

