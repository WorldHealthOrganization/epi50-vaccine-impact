###########################################################
# LAUNCH
#
# Main launch function for WHO EPI50 analysis.
#
# Note that this software requires an internet connection.
#
# Authors: A.J.Shattock & H.C.Johnson
###########################################################

# Set working directory to sourced file
if (interactive()) setwd(getSrcDirectory(function() {}))

# Load all required packages and functions
source("dependencies.R")

message("Running EPI50 pipeline")

# Set options (see options.R)
o = set_options(do_step = c(8))

# Step 1) Prepare all inputs (only needs to be done once)
run_prepare()  # See prepare.R

# Step 2) Interface with external polio and measles models
run_external()  # See external.R

# Step 3) Estimate impact for static models using GBD
run_static()  # See static.R

# Step 4) Impute missing countries for VIMC-modelled pathogens
run_regression("impute", "deaths")  # See regression.R
#run_regression("impute", "dalys")   # See regression.R

# Step 5) Fit and select impact-FVP functions
#run_impact("deaths")  # See impact.R  
run_impact("dalys")   # See impact.R

# Step 6) Apply impact functions to historical coverage
#run_history("deaths")  # See history.R
run_history("dalys")   # See history.R

# Step 7) Re-fit time series regression models to infer predictors
# run_regression("infer", "deaths")  # See regression.R

# Step 8) Produce results
run_results()  # See results.R

# Finish up
message("* Finished!")

