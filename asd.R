############################################################
# ADAPTIVE STOCHASTIC DESCENT
#
# Simple implementation of Adaptive Stochastic Descent; a global
# optimisation algorithm that attempts to quickly locate a high 
# dimensional minimiser for all types of objective. For best results
# use with monte-carlo initialisation (in parallel if necessary).
#
# INPUT:
#   fn* := Ojective function to evaluate and minimise.
#   x0* := Vector of initial paramter points.
# args  := Additional variables to pass into objective function.
#   lb  := Vector of lower bound constraints (length of x0).
#   ub  := Vector of upper bound constraints (length of x0).
#
#  * Denotes mandatory inputs.
#
# OPTIONAL SETTINGS:
#   iters := Maximum number of ASD iterations.
# verbose := Display iteration-by-iteration progress in console.
# 
# OUTPUT:
#     x := Vector of parameter values that minimise fn.
#     y := Minimised value of fn (that is, y = fn(x)).
# y_vec := Vector of y progress over ASD iterations.
#
# EXAMPLE USAGE:
#   asd(fn = obj_function, x0 = runif(n), lb = rep(0, n), ub = rep(1, n))
#
# NOTES:
# - The objective function, fn, must return a list containing a 
#   variable, y, denoting the values of the function at input x. 
#   That is, y = fn(x).
# - Set 'verbose' option to FALSE to turn off console messages.
#
# Credit for the algorthim goes to C.C.Kerr and colleagues:
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0192944
#
# This R implementation written by A.J.Shattock
############################################################

# ---------------------------------------------------------
# Main function: Adaptive Stochastic Descent algorithm.
# ---------------------------------------------------------
asd = function(fn, x0, args = NULL, lb = NULL, ub = NULL, iters = 100) {
  
  # ---- Algorithm settings ----
  
  # Factor for scaling step sizes and selection probabilities
  learn_punish_rate = 2
  
  # Starting step size as a fraction of initial value
  initial_step_fraction = 10
  
  # Maximum selection probability for any parameter
  max_param_prob = 0.8
  
  # ---- Algorithm set up ----
  
  # Number of optimisation parameters
  n_params = length(x0)
  
  # Ensure bounds are appropriate (ub > lb)
  if (!all(ub - lb > 0)) stop("Bounds must be non-trivial")
  
  # Initiate normalised vector of parameter selection probabilities
  param_probs = rep(1, n_params) / n_params
  
  # Initiate change in parameter value for each parameter
  change_val = init_change_val = x0 / initial_step_fraction
  change_dir = (as.numeric(runif(n_params) < 0.5) * 2) - 1  # Yep, pretty ugly! More elegant suggestions welcome :)
  
  # ---- Initiate outcomes ----
  
  # Initial function evaluation
  fn_out = fn(x0, args)
  y0 = y = fn_out$y  # Extract obj fn value
  
  # Vector to store all minimised objective function values
  y_vec = rep(NA, iters + 1)
  y_vec[1] = y0
  
  # Matrix to store history of 'current best' parameter value
  x_mat = matrix(NA, nrow = iters + 1, ncol = n_params)
  x_mat[1, ] = x = x0
  
  # ---- Main algorithm loop ----
  
  # Main optimisation loop
  for (i in 1 : iters) {
    
    # Randomly select a parameter to vary
    param_idx = min(which(runif(1) < cumsum(param_probs)))
    
    # Amount (and direction) to change for chosen parameter
    param_change = change_dir[param_idx] * change_val[param_idx]
    
    # Initiate and and alter a new vector representing x
    x_change = x
    x_change[param_idx] = x[param_idx] + param_change
    
    # Ensure bound constraints are satisfied
    x_change = pmax(pmin(x_change, ub), lb)
    
    # Run the model with this altered vector of parameters
    fn_out   = fn(x_change, args)
    y_change = fn_out$y  # Extract obj fn value
    
    # Check whether we have made any improvement
    if (y_change < y) {
      
      # Increase chance of choosing this parameter 
      param_probs[param_idx] = param_probs[param_idx] * learn_punish_rate
      param_probs = pmin(param_probs, max_param_prob)
      
      # Increase amount to vary this paramter by
      change_val[param_idx]  = change_val[param_idx]  * learn_punish_rate
      
      # NOTE: We are improving, so no need to change direction here
      
      # Update x and y with latest best values
      x = x_change  # Update parameter vector with latest minimizer
      y = y_change  # Update y with latest minimum
      
    } else {  # No improvement either up or down
      
      # Decrease change of choosing this parameter and decrease amount to change
      param_probs[param_idx] = param_probs[param_idx] / learn_punish_rate
      change_val[param_idx]  = change_val[param_idx]  / learn_punish_rate
      
      # If change value has reduced past it's initial value, change direction of search
      if (change_val[param_idx] < init_change_val[param_idx]) {
        change_dir[param_idx] = change_dir[param_idx] * -1
      }
    }
    
    # Store latest minimised objective function value
    y_vec[i + 1]   = y
    x_mat[i + 1, ] = x
    
    # Re-normalise probability vector after it has been altered
    param_probs = param_probs / sum(param_probs)
  }
  
  # ---- Prepare output ----
  
  # Initiate list
  output = list()
  
  # Append optimal parameters, x, and optimised objetive function value, y.	
  output$x = x
  output$y = y
  output$y_vec = y_vec
  
  return(output)
}

