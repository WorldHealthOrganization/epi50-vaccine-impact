###########################################################
# AUXILIARY FUNCTIONS
#
# A series of helpful R functions.
#
# Written by A.J.Shattock
###########################################################


# ---------------------------------------------------------
# Simple wrapper for selecting all names of an object
# ---------------------------------------------------------
all_names = function(x) all_of(names(x))

# ---------------------------------------------------------
# Simple wrapper for selecting any names of given object
# ---------------------------------------------------------
any_names = function(x) any_of(names(x))

# ---------------------------------------------------------
# Set as datatable and rename columns in one line
# ---------------------------------------------------------
as_named_dt = function(x, new_names) {
  
  # Convert to datatable
  dt = as.data.table(x)
  
  # Check new names are correct length
  old_names = names(dt)
  if (length(old_names) != length(new_names))
    stop("Inconsistent number of column names provided")
  
  # Set new column names
  named_dt = setnames(dt, old_names, new_names)
  
  return(named_dt)
}

# ---------------------------------------------------------
# Clear the console
# ---------------------------------------------------------
clc = function() cat("\014")

# ---------------------------------------------------------
# Clear all figures
# ---------------------------------------------------------
clf = function() graphics.off()

# ---------------------------------------------------------
# Create colour scheme
# ---------------------------------------------------------
colour_scheme = function(map, pal = NULL, n = 1, ...) {
  
  # Has colour palette been defined
  if (is.null(pal)) {
    
    # That's ok as long as it's defined within the map argument
    if (!grepl("::", map))
      stop("Palette not defined - Use 'pal = my_pal' or 'map = my_map::my_pal'")
    
    # Seperate out the map and the palette
    pal = stringr::str_remove(map, ".*\\::")
    map = stringr::str_remove(map, "\\::.*")
  }
  
  # Initiate colours variable
  colours = NULL
  
  # Built in colour schemes
  if (map == "base")
    colours = get(pal)(n, ...)	
  
  # A load of colour maps from the pals package
  #
  # See: https://www.rdocumentation.org/packages/pals/versions/1.6
  if (map == "pals")
    colours = get(pal)(n, ...)
  
  # Stylish HCL-based colour maps
  #
  # See: https://colorspace.r-forge.r-project.org/articles/hcl_palettes.html
  if (grepl("_hcl$", map))
    colours = get(map)(palette = pal, n = n, ...)
  
  # Colour Brewer colour schemes
  if (map == "brewer")
    colours = brewer_pal(palette = first_cap(pal), ...)(n)
  
  # Viridis colour schemes
  if (map == "viridis")
    colours = viridis_pal(option = pal, ...)(n)
  
  # Throw an error if colours not yetr defined
  if (is.null(colours))
    stop("Colour map '", map, "' not recognised (supported: base, pals, hcl, brewer, viridis)")
  
  return(colours)
}

# ---------------------------------------------------------
# Load data from package and store in self named list
# ---------------------------------------------------------
data_package = function(..., package = NULL) {
  
  # Initiate list to store loaded data
  data_list = list()
  
  # Unpack variable input arguments - these are files to load
  load_files = unlist(list(...))
  
  # Iterate through files to load
  for (load_file in load_files) {
    
    # Load data from specified package into local environment
    eval_str("data(", load_file, ", ", 
             "package = '", package, "', ", 
             "envir = environment())")
    
    # Append data set to list using same name
    data_list[[load_file]] = get(load_file)
    
    # Remove data from local environment
    rm(list = load_file, envir = environment())
  }
  
  return(data_list)
}

# ---------------------------------------------------------
# Apply lappy to each row of a dataframe or datatable
# ---------------------------------------------------------
dtapply = function(dt, fn, ...) {
  y = lapply(seq_row(dt), function(i) fn(dt[i, ], ...))
  return(y)
}

# ---------------------------------------------------------
# Reset R's most annoying default options
# ---------------------------------------------------------
default_R_options = function() {
  options(dplyr.summarise.inform = FALSE, 
          stringsAsFactors = FALSE, 
          scipen = 999)
}

# ---------------------------------------------------------
# Evaluate a string (in calling function environment) using eval
# ---------------------------------------------------------
eval_str = function(...)
  eval(parse(text = paste0(...)), envir = parent.frame(n = 1))

# ---------------------------------------------------------
# Exponential growth that goes through the origin
# ---------------------------------------------------------
exponential_growth = function(x, a, b) {
  y = a * exp(b * x) - a
  return(y)
}

# ---------------------------------------------------------
# Biphasic exponential function
# ---------------------------------------------------------
exp_biphasic = function(x, peak, p, d1, d2, vmax, alpha, beta) {
  
  # Both exponential 'phases': short and long
  exp1 = exp(-x * log(2)/d1) * p
  exp2 = exp(-x * log(2)/d2) * (1-p)
  
  # Combine the phases
  bi_exp = peak * (exp1 + exp2)
  
  # Bound above and below
  bi_exp_scaled = vmax * (1 - 1 / (1 + (bi_exp / beta)^alpha))
  
  return(bi_exp_scaled)
}

# ---------------------------------------------------------
# Double exponential function
# ---------------------------------------------------------
exp_double = function(x, b, g, h, d) {
  y = (h * b * exp(-d * x)) / ((h - b) * exp(-g * x) + b)
  return(y)
}

# ---------------------------------------------------------
# Platform specific file separator - for readability
# ---------------------------------------------------------
file_sep = function() {
  platform_file_sep = .Platform$file.sep
  return(platform_file_sep)
}

# ---------------------------------------------------------
# Capitalise first letter of a string (or vector of strings)
# ---------------------------------------------------------
first_cap = function(string) {
  string_cap = paste0(toupper(substring(string, 1, 1)), substring(string, 2))
  return(string_cap)
}

# ---------------------------------------------------------
# Format heterogeneous styles of dates
# ---------------------------------------------------------
format_date = function(dates, convert = "ymd") {
  styles = c("dmy", "dmY", "ymd", "Ymd")
  dates = parse_date_time(dates, styles)
  dates = get(convert)(dates)
  return(dates)
}

# ---------------------------------------------------------
# Interpolate time series trends
# ---------------------------------------------------------
interp_ts_trend = function(dt) {
  interp_dt = dt %>%
    model(lm = TSLM(log(value) ~ trend())) %>%
    interpolate(dt)
  
  return(interp_dt)
}

# ---------------------------------------------------------
# Convert list to datatable
# ---------------------------------------------------------
list2dt = function(x, ...) {
  dt = rbindlist(lapply(x, as.data.table), ...)
  return(dt)
}

# ---------------------------------------------------------
# Logistic growth that goes through the origin
# ---------------------------------------------------------
logarithmic_growth = function(x, a, b) {
  y = a / (1 + exp(-b * x)) - a / 2
  return(y)
}

# ---------------------------------------------------------
# Logistic function
# ---------------------------------------------------------
logistic = function(x, slope, mid, lower = 0, upper = 1) {
  y = lower + (upper - lower) / (1 + (x / mid) ^ slope)
  return(y)
}

# ---------------------------------------------------------
# Simple wrapper for number of unique observations
# ---------------------------------------------------------
n_unique = function(x) length(unique(x))

# ---------------------------------------------------------
# Wrapper for lapply that also extracts element name
# ---------------------------------------------------------
napply = function(x, fn, ...) {
  y = lapply(seq_along(x), function(i) fn(x[[i]], name = names(x)[i], ...))
  return(y)
}

# ---------------------------------------------------------
# Normalise a vector of values to between 0 and 1
# ---------------------------------------------------------
normalise_0to1 = function(x, x_min = NULL, x_max = NULL, direction = "forward") {
  
  if (!tolower(direction) %in% c("forward", "backward"))
    stop("Normalisation direction must be either 'forward' or 'backward'")
  
  # Forward normalisation
  if (tolower(direction) == "forward") {
    
    # Take bounds from data unless given
    if (is.null(x_min)) x_min = min(x)
    if (is.null(x_max)) x_max = max(x)
    
    # Normalisation equation
    y = (x - x_min) / (x_max - x_min)
    
    # Append original min and max values, needed to backtransform
    attributes(y)$x_min = x_min
    attributes(y)$x_max = x_max
  }
  
  # Backward normalisation
  if (tolower(direction) == "backward") {
    
    # Take bounds from attitubutes of pre-normalised data unless given
    if (is.null(x_min)) x_min = attributes(x)$x_min
    if (is.null(x_max)) x_max = attributes(x)$x_max
    
    # Rearrange equation to solve for x
    #
    # NOTE: as.vector removes all attributes
    y = as.vector(x * (x_max - x_min) + x_min)
  }
  
  return(y)
}

# ---------------------------------------------------------
# Equivalent of paste, but with an underscore instead of space
# ---------------------------------------------------------
paste1 = function(...) paste(..., sep = "_")

# ---------------------------------------------------------
# Convenience wrapper for readRDS
# ---------------------------------------------------------
read_rds = function(pth, ..., err = TRUE) {
  
  # Special use case: pth is the full .rds file path
  if (grepl(".*\\.rds$", pth)) {
    full_path = pth
    
  } else { # Otherwise standard use case
    
    # Construct path and file name using inputs
    file_path = o$pth[[pth]]
    file_name = paste(unlist(list(...)), collapse = "_")
    
    # Concatenate full .rds file path
    full_path = paste0(file_path, file_name, ".rds")
  }
  
  # Check whether file exists
  exists = file.exists(full_path)
  
  # If file exists, load it
  if (exists)
    x = readRDS(file = full_path)
  
  # If file does not exist
  if (!exists) {
    
    # Construct error / warning message
    err_message = paste("Unable to load file '", full_path, "'")
    
    # Either throw an error or warning depending on err argument
    if (err) stop(err_message)
    if (!err) warning(err_message)
    
    # Return out trivial result
    x = NULL
  }
  
  return(x)
}

# ---------------------------------------------------------
# Load Excel files from URL
# ---------------------------------------------------------
read_url_xls = function(url, sheet = 1) {
  
  # Create temporary file
  xls = tempfile()
  
  # Download from URL to temporary file
  download.file(url, xls, quiet = TRUE, mode = 'wb')
  
  # Read the xls file (xlsx also handled)
  url_dt = readxl::read_excel(
    path  = xls, 
    sheet = sheet) %>%
    as.data.table()
  
  # Delete temporary file
  file.remove(xls)
  
  return(url_dt)
}

# ---------------------------------------------------------
# Inverse of cumsum - use to extract the vector which created a cumsum
# ---------------------------------------------------------
rev_cumsum = function(x) {
  
  # Return input if single value
  if (length(x) == 1)
    return(x)
  
  # Take difference of x with a lag of one
  y = x - c(0, x[1 : (length(x) - 1)])
  
  return(y)
}

# ---------------------------------------------------------
# Wrapper for consistent behaviour of base::sample when length(x) is one
# ---------------------------------------------------------
sample_vec = function(x, ...) 
  x[sample.int(length(x), ...)]

# ---------------------------------------------------------
# Convenience wrapper for saveRDS
# ---------------------------------------------------------
save_rds = function(x, pth, ...) {
  
  # Special use case: pth is the full .rds file path
  if (grepl(".*\\.rds$", pth)) {
    full_path = pth
    
  } else { # Otherwise standard use case
    
    # Construct path and file name using inputs
    file_path = o$pth[[pth]]
    file_name = paste(unlist(list(...)), collapse = "_")
    
    # Concatenate full .rds file path
    full_path = paste0(file_path, file_name, ".rds")
  }
  
  # Save as an RDS
  saveRDS(x, file = full_path)
}

# ---------------------------------------------------------
# Simple wrapper for sequence along dataframe rows
# ---------------------------------------------------------
seq_row = function(x) seq_len(nrow(x))

# ---------------------------------------------------------
# Sigmodial growth going through the origin (inverse of logistic)
# ---------------------------------------------------------
sigmoidal_growth = function(x, slope, mid, max) {
  y = 1 - logistic(x, slope, mid, lower = 1 - max, upper = 1)
  return(y)
}

# ---------------------------------------------------------
# Initiate progress bar with normal-use options
# ---------------------------------------------------------
start_progress_bar = function(n) {
  
  # Initiate progress bar from progress package
  pb = progress_bar$new(
    format     = " [:bar] :percent (remaining: :eta)",
    total      = n,     # Number of tasks to complete
    complete   = "-",   # Completion bar character
    incomplete = " ",   # Incomplete bar character
    current    = ">",   # Current bar character
    clear      = TRUE,  # If TRUE, clears the bar when finish
    width      = 125)   # Width of the progress bar
  
  return(pb)
}

# ---------------------------------------------------------
# Bi-directional setdiff - elements not in both x and y
# ---------------------------------------------------------
symdiff = function(x, y) setdiff(union(x, y), intersect(x, y))

# ---------------------------------------------------------
# Format a number with thousand mark separators
# ---------------------------------------------------------
thou_sep = function(val) {
  format_val = format(val, scientific = FALSE,
                      trim = TRUE, 
                      drop0trailing = TRUE, 
                      big.mark = ",")
  return(format_val)
}

# ---------------------------------------------------------
# Load an file if it exists, throw an error if not
# ---------------------------------------------------------
try_load = function(pth, file, msg = NULL, type = "rds", throw_error = TRUE, sep = FALSE) {
  
  # Initiate trivial output
  file_contents = NULL
  
  # Set default error message
  if (is.null(msg))
    msg = "Cannot load file"
  
  # Switch case for loading function
  loading_fnc = switch(
    tolower(type), 
    
    # Support both RDS and CSV
    "rds" = "readRDS", 
    "csv" = "read.csv",
    
    # Throw an error if anything else requested
    stop("File type '", type, "' not supported")
  )
  
  # Concatenate path and file name
  file_name = paste0(pth, ifelse(sep, file_sep(), ""), file, ".", type)
  
  # If file doesn't exist, throw an error if desired
  if (!file.exists(file_name) && throw_error == TRUE)
    stop(msg, " [missing: ", file_name, "]")
  
  # If file exists, try to load it
  if (file.exists(file_name)) {
    
    # Get the loading function and attempt to load file
    file_contents = tryCatch(
      get(loading_fnc)(file_name),
      
      # Catch the error - we may not want to throw it
      error = function(e) {
        
        # Throw descriptive error if desired
        if (throw_error == TRUE) 
          stop(msg, " [unreadable: ", file_name, "]")
      }
    )
  }
  
  return(file_contents)
}

