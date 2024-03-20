###########################################################
# COLOURS
#
# Sample WHO colour schemes.
#
# SOURCE:
# https://apps.who.int/gho/data/design-language/design-system/colors/
#
###########################################################

# ---------------------------------------------------------
# Primary colour creation and indexing function
# ---------------------------------------------------------
colours_who = function(type, n) {
  
  # Call colour map according to type
  colour_fn = paste1("colours_who", type)
  
  # Index sdesired number of colours 
  colours = get(colour_fn)()[seq_len(n)]
  
  return(colours)
}

# ---------------------------------------------------------
# Colour palette: categorical
# ---------------------------------------------------------
colours_who_logo = function() {
  
  # WHO logo colours
  colours = c(
    "#009CDE",  # Vibrant blue
    "#001F58")  # Navy blue
  
  return(colours)
}

# ---------------------------------------------------------
# Colour palette: categorical
# ---------------------------------------------------------
colours_who_category = function() {
  
  # 6 colours and 1 grey
  colours = c(
    "#f4a81d", 
    "#f26829", 
    "#bd53bd", 
    "#6363c0", 
    "#008dc9", 
    "#40bf73",
    "#cccccc")
  
  return(colours)
}

# ---------------------------------------------------------
# Colour palette: categorical
# ---------------------------------------------------------
colours_who_region = function() {
  
  # 1 colour per region
  colours = c(
    AFRO  = "#6363c0",
    EMRO  = "#bd53bd",
    EURO  = "#008dc9",
    PAHO  = "#f26829",
    SEARO = "#40bf73",
    WPRO  = "#f4a81d")
  
  return(colours)
}

