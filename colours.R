

colours_who = function(type, n) {
  
  colour_fn = paste1("colours_who", type)
  
  colours = get(colour_fn)()[seq_len(n)]
  
  return(colours)
}

colours_who_category = function() {
  
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