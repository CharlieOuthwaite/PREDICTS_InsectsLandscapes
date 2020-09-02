##%######################################################%##
#                                                          #
#### Function to unscale variables using saved scalers  ####
#                                                          #
##%######################################################%##


unscale <- function(x, scale, centre, logval = FALSE){
  if(logval == FALSE){
    ((x*scale) + centre)
    
  }else{
    #((log(x)*scale) + centre) # this needs correcting
    
  }
  
  
}
