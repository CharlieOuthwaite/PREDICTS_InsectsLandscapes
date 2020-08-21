##%######################################################%##
#                                                          #
#### Function to rescale variables using saved scalers  ####
#                                                          #
##%######################################################%##


rescale <- function(x, scale, centre, logval = FALSE){
  if(logval == FALSE){
    (x-centre)/scale
  }else{
    (log(x+1)-centre)/scale
  }
  
  
  
}
