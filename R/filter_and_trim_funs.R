

##' function to add default parameters
##'
##' @param param_name name of the parameter to look for
##' @param param_frame data.frame with the chosen parameters
##' @return the option in param_frame of the default value for
##'     \code{param_name}
##' @export
get_param_filter_trim <- function(param_name, param_frame){

  out <- param_frame[[param_name]]

  if(is.null(out)){

    if(param_name %in% c("truncLen","trimLeft","trimRight",
                  "maxN","minQ")) out <- 0
    if(param_name %in% c("maxLen", "maxEE")) out <- Inf
    if( param_name == "minLen") out <- 20
    if( param_name == "truncQ") out <- 2

  }else{
    if(class(out) == "list"){
      out <- out[[1]]
    }else{
      if(out == "Inf") out <- Inf
    }
  }
  out
}
