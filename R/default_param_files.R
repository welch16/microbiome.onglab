##' @importFrom jsonlite toJSON
##' @importFrom jsonlite fromJSON
##' @importFrom stringr str_replace_all str_c
NULL


move_file <- function(template,paramfile)
{
  warning("recall: edit ",paramfile, " manually")
  params <- jsonlite::fromJSON(template)
  if(any(is.infinite(unlist(params)))){
    inf_var <- which(is.infinite(unlist(params)))
    params[[names(inf_var)]] = as.character(params[[names(inf_var)]])
  }
  readr::write_lines(jsonlite::toJSON(params), paramfile)

}



##' creates a default param file for the filter and trim process
##' @param paramfile name of the output json file with the filter and trim parameters
##' @export
dada2_param_filter_trim <- function(paramfile)
{
  template <- system.file("template_params/filter_and_trim.json",
                          package = "microbiome.onglab")
  move_file(template,paramfile)
}

##' creates a default param file for the learn error rates parameter
##' @param paramfile name of the output json file with the learn error rates parameters
##' @export
dada2_param_learn_error_rates <- function(paramfile)
{
  template <- system.file("template_params/error_rates.json",
                          package = "microbiome.onglab")
  move_file(template,paramfile)
}

##' creates a default param file for the merge pairs process
##' @param paramfile name of the output json file with the merge pairs parameters
##' @export
dada2_param_merge_pairs <- function(paramfile)
{
  template <- system.file("template_params/merge_pairs.json",
                          package = "microbiome.onglab")
  move_file(template,paramfile)
}


##' creates a default param file for the remove chimera process
##' @param paramfile name of the output json file with the remove chimera parameters
##' @export
dada2_param_remove_chimera <- function(paramfile)
{
  template <- system.file("template_params/remove_chimera.json",
                          package = "microbiome.onglab")
  move_file(template,paramfile)
}


##' copies all the parameter templates to a directory, and
##' uses a prefix to identify all
##' @param outdir directory where all the parameters are saved
##' @param prefix prefix for all the files
##' @export
dada2_param_copy_all <- function(outdir,prefix)
{
  my_date <- stringr::str_replace_all(Sys.Date(),"-","_")
  suffix <- stringr::str_c(my_date,".json")

  dada2_param_filter_trim(file.path(outdir,str_c(prefix,"_filter_and_trim_",suffix)))
  dada2_param_learn_error_rates(file.path(outdir,str_c(prefix,"_learn_error_rates_",suffix)))
  dada2_param_merge_pairs(file.path(outdir,str_c(prefix,"_merge_pairs_",suffix)))
  dada2_param_remove_chimera(file.path(outdir,str_c(prefix,"_remove_chimera_",suffix)))
}

##' creates a default param file for the source tracker check
##' @param paramfile name of the output json file with the sourcetracker parameters
##' @export
sourcetracker_parameters <- function(paramfile)
{
  template <- system.file("templates_params/sourcetracker.json",
                          package = "microbiome.onglab")
  move_file(template, paramfile)
}




