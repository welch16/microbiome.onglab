##' @importFrom jsonlite toJSON
##' @importFrom jsonlite fromJSON
NULL


move_file <- function(template,paramfile)
{
  warning("recall: edit ",paramfile, " manually")
  params <- jsonlite::fromJSON(template)
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

