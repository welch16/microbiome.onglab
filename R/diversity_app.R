##' @import shiny
NULL


##' Runs a shiny app to perform a diversity analysis
##' @export
run_diversity_app <- function()
{
  app_directory <- system.file("apps","diversity", package = "microbiome.onglab")
  if(app_directory == ""){
    stop("Could not find app diversity app directory. Try re-installing `microbiome.onglab`.",
         call. = FALSE)
  }

  shiny::runApp(app_directory, display.mode = "normal")

}
