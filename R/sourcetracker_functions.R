##' @importFrom RCurl getURL
NULL


##' function to add default parameters for learn error rates script
##'
##' @param param_name name of the parameter to look for
##' @param param_frame data.frame with the chosen parameters
##' @return the option in param_frame of the default value for
##'     \code{param_name}
##' @export
get_param_sourcetracker <- function(param_name, param_frame){

  out <- param_frame[[param_name]]

  if(is.null(out)){

    if(param_name %in% c("alpha1","alpha2")) out <- 0.001
    if(param_name == "beta") out <- 10

  }else{
    if(class(out) == "list"){
      out <- out[[1]]
    }
  }
  out
}



##' generates the condor file to performing quality control with sourcertracker
##'
##' sourcetracker (Knights et al 2013, https://dx.doi.org/10.1038%2Fnmeth.1650) fits a mixture model of the samples with sources as components, in our
##' case we use the negative controls as such.
##'
##' @param sourcetracker_asvfile name of the rds file with the asv matrix
##' @param sourcetracker_negcontrol name of the rds file with the relationship between samples and negative controls
##' @param sourcetracker_paramfile json file with sourcetracker's parameters
##' @param sourcetracker_outdir directory where the results are saved
##' @param top_asv number of ASVs to be used with sourcetracker to accelerate performance
##' @param max_split_size max. number of samples to do by a cpu
##' @param rarefaction_depth parameter to control for ASVs (OTUs) with higher abundance. This makes all taxa to have at most this abundance
##'   (by default NULL skips this step)
##' @param sourcetracker_prefix prefix of the file with the results
##' @param condor_file file where the condor commands are saved
##' @param batch_name name of the batch
##' @param request_cores Number of cpus to request by condor
##' @param request_mem Memory request by condor
##' @export
condor_sourcetracker <- function(

  sourcetracker_asvfile = "",
  sourcetracker_negcontrol = "",
  sourcetracker_paramfile = "",
  sourcetracker_outdir = ".",
  top_asv = 50, max_split_size = 30, rarefaction_depth = NULL,
  sourcetracker_prefix = "sourcetracker",
  condor_file = "./condor_sourcetracker",
  batch_name = "sourcetracker_run",
  request_cores = 8,
  request_mem = "4 GB")
{
  str_c <- stringr::str_c

  if(!file.exists(sourcetracker_asvfile))warning(sourcetracker_asvfile,"Need to define asv file")
  if(!file.exists(sourcetracker_negcontrol))warning(sourcetracker_negcontrol,"Need to define negative control file")
  if(!file.exists(sourcetracker_paramfile))warning(sourcetracker_paramfile,"Need to define parameters file")

  rscript <- system("which Rscript", intern = TRUE)

  file_connection <- file(condor_file)

  if(is.null(rarefaction_depth)) rarefaction_depth <- "Inf"
  outfile <- file.path(sourcetracker_outdir, str_c(sourcetracker_prefix,".rds"))

  writeLines(
    c(
      "universe         = vanilla",
      str_c("batch_name       = ", batch_name),
      str_c("executable       = ", rscript),
      str_c("args             = $(script_r) --asv_file $(asvfile) --neg_control_file $(negcontrol) --param_file $(param) --top_asvs $(topasv) --max_split_size $(maxs) --rarefaction_depth $(rd) --outfile $(outfile) --cores ", request_cores),
      str_c("request_cpus     = ", request_cores),
      str_c("request_memory   = ", request_mem),
      "on_exit_hold     = (ExitBySignal == True) || (ExitCode != 0)",
      "periodic_release = (NumJobStarts < 5) && ((CurrentTime - EnteredCurrentStatus) > 180)",
      str_c("script_r         = ", system.file("scripts/sourcetracker_run.R", package = "microbiome.onglab")),
      str_c("outdir           = ", sourcetracker_outdir),
      "output           = $(outdir)/out/sourcetracker_$(cluster).$(process).out",
      "error            = $(outdir)/err/sourcetracker_$(cluster).$(process).err",
      "log              = $(outdir)/log/sourcetracker_$(cluster).$(process).log",
      str_c("asvfile         = ", sourcetracker_asvfile),
      str_c("negcontrol      = ", sourcetracker_negcontrol),
      str_c("param           = ", sourcetracker_paramfile),
      str_c("topasv          = ", top_asv),
      str_c("maxs            = ", max_split_size),
      str_c("rd              = ", rarefaction_depth),
      str_c("outfile         = ", outfile),
      "queue 1"), file_connection)

  close(file_connection)

}
