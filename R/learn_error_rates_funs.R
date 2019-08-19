

##' function to add default parameters for learn error rates script
##'
##' @param param_name name of the parameter to look for
##' @param param_frame data.frame with the chosen parameters
##' @return the option in param_frame of the default value for
##'     \code{param_name}
##' @export
get_param_error_rates <- function(param_name, param_frame){

  out <- param_frame[[param_name]]

  if(is.null(out)){

    if(param_name %in% c("nbases")) out <- 1e8

  }else{
    if(class(out) == "list"){
      out <- out[[1]]
    }else{
      if(out == "Inf") out <- Inf
    }
  }
  out
}


##' generate the input to learn the error rates with condor
##'
##' @param error_queue_files vector with names of the file with the c(R1,R2) fastq files
##' @param error_prefixes vector with the names of the output error rates prefix for the c(R1,R2) files
##' @param error_param_file name of the file with the learn error rates parameters
##' @param error_outdir directory where the output is going to be saved
##' @param condor_file name of the file with condor instructions
##' @param batch_name string with the name of the batch
##' @param request_cores number of cpus per machine
##' @export
condor_error_rates <- function(
  error_queue_files = c("",""),
  error_prefixes = c("R1","R2"),
  error_param_file = "",
  error_outdir = ".",
  condor_file = "./condor_error_rates",
  batch_name = "dada2_error_rates",
  request_cores = 4)
{
  str_c <- stringr::str_c


  warning_file(error_queue_files[1],"Need to define error rates R1 queue file")
  warning_file(error_queue_files[2],"Need to define error rates R2 queue file")
  warning_file(error_param_file,"Need to define error rates parameters file")

  rscript <- system("which Rscript", intern = TRUE)

  file_connection <- file(condor_file)

  writeLines(
    c(
      "universe         = vanilla",
      str_c("batch_name       = ",batch_name),
      str_c("executable       = ", rscript),
      str_c("args             = $(script_r) --rate_files $(rate_file) --param_file $(param_file) --outprefix $(outprefix) --outdir $(outdir) --cores ",request_cores),
      str_c("request_cpus     = ",request_cores),
      "on_exit_hold     = (ExitBySignal == True) || (ExitCode != 0)",
      "periodic_release = (NumJobStarts < 5) && ((CurrentTime - EnteredCurrentStatus) > 180)",
      str_c("script_r         = ", system.file("scripts/learn_error_rates.R",package = "microbiome.onglab")),
      str_c("param_file       = ", error_param_file),
      str_c("outdir           = ", error_outdir),
      "output           = $(outdir)/out/dada2_error_rates_$(outprefix).$(cluster).$(process).out",
      "error            = $(outdir)/err/dada2_error_rates_$(outprefix).$(cluster).$(process).err",
      "log              = $(outdir)/log/dada2_error_rates_$(outprefix).$(cluster).$(process).log",
      str_c("rate_file        = ", error_queue_files[1]),
      str_c("outprefix        = ", error_prefixes[1]),
      "queue 1",
      str_c("rate_file        = ", error_queue_files[2]),
      str_c("outprefix        = ", error_prefixes[2]),
      "queue 1"),file_connection)

  close(file_connection)

}
