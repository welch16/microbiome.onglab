#' function to add default parameters for filter and trim script
#'
#' @param param_name name of the parameter to look for
#' @param param_frame data.frame with the chosen parameters
#' @return the option in param_frame of the default value for
#'     \code{param_name}
#' @export
get_param_filter_trim <- function(param_name, param_frame) {

  out <- param_frame[[param_name]]

  if (is.null(out)) {

    if (param_name %in% c("trimLeft", "trimRight", "maxN", "minQ")) out <- 0
    if (param_name == "truncLen") out <- c(0, 0)
    if (param_name %in% c("maxLen")) out <- Inf
    if (param_name == "maxEE") out <- c(Inf, Inf)
    if (param_name == "minLen") out <- 20
    if (param_name == "truncQ") out <- 2

  }else{
    if (class(out) == "list") {
      out <- out[[1]]
    }else{
      if (out == "Inf") out <- Inf
    }
  }
  out
}


#' generate the input to filter and trim the fastq files
#' with condor
#'
#' @param trim_queue_file name of the file with the queue,
#'  it has 3 columns sample_name | R1.fastq | R2.fastq
#' @param trim_param_file name of the file with the filter and trim parameters
#' @param trim_outdir directory where the output is going to be saved
#' @param condor_file name of the file with condor instructions
#' @param batch_name string with the name of the batch
#' @param request_cores number of cpus per machine
#' @param request_mem number of GB required as memory
#' @export
condor_filter_trim <- function(
  trim_queue_file = "",
  trim_param_file = "",
  trim_outdir = ".",
  condor_file = "./condor_filter_trim",
  batch_name = "dada2_trim_sequences",
  request_cores = 4,
  request_mem = "4 GB") {

  warning_file(trim_queue_file,
    "Need to define filter and trim queue file")
  warning_file(trim_param_file,
    "Need to define filter and trim parameters file")

  rscript <- system("which Rscript", intern = TRUE)

  file_connection <- file(condor_file)

  writeLines(
    c( "universe         = vanilla",
      stringr::str_c("batch_name       = ", batch_name),
      stringr::str_c("executable       = ", rscript),
      stringr::str_c("args             = $(script_r) --sample_name $(sample)",
        "--fastq1 $(fastq1) --fastq2 $(fastq2) --param_file $(param_file)",
        "--outdir $(outdir) --cores", request_cores, sep = " "),
      stringr::str_c("request_cpus     = ", request_cores),
      stringr::str_c("request_memory   = ", request_mem),
      "on_exit_hold     = (ExitBySignal == True) || (ExitCode != 0)",
      stringr::str_c("periodic_release = (NumJobStarts < 10) &&",
        "((CurrentTime - EnteredCurrentStatus) > 300)", sep = " "),
      stringr::str_c("script_r         = ",
        system.file("scripts/filter_and_trim_pairs.R",
        package = "microbiome.onglab")),
      stringr::str_c("queue_file       = ", trim_queue_file),
      stringr::str_c("param_file       = ", trim_param_file),
      stringr::str_c("outdir           = ", trim_outdir),
      stringr::str_c("output           = $(outdir)/out/dada2_filter_trim_",
        "$(sample_name).$(cluster).$(process).out"),
      stringr::str_c("error            = $(outdir)/err/dada2_filter_trim_",
        "$(sample_name).$(cluster).$(process).err"),
      stringr::str_c("log              = $(outdir)/log/dada2_filter_trim_",
        "$(sample_name).$(cluster).$(process).log"),
      "queue sample, fastq1, fastq2 from $(queue_file)"), file_connection)

  close(file_connection)
}
