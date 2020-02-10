
#' function to add default parameters for merge_pairs script
#'
#' @param param_name name of the parameter to look for
#' @param param_frame data.frame with the chosen parameters
#' @return the option in param_frame of the default value for
#'     \code{param_name}
#' @export
get_param_merge_pairs <- function(param_name, param_frame) {

  out <- param_frame[[param_name]]

  if (is.null(out)) {

    if (param_name == "minOverlap") out <- 12
    if (param_name == "maxMismatch") out <- 0

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
#' @param merge_queue_file name of the file with the queue, it has 3 columns
#'   sample_name | R1.fastq | R2.fastq
#' @param merge_param_file name of the file with the merge pairs parameters
#' @param merge_outdir directory where the output is going to be saved
#' @param error_rate_files vector with names of the learned error rates
#' @param condor_file name of the file with condor instructions
#' @param batch_name string with the name of the batch
#' @param request_cores number of cpus per machine
#' @param request_mem number of GB required as memory
#' @export
condor_merge_pairs <- function(
    merge_queue_file = "",
    merge_param_file = "",
    merge_outdir = ".",
    error_rate_files = c("R1.rds", "R2.rds"),
    condor_file = "./condor_merge_pairs",
    batch_name = "dada2_merge_pairs",
    request_cores = 4,
    request_mem = "4 GB"
) {

  warning_file(merge_queue_file, "Need to define merge pairs file")
  warning_file(merge_param_file, "Need to define merge pairs parameters file")

  file_connection <- file(condor_file)

  rscript <- system("which Rscript", intern = TRUE)
  writeLines(
    c( "universe         = vanilla",
      stringr::str_c("batch_name       = ", batch_name),
      stringr::str_c("executable       = ", rscript),
      stringr::str_c("args             = $(script_r) --sample_name $(sample)",
        "--fastq1 $(fastq1) --fastq2 $(fastq2) --learned_rates1 $(rates1)",
        "--learned_rates2 $(rates2) --param_file $(param_file)",
        "--outdir $(outdir) --cores", request_cores, sep = " "),
      stringr::str_c("request_cpus     = ", request_cores),
      stringr::str_c("request_memory   = ", request_mem),
      "on_exit_hold     = (ExitBySignal == True) || (ExitCode != 0)",
      stringr::str_c("periodic_release = (NumJobStarts < 5) &&",
        "((CurrentTime - EnteredCurrentStatus) > 180)", sep = " "),
      stringr::str_c("script_r         = ",
        system.file("scripts/merge_pairs.R", package = "microbiome.onglab")),
      stringr::str_c("param_file       = ", merge_param_file),
      stringr::str_c("outdir           = ", merge_outdir),
      stringr::str_c("queue_file       = ", merge_queue_file),
      stringr::str_c("output           = $(outdir)/out/",
        "dada2_merge_pairs_$(cluster).$(process)_$(sample_name).out"),
      stringr::str_c("error            = $(outdir)/err/",
        "dada2_merge_pairs_$(cluster).$(process)_$(sample_name).err"),
      stringr::str_c("log              = $(outdir)/log/",
        "dada2_merge_pairs_$(cluster).$(process)_$(sample_name).log"),
      stringr::str_c("rates1           = ", error_rate_files[1]),
      stringr::str_c("rates2           = ", error_rate_files[2]),
      "queue sample, fastq1, fastq2 from $(queue_file)"), file_connection)

  close(file_connection)

}
