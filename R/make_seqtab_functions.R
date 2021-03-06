
#' generate the input to make the sequence table files with condor
#'
#' @param make_queue_file name of the file with the queue,
#'  it has 3 columns sample_name | R1.fastq | R2.fastq
#' @param make_prefix prefix name for the sequence table output file
#' @param make_outdir directory where the output is going to be saved
#' @param condor_file name of the file with condor instructions
#' @param batch_name string with the name of the batch
#' @param request_cores number of cpus per machine
#' @param request_mem number of GB required as memory
#' @export
condor_make_sequence_table <- function(
  make_queue_file = "",
  make_prefix = "seqtab",
  make_outdir = ".",
  condor_file = "./condor_make_seqtab",
  batch_name = "dada2_make_seqtab",
  request_cores = 4,
  request_mem = "4 GB"
) {

  warning_file(make_queue_file, "Need to define make sequence tab queue file")

  file_connection <- file(condor_file)
  rscript <- system("which Rscript", intern = TRUE)

  writeLines(
    c("universe         = vanilla",
    stringr::str_c("batch_name       = ", batch_name),
    stringr::str_c("executable       = ", rscript),
    stringr::str_c("args             = $(script_r) --queue_file $(queue_file)",
      "--outprefix $(outprefix) --outdir $(outdir) --cores ", request_cores,
      sep = " "),
    stringr::str_c("request_cpus     = ", request_cores),
    stringr::str_c("request_memory   = ", request_mem),
    "on_exit_hold     = (ExitBySignal == True) || (ExitCode != 0)",
    stringr::str_c("periodic_release = (NumJobStarts < 5) &&",
      "((CurrentTime - EnteredCurrentStatus) > 180)", sep = " "),
    stringr::str_c("script_r         = ",
      system.file("scripts/make_sequence_table.R",
      package = "microbiome.onglab")),
    stringr::str_c("outdir           = ", make_outdir),
    stringr::str_c("output           = ",
      "$(outdir)/out/dada2_seqtab_$(outprefix).$(cluster).$(process).out"),
    stringr::str_c("error            = ",
      "$(outdir)/err/dada2_seqtab_$(outprefix).$(cluster).$(process).err"),
    stringr::str_c("log              = ",
      "$(outdir)/log/dada2_seqtab_$(outprefix).$(cluster).$(process).log"),
    stringr::str_c("queue_file       = ", make_queue_file),
    stringr::str_c("outprefix        = ", make_prefix),
    "queue 1"), file_connection)

  close(file_connection)
}

#' function to add default parameters to remove_chimeras script
#'
#' @param param_name name of the parameter to look for
#' @param param_frame data.frame with the chosen parameters
#' @return the option in param_frame of the default value for
#'     \code{param_name}
#' @export
get_param_chimeras <- function(param_name, param_frame) {

  out <- param_frame[[param_name]]

  if (is.null(out)) {

    if (param_name %in% c("minSampleFraction")) out <- 0.9
    if (param_name %in% c("ignoreNNegatives")) out <- 1
    if (param_name %in% c("minFoldParentOverAbundance")) out <- 1.5
    if (param_name %in% c("allowOneOf")) out <- FALSE
    if (param_name %in% c("minOneOffParentDistance")) out <- 4
    if (param_name %in% c("maxShift")) out <- 16

  }else{
    if (class(out) == "list") {
      out <- out[[1]]
    }else{
      if (out == "Inf") out <- Inf
      if (out == "FALSE") out <- FALSE
    }
  }
  out
}

#' generate the input file to remove chimeras from the sequence
#'  table with condor
#'
#' @param sequence_table_file name of the file with the sequence
#'  tables to merge (one p/line)
#' @param chim_param_file name of the remove chimeras parameter file
#' @param chim_prefix prefix name for the sequence table output file
#' @param chim_outdir directory where the output is going to be saved
#' @param condor_file name of the file with condor instructions
#' @param batch_name string with the name of the batch
#' @param request_cores number of cpus per machine
#' @param request_mem number of GB required as memory
#' @export
condor_remove_chimeras <- function(
  sequence_table_file = "",
  chim_param_file = "",
  chim_prefix = "remove_chimeras",
  chim_outdir = ".",
  condor_file = "./condor_remove_chimeras",
  batch_name = "dada2_remove_chimeras",
  request_cores = 4,
  request_mem = "4 GB"
) {

  warning_file(sequence_table_file, "Need to define the remove chimeras file")
  warning_file(chim_param_file,
    "Need to define the parameter file for removing chimeras")

  file_connection <- file(condor_file)
  rscript <- system("which Rscript", intern = TRUE)

  writeLines(
    c(
      "universe         = vanilla",
      stringr::str_c("batch_name       = ", batch_name),
      stringr::str_c("executable       = ", rscript),
      stringr::str_c("args             = $(script_r) --dada_file $(infile)",
        "--param_file $(param) --outprefix $(outprefix)",
        "--outdir $(outdir) --cores ", request_cores, sep = " "),
      stringr::str_c("request_cpus     = ", request_cores),
      stringr::str_c("request_memory   = ", request_mem),
      "on_exit_hold     = (ExitBySignal == True) || (ExitCode != 0)",
      stringr::str_c("periodic_release = (NumJobStarts < 5) &&",
        "((CurrentTime - EnteredCurrentStatus) > 180)", sep = " "),
      stringr::str_c("script_r         = ",
        system.file("scripts/remove_chimeras.R",
        package = "microbiome.onglab")),
      stringr::str_c("infile           = ", sequence_table_file),
      stringr::str_c("param            = ", chim_param_file),
      stringr::str_c("outdir           = ", chim_outdir),
      stringr::str_c("outprefix        = ", chim_prefix),
      stringr::str_c("output           = ",
        "$(outdir)/out/dada2_remove_chimera_$(cluster).$(process)_$(outprefix)",
        ".out"),
      stringr::str_c("error            = ",
        "$(outdir)/err/dada2_remove_chimera_$(cluster).$(process)_$(outprefix)",
        ".err"),
      stringr::str_c("log              = ",
        "$(outdir)/log/dada2_remove_chimera_$(cluster).$(process)_$(outprefix)",
        ".log"),
      "queue 1"), file_connection)

  close(file_connection)

}

