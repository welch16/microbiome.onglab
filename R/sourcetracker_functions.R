#' @importFrom RCurl getURL
#' @importFrom abind abind
NULL

#' loads sourcetracker code from the raw github repository
#' @param env the environment where sourcetracker is going to be loaded. By
#'   default, the sourcetracker functions are loaded to the global environment
#' @export
load_sourcetracker <- function(env = globalenv()) {

  sourcetracker_url <- stringr::str_c("https://raw.githubusercontent.com/",
    "danknights/sourcetracker/master/src/SourceTracker.r")
  script <- RCurl::getURL(sourcetracker_url)
  eval(parse(text = script), envir = env)
}

#' gets the ASV matrix for some samples
#' @param asv the full asv matrix
#' @param samples a vector of samples to select
#' @export
aux_get_matrix <- function(asv, samples) {
  mat <- asv[samples, ]
  if (length(samples) == 1) {
    mat <- matrix(mat, nrow = 1)
  }
  mat
}

#' splits the ASV matrix according to the negative control structure
#' @param group negative control structure indicating the samples for each
#'  group
#' @param asv ASV matrix, nsamples x nASVs
#' @param split_size maximum split size when partitioning the matrix
#' @export
aux_split_run <- function(group, asv, split_size) {

  samples <- dplyr::pull(group, sample)
  sample_table <- aux_get_matrix(asv, samples)

  if (nrow(sample_table) > split_size) {

    splits <- seq_len(ceiling(nrow(sample_table) / split_size))
    splits <- rep(splits, each = split_size)
    splits <- splits[seq_len(nrow(sample_table))]
    split_samples <- split(samples, splits)
    split_samples <- purrr::map(split_samples, ~ sample_table[., ])
  }else{
    split_samples <- list(sample_table)
  }

  split_samples
}

#' merge sourcetracker results into a big sourcetracker object
#' @param split_results a list of sourcetracker objects
#' @return a sourcetracker object
#' @export
merge_sourcetracker_results <- function(split_results) {

  draws <- purrr::map(split_results, "draws")
  draws <- do.call(abind::abind, draws, 3)

  proportions <- purrr::map(split_results, "proportions")
  proportions <- do.call(rbind, proportions)

  proportions_sd <- purrr::map(split_results, "proportions_sd")
  proportions_sd <- do.call(rbind, proportions_sd)

  train_envs <- purrr::map(split_results, "train.envs")
  train_envs <- train_envs[[1]]

  samplenames <- purrr::map(split_results, "samplenames")
  samplenames <- do.call(c, samplenames)

  out <- list(
    "draws" = draws,
    "proportions" = proportions,
    "proportions_sd" = proportions_sd,
    "train.envs" = train_envs,
    "samplenames" = samplenames)

  class(out) <- "sourcetracker.fit"

  invisible(out)
}

#' function to add default parameters for learn error rates script
#'
#' @param param_name name of the parameter to look for
#' @param param_frame data.frame with the chosen parameters
#' @return the option in param_frame of the default value for
#'     \code{param_name}
#' @export
get_param_sourcetracker <- function(param_name, param_frame) {

  out <- param_frame[[param_name]]

  if (is.null(out)) {

    if (param_name %in% c("alpha1", "alpha2")) out <- 0.001
    if (param_name == "beta") out <- 10

  }else{
    if (class(out) == "list") {
      out <- out[[1]]
    }
  }
  out
}



#' generates the condor file to performing quality control with sourcertracker
#'
#' sourcetracker (https://dx.doi.org/10.1038%2Fnmeth.1650)
#' fits a mixture model of the samples with sources as components, in our
#' case we use the negative controls as such.
#'
#' @param sourcetracker_asvfile name of the rds file with the asv matrix
#' @param sourcetracker_negcontrol name of the rds file with the relationship
#'  between samples and negative controls
#' @param sourcetracker_paramfile json file with sourcetracker's parameters
#' @param sourcetracker_outdir directory where the results are saved
#' @param top_asv number of ASVs to be used with sourcetracker to
#'  accelerate performance
#' @param max_split_size max. number of samples to do by a cpu
#' @param rarefaction_depth parameter to control for ASVs (OTUs) with higher
#'  abundance. This makes all taxa to have at most this abundance
#'   (by default NULL skips this step)
#' @param sourcetracker_prefix prefix of the file with the results
#' @param condor_file file where the condor commands are saved
#' @param batch_name name of the batch
#' @param request_cores Number of cpus to request by condor
#' @param request_mem Memory request by condor
#' @export
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
  request_mem = "4 GB") {

  if (!file.exists(sourcetracker_asvfile))
    warning(sourcetracker_asvfile, "Need to define asv file")
  if (!file.exists(sourcetracker_negcontrol))
    warning(sourcetracker_negcontrol, "Need to define negative control file")
  if (!file.exists(sourcetracker_paramfile))
    warning(sourcetracker_paramfile, "Need to define parameters file")

  rscript <- system("which Rscript", intern = TRUE)

  file_connection <- file(condor_file)

  if (is.null(rarefaction_depth)) rarefaction_depth <- "Inf"
  outfile <- file.path(sourcetracker_outdir,
    stringr::str_c(sourcetracker_prefix, ".rds"))

  writeLines(
    c(
      "universe         = vanilla",
      stringr::str_c("batch_name       = ", batch_name),
      stringr::str_c("executable       = ", rscript),
      stringr::str_c("args             = $(script_r) --asv_file $(asvfile)",
        "--neg_control_file $(negcontrol) --param_file $(param)",
        "--top_asvs $(topasv) --max_split_size $(maxs)",
        "--rarefaction_depth $(rd) --outfile $(outfile) --cores",
        request_cores, sep = " "),
      stringr::str_c("request_cpus     = ", request_cores),
      stringr::str_c("request_memory   = ", request_mem),
      "on_exit_hold     = (ExitBySignal == True) || (ExitCode != 0)",
      stringr::str_c("periodic_release = (NumJobStarts < 5) && ",
        "((CurrentTime - EnteredCurrentStatus) > 180)"),
      stringr::str_c("script_r         = ",
        system.file("scripts/sourcetracker_run.R",
        package = "microbiome.onglab")),
      stringr::str_c("outdir           = ", sourcetracker_outdir),
      stringr::str_c("output           = $(outdir)/out/",
        "sourcetracker_$(cluster).$(process).out"),
      stringr::str_c("error            = $(outdir)/err/",
        "sourcetracker_$(cluster).$(process).err"),
      stringr::str_c("log              = $(outdir)/log/",
        "sourcetracker_$(cluster).$(process).log"),
      stringr::str_c("asvfile         = ", sourcetracker_asvfile),
      stringr::str_c("negcontrol      = ", sourcetracker_negcontrol),
      stringr::str_c("param           = ", sourcetracker_paramfile),
      stringr::str_c("topasv          = ", top_asv),
      stringr::str_c("maxs            = ", max_split_size),
      stringr::str_c("rd              = ", rarefaction_depth),
      stringr::str_c("outfile         = ", outfile),
      "queue 1"), file_connection)

  close(file_connection)

}
