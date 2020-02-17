#' @import magrittr
#' @import stringr
#' @import dplyr
#' @importFrom purrr reduce
#' @import readr
#' @import optparse
#' @import dada2
#' @import DECIPHER
#' @import furrr
#' @import abind
NULL

#' converts from input files into the filtered file format for dada2 pipeline
#'
#' @param fastq a character vector with the name of a file
#' @param outdir a character with the name of the output directory
#'
#' @return a character vector with modified
#' @export
parse_filtered_file <- function(fastq, outdir) {

  file.path(outdir,
    stringr::str_replace(
      basename(fastq), ".fastq", "_filtered.fastq.gz"))

}

#' write csv files structured to make the dada2 pipeline
#' run with condor
#'
#' @param all_samples a \code{tibble} with 3 columns:
#'    name, end1, end2
#' @param outfile the file where the full tibble is going to be saved
#' @export
save_files <- function(all_samples, outfile) {

  all_samples %>% readr::write_csv(outfile, col_names = FALSE)

}

#' creates a warning if the file doesn't exists
#' @param file a string with the name of the file
#' @param warn_text a string with the text to be added as warning
#' @export
warning_file <- function(file, warn_text) {
  if (!file.exists(file)) {
    if (file == "") {
      warning(warn_text)
    }else{
      warning(file, " doesn't exists")
    }
  }
}

#' creates the file structure for running dada2 with condor
#' @param outdir a string with the name of the base output directoryn
#' @export
create_file_structure_group <- function(outdir) {
  dir.create(outdir, showWarnings = FALSE)

  ### condor stuff
  dir.create(file.path(outdir, "err"), showWarnings = FALSE)
  dir.create(file.path(outdir, "log"), showWarnings = FALSE)
  dir.create(file.path(outdir, "out"), showWarnings = FALSE)

  dir.create(file.path(outdir, "summary"), showWarnings = FALSE)
  dir.create(file.path(outdir, "filter_fastq"), showWarnings = FALSE)
  dir.create(file.path(outdir, "filter_fastq_summary"), showWarnings = FALSE)
  dir.create(file.path(outdir, "error_rates"), showWarnings = FALSE)
  dir.create(file.path(outdir, "merged_pairs"), showWarnings = FALSE)
  dir.create(file.path(outdir, "merged_pairs_summary"), showWarnings = FALSE)
  dir.create(file.path(outdir, "ASV_tables"), showWarnings = FALSE)
}


#' create the full file structure that integrated across different groups
#' @param outdir a string with the name of the base output directoryn
#' @param groups a string vector with the names of samples groups
#' @export
#' @examples
#' \dontrun{
#'  groups <- c("group1", "group2")
#'  create_file_structure(outdir, groups)
#' }
create_file_structure <- function(outdir, groups) {

  dir.create(outdir, showWarnings = FALSE)

  # condor stuff

  dir.create(file.path(outdir, "ASVs"))
  dir.create(file.path(outdir, "prevalence"))
  dir.create(file.path(outdir, "kraken"))

  purrr::map(file.path(outdir, groups), create_file_structure_group)

}



#' generates all condor files
#' @param condordir directory where all condor files are going to be saved
#' @param prefix prefix to be used in all files
#' @export
condor_generate_all <- function(condordir, prefix) {

  my_date <- Sys.Date()
  my_date <- stringr::str_replace_all(my_date, "-", "_")
  aux_prefix <- prefix
  prefix <- file.path(condordir, stringr::str_c(prefix, "_dada2"))

  filter_trim_file <- stringr::str_c(prefix, "_filter_and_trim_", my_date)
  error_rates_file <- stringr::str_c(prefix, "_learned_error_rates_", my_date)
  merge_pairs_file <- stringr::str_c(prefix, "_merge_pairs_", my_date)
  seque_table_file <- stringr::str_c(prefix, "_make_seqtab_", my_date)

  condor_filter_trim(
    condor_file = filter_trim_file)
  condor_error_rates(
    condor_file = error_rates_file)
  condor_merge_pairs(
    condor_file = merge_pairs_file)
  condor_make_sequence_table(
    condor_file = seque_table_file)

   condor_generate_group_dag(condordir, basename(aux_prefix),
    filter_trim_file, error_rates_file,
    merge_pairs_file, seque_table_file)

  # condor_remove_chimeras(
  # condor_file = stringr::str_c(prefix, "_remove_chimeras_", my_date))
  # condor_label_taxa(
  #   condor_file = stringr::str_c(prefix, "_label_taxa_", my_date))

}


#' generate the dagman file for running the pipeline with condor
#' @param condordir directory where all condor files are going to be saved
#' @param prefix prefix to be used in all files
#' @param filter_trim_file name of the filter and trim condor file
#' @param error_rates_file name of the learning error rates condor file
#' @param merge_pairs_file name of the merge pairs condor file
#' @param seque_table_file name of the make sequence table condor file
#' @export
condor_generate_group_dag <- function(
  condordir, prefix,
  filter_trim_file,
  error_rates_file,
  merge_pairs_file,
  seque_table_file) {

  my_date <- Sys.Date()
  my_date <- stringr::str_replace_all(my_date, "-", "_")

  condor_file <- file.path(condordir, stringr::str_c(prefix, my_date,
    "dada2.dag", sep = "_"))
  file_connection <- file(condor_file)

  condor_node <- function(parent, child) {

    parent <- purrr::reduce(parent, stringr::str_c, sep = " ")
    child <- purrr::reduce(child, stringr::str_c, sep = " ")

    stringr::str_c("PARENT", parent, "CHILD", child, sep = " ")
  }

  writeLines(
    c(
      stringr::str_c("JOB", "filter_trim", filter_trim_file, sep = " "),
      stringr::str_c("JOB", "error_rates", error_rates_file, sep = " "),
      stringr::str_c("JOB", "merge_pairs", merge_pairs_file, sep = " "),
      stringr::str_c("JOB", "seq_tab", seque_table_file, sep = " "),
      condor_node("filter_trim", "error_rates"),
      condor_node("error_rates", "merge_pairs"),
      condor_node("merge_pairs", "seq_tab")), file_connection)

  close(file_connection)
}