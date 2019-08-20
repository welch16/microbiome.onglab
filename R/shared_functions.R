##' @import magrittr
##' @import stringr
##' @import dplyr
##' @import readr
##' @import optparse
##' @import dada2
##' @import DECIPHER
##' @import furrr
NULL


##' converts from input files into the filtered file format for dada2 pipeline
##'
##' @param fastq a character vector with the name of a file
##' @param outdir a character with the name of the output directory
##'
##' @return a character vector with modified
##' @export
parse_filtered_file <- function(fastq,outdir){

  file.path(outdir,
            stringr::str_replace(basename(fastq), ".fastq", "_filtered.fastq.gz"))

}


##' write csv files structured to make the dada2 pipeline
##' run with condor
##'
##' @param all_samples a \code{tibble} with 3 columns:
##'    name, R1, R2
##' @param outfile the file where the full tibble is going to be saved,
##'   the R1 and R2 files, are going to share almost the same name
##'   but replacing `.csv` with `_R1.csv` and `_R2.csv`, respectively.
##' @export
save_files <- function(all_samples, outfile )
{
  R1 <- R2 <- NULL

  all_samples %>% write_csv( outfile, col_names = FALSE)
  all_samples %>% select(R1) %>%
    write_csv(str_replace(outfile,".csv","_R1.csv"), col_names = FALSE)
  all_samples %>% select(R2) %>%
    write_csv(str_replace(outfile,".csv","_R2.csv"), col_names = FALSE)

}

warning_file <- function(file,warn_text)
{
  if(!file.exists(file)){
    if(file == ""){
      warning(warn_text)
    }else{
      warning(file," doesn't exists")
    }
  }
}

##' creates the file structure for running dada2 with condor
##' @param outdir Name of the base output directory
##' @export
create_file_structure <- function(outdir)
{
  dir.create(outdir, showWarnings = FALSE)

  ### condor stuff

  dir.create(file.path(outdir,"err"), showWarnings = FALSE)
  dir.create(file.path(outdir,"log"), showWarnings = FALSE)
  dir.create(file.path(outdir,"out"), showWarnings = FALSE)

  dir.create(file.path(outdir,"summary"), showWarnings = FALSE)
  dir.create(file.path(outdir,"filter_fastq"), showWarnings = FALSE)
  dir.create(file.path(outdir,"filter_fastq_summary"), showWarnings = FALSE)
  dir.create(file.path(outdir,"error_rates"), showWarnings = FALSE)
  dir.create(file.path(outdir,"merged_pairs"), showWarnings = FALSE)
  dir.create(file.path(outdir,"merged_pairs_summary"), showWarnings = FALSE)
  dir.create(file.path(outdir,"ASV_tables"), showWarnings = FALSE)
  dir.create(file.path(outdir,"ASV_tables","idtaxa"), showWarnings = FALSE)

}

##' generates all condor files
##' @param condordir directory where all condor files are going to be saved
##' @param prefix prefix to be used in all files
##' @export
condor_generate_all <- function(condordir, prefix)
{
  my_date <- Sys.Date()
  my_date <- stringr::str_replace_all(my_date,"-","_")

  prefix <- file.path(condordir,str_c(prefix,"_dada2_"))

  condor_filter_trim(condor_file = str_c(prefix,"_filter_and_trim_",my_date))
  condor_error_rates(condor_file = str_c(prefix,"_learned_error_rates_",my_date))
  condor_merge_pairs(condor_file = str_c(prefix,"_merge_pairs_",my_date))
  condor_make_sequence_table(condor_file = str_c(prefix,"_make_seqtab_",my_date))
  condor_remove_chimeras(condor_file = str_c(prefix,"_remove_chimeras_",my_date))
  condor_label_taxa(condor_file = str_c(prefix,"_label_taxa_",my_date))

}
