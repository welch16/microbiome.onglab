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
