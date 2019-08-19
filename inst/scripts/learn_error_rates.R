#!/usr/bin/env Rscript

## dada2::learnErrors wrap for condor based pipeline
##
## Assumes that we have paired end files, therefore, we
## learn two error rates matrices, i.e. one for each end,
## and that the parameters for dada2 are stored in a json file.
##
## If any of the parameters is not present in the json file,
## or there is not json file, then the code will use the
## dada2's default parameters.

library(optparse)

opt_list <- list(
  make_option("--rate_files", action = "store_true", type = "character",
              help = "Name of a text file with the location of the
        		filtered and trimmed files used to learn the error rates"),
  make_option("--param_file", action = "store_true", type = "character",
              help = "Location of the parameter file in json format, if not provided
        	the script will do dada2::learnErrors with all the default parameters"),
  make_option("--outprefix", action = "store_true", type = "character",
              default = "learned_error_rates",
              help = "Name of the output file with learned error rates, the full file name will
        	be {outdir}/{outprefix}_learned_error_rates.rds"),
  make_option("--outdir", action = "store_true", type = "character",
              default = tempdir(),
              help = "Location of the output directory"),
  make_option("--cores", action = "store_true", type = "numeric",
              default = 4,
              help = "Number of parallel cpus to use")
)


opt <- parse_args(OptionParser(option_list = opt_list))

options(mc.cores = opt$cores)

stopifnot(file.exists(opt$rate_files))

out_file <- file.path(opt$outdir,"error_rates",paste0(opt$outprefix,"_learned_error_rates.rds"))
plot_file <- file.path(opt$outdir,"error_rates",paste0(opt$outprefix,"_diagnostics_error_rates.pdf"))

library(magrittr)
library(tidyverse)
library(dada2)
library(microbiome.onglab)
library(jsonlite)

if(!file.exists(opt$param_file)){
  params <- data.frame()
}else{
  params <- fromJSON(opt$param_file, flatten = TRUE)
}

dada2_params <- c("nbases")

stopifnot( all( names(params) %in% dada2_params))

message("Learning error rates")
message("Input files in: ", opt$rate_files)
message("Output file: ", out_file)

filtered_fastq <- read_tsv(opt$rate_files, col_names = FALSE)
filtered_fastq %<>% pluck(names(filtered_fastq)[1])

filtered_fastq %<>% parse_filtered_file( file.path(opt$outdir,"filter_fastq"))

if(!file.exists(out_file)){
  error_rates <- learnErrors( filtered_fastq, nbases = get_param_error_rates("nbases", params),
                              multithread = TRUE)


  saveRDS(error_rates, file = out_file)

  error_plot <- plotError(error_rates,nominalQ = TRUE)

  ggsave(
    filename = plot_file,
    plot = error_plot,
    width = 20,
    height = 20,
    units = "cm"
  )

  message("Done! error rates saved at ", out_file)
}
