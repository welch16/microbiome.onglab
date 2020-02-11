
# dada2::learnErrors wrap for condor based pipeline
#
# Assumes that we have paired end files, therefore, we learn two error rates
# matrices, i.e. one for each end, and that the parameters for dada2 are stored
# in a json file.
#
# If any of the parameters is not present in the json file, or there is not
# json file, then the code will use the dada2's default parameters.
#
# We changed the file to use the same queue file than the other scripts, and
# included another parameter end to select the column


info <- Sys.info();
message(stringr::str_c(names(info), " : ", info, "\n"))

library(optparse)

opt_list <- list(
  optparse::make_option("--rate_files", action = "store_true",
    type = "character",
    help = "Name of the file with the queue, it has 3 columns:
      sample_name | R1.fastq | R2.fastq
      The end parameter will select between R1 and R2 ends to learn the
      error rates"),
  optparse::make_option("--R1_end", action = "store_true",
    type = "logical", default = TRUE,
    help = "Logical value indicating if the R1 ends are going to be used to
      learn the error rates matrices. In the FALSE case, the program will
      select the R2 ends"),
  optparse::make_option("--param_file", action = "store_true",
    type = "character",
    help = "Location of the parameter file in json format, if not provided
      the script will do dada2::learnErrors with all the default parameters"),
  optparse::make_option("--outprefix", action = "store_true",
    type = "character", default = "learned_error_rates",
    help = "Name of the output file with learned error rates, the full file
      name will be {outdir}/{outprefix}_learned_error_rates.rds"),
  optparse::make_option("--outdir", action = "store_true", type = "character",
    default = tempdir(),
    help = "Location of the output directory"),
  optparse::make_option("--cores", action = "store_true", type = "numeric",
    default = 4,
    help = "Number of parallel cpus to use")
)

opt <- optparse::parse_args(
  optparse::OptionParser(option_list = opt_list))

options(mc.cores = opt$cores)
stopifnot(file.exists(opt$rate_files))

out_file <- file.path(opt$outdir, "error_rates",
  stringr::str_c(opt$outprefix, "learned_error_rates.rds", sep = "_"))
plot_file <- file.path(opt$outdir, "error_rates",
  stringr::str_c(opt$outprefix, "diagnostics_error_rates.pdf", sep = "_"))

library(magrittr)
library(tidyverse)
library(dada2)
library(microbiome.onglab)
library(jsonlite)

if (!file.exists(opt$param_file)) {
  params <- data.frame()
}else{
  params <- jsonlite::fromJSON(opt$param_file, flatten = TRUE)
}

dada2_params <- c("nbases")

stopifnot(all(names(params) %in% dada2_params))
message("Learning error rates")
message("Input files in: ", opt$rate_files)
message("Output file: ", out_file)

filtered_fastq <- readr::read_csv(opt$rate_files, col_names = FALSE)

# selects column based on the R1_end parameter
end_select <- dplyr::if_else(opt$R1_end, 2, 3)
filtered_fastq %<>% purrr::pluck(names(filtered_fastq)[end_select])
filtered_fastq %<>% microbiome.onglab::parse_filtered_file(
  file.path(opt$outdir, "filter_fastq"))
# removes non-existant files from the table
filtered_fastq <- filtered_fastq[file.exists(filtered_fastq)]

if (!file.exists(out_file)) {

  # learn error rates
  error_rates <- dada2::learnErrors(
    filtered_fastq,
    nbases = microbiome.onglab::get_param_error_rates("nbases", params),
    multithread = TRUE)
  saveRDS(error_rates, file = out_file)

  error_plot <- dada2::plotErrors(error_rates, nominalQ = TRUE)
  ggplot2::ggsave(
    filename = plot_file,
    plot = error_plot,
    width = 20,
    height = 20,
    units = "cm"
  )

  message("Done! error rates saved at ", out_file)
}
