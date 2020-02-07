##  dada2::filterAndTrim wrap for condor based pipeline
##
##  Assumes that we have paired end files,
##  and that the parameters for dada2 are stored in a json file.
##
## If any of the parameters is not present in the json file,
## or there is not json file, then the code will use the
## dada2's default parameters.

library(optparse)

opt_list <- list(
  make_option("--sample_name", action = "store_true", type = "character",
              default = "dada2_sample",
              help = "Name of this sample"),
  make_option("--fastq1", action = "store_true", type = "character",
              help = "Location of the R1 fastq.gz file"),
  make_option("--fastq2", action = "store_true", type = "character",
              help = "Location of the R2 fastq.gz file"),
  make_option("--param_file", action = "store_true", type = "character",
              help = "Location of the parameter file in json format, if not provided
        	the script will do dada2::filterAndTrim with all the default parameters"),
  make_option("--outdir", action = "store_true", type = "character",
              default = tempdir(),
              help = "Location of the output directory"),
  make_option("--cores", action = "store_true", type = "numeric",
              default = 4,
              help = "Number of parallel cpus to use")
)


opt <- parse_args(OptionParser(option_list = opt_list))

options(mc.cores = opt$cores)

if(!file.exists(opt$fastq1)){
  stop("file ", opt$fastq1, "doesn't appear to exist")
}

if(!file.exists(opt$fastq2)){
  stop("file ", opt$fastq2, "doesn't appear to exist")
}

stopifnot(
  file.exists(opt$fastq1),
  file.exists(opt$fastq2))

dir.create(opt$outdir, showWarnings = FALSE)

library(magrittr)
library(tidyverse)
library(microbiome.onglab)
library(dada2)
library(jsonlite)


if(!file.exists(opt$param_file)){
  params <- data.frame()
}else{
  params <- jsonlite::fromJSON(opt$param_file, flatten = TRUE)
}

dada2_params <- c("truncQ","truncLen","trimLeft","trimRight","maxLen",
                  "minLen","maxN","minQ","maxEE")

stopifnot( all( names(params) %in% dada2_params))

fastq1_filtered <- parse_filtered_file(opt$fastq1, file.path(opt$outdir,"filter_fastq"))
fastq2_filtered <- parse_filtered_file(opt$fastq2, file.path(opt$outdir,"filter_fastq"))

message("filter and trim for files:")
message("fastq1: ", opt$fastq1)
message("fastq2: ", opt$fastq2)

message("output files")
message("filtered fastq1: ", fastq1_filtered)
message("filtered fastq2: ", fastq2_filtered)

out_file <- file.path(opt$outdir,"filter_fastq_summary",paste0(opt$sample_name,"_trim_summary.rds"))



if(!file.exists(out_file)){
  out <- filterAndTrim(
    opt$fastq1, fastq1_filtered,
    opt$fastq2, fastq2_filtered,
    truncQ = get_param_filter_trim("truncQ", params),
    truncLen = get_param_filter_trim("truncLen", params),
    trimLeft = get_param_filter_trim("trimLeft", params),
    trimRight = get_param_filter_trim("trimRight", params),
    maxLen = get_param_filter_trim("maxLen", params),
    minLen = get_param_filter_trim("minLen", params),
    maxN = get_param_filter_trim("maxN", params),
    minQ = get_param_filter_trim("minQ", params),
    maxEE = get_param_filter_trim("maxEE", params),
    rm.phix = TRUE,
    compress = TRUE,
    multithread = TRUE)

  summary_out <- dplyr::as_tibble(out,rownames = "sample") %>%
    mutate(
      fold_change = reads.out / reads.in )

  saveRDS(summary_out, file = out_file)
  message("Done! Summary file saved at ", out_file)
}
