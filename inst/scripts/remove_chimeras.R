## merges different ASV tables, and removes chimeras

info <- Sys.info();
message(stringr::str_c(names(info), " : ", info, "\n"))

library(optparse, quietly = TRUE)

opt_list <-  list(
  optparse::make_option("--dada_file", action = "store_true",
    type = "character",
    help = "A file with the name of all the dada2 tables to merge"),
  optparse::make_option("--param_file", action = "store_true",
    type = "character",
    help = "json file with the parameters used to remove bimeras"),
  optparse::make_option("--outprefix", action = "store_true",
    type = "character", default = "dada2",
    help = "Name of the output file with the labelled ASVs after
      removing the bimeras, the full file name will
      be {outdir}/ASV_tables/{outprefix}_asv_wo_bimeras.rds"),
  optparse::make_option("--outdir", action = "store_true", type = "character",
      default = tempdir(),
      help = "Location of the output directory"),
  optparse::make_option("--cores", action = "store_true", type = "numeric",
      default = 4,
      help = "Number of parallel cpus to use"))

opt <- optparse::parse_args(optparse::OptionParser(option_list = opt_list))

options(mc.cores = opt$cores)

out_file <- file.path(opt$outdir, "ASVs",
  stringr::str_c(opt$outprefix, "_asv_wo_bimeras.rds"))

stopifnot(file.exists(opt$dada_file))

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

dada2_params <- c(
  "minSampleFraction",
  "ignoreNNegatives",
  "minFoldParentOverAbundance",
  "allowOneOf",
  "minOneOffParentDistance",
  "maxShift")

stopifnot(all(names(params) %in% dada2_params))

message("Removing bimeras")
message("Input files in: ", opt$dada_file)
message("Output file: ", out_file)

if (!file.exists(out_file)) {

  asv_tables <- readr::read_delim(opt$dada_file, " ", col_names = FALSE) %>%
    dplyr::pull(X1)
  asv_tables %<>% purrr::map(readRDS)

  if (length(asv_tables) > 1) {
    asv_table <- dada2::mergeSequenceTables(tables = asv_tables)
  }else{
    asv_table <- asv_tables[[1]]
  }

  asv_table <- dada2::removeBimeraDenovo(
    asv_table, method = "consensus",
    minSampleFraction =
      microbiome.onglab::get_param_chimeras("minSampleFraction", params),
    ignoreNNegatives =
      microbiome.onglab::get_param_chimeras("ignoreNNegatives", params),
    minFoldParentOverAbundance =
      microbiome.onglab::get_param_chimeras(
        "minFoldParentOverAbundance", params),
    allowOneOf = microbiome.onglab::get_param_chimeras("allowOneOf", params),
    minOneOffParentDistance = microbiome.onglab::get_param_chimeras(
      "minOneOffParentDistance", params),
    maxShift = microbiome.onglab::get_param_chimeras("maxShift", params),
    multithread = TRUE)

  saveRDS(asv_table, out_file)
} else {
  message("Output file ", out_file,
  " already exists. If you want to rerun this script, delete that file first.")
}
