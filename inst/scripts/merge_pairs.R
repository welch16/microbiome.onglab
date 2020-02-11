#!/usr/bin/env Rscript

# dada2::dada and dada2::mergePairs wrap
#
# This script is a bit more elaborated, because instead of wrapping only one
# dada2 method, this wraps the following:
#
# * dada2::derepFastq
# * dada2::dada
# * dada2::mergePairs
#
# The first two methods barely require a parameter, while the third one
# requires a whole set of parameters that tune the behaviour to match the
# ASVs obtained for each pair

info <- Sys.info();
message(stringr::str_c(names(info), " : ", info, "\n"))

library(optparse)

opt_list <- list(
  optparse::make_option("--sample_name", action = "store_true",
    type = "character", default = "dada2_sample_",
    help = "Name of this sample"),
  optparse::make_option("--fastq1", action = "store_true", type = "character",
    help = "Location of the filtered and trimmed R1 fastq.gz file"),
  optparse::make_option("--fastq2", action = "store_true", type = "character",
    help = "Location of the filtered and trimmed R2 fastq.gz file"),
  optparse::make_option("--learned_rates1", action = "store_true",
  type = "character",
  help = "Name of the error rates learned with the fastq R1 files"),
  optparse::make_option("--learned_rates2", action = "store_true",
    type = "character",
    help = "Name of the error rates learned with the fastq R2 files"),
  optparse::make_option("--param_file", action = "store_true",
    type = "character",
    help = "Location of the parameter file in json format, if not provided
      the script will do dada2::mergePairs with all the default parameters"),
  optparse::make_option("--outdir", action = "store_true", type = "character",
    default = tempdir(),
    help = "Location of the output directory"),
  optparse::make_option("--cores", action = "store_true", type = "numeric",
    default = 4,
    help = "Number of parallel cpus to use")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = opt_list))


options(mc.cores = opt$cores)

stopifnot(
  file.exists(opt$fastq1),
  file.exists(opt$fastq2),
  file.exists(opt$learned_rates1),
  file.exists(opt$learned_rates2)
)

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

dada2_params <- c("minOverlap", "maxMismatch")

stopifnot(all(names(params) %in% dada2_params))

fastq1 <- microbiome.onglab::parse_filtered_file(opt$fastq1,
  file.path(opt$outdir, "filter_fastq"))
fastq2 <- microbiome.onglab::parse_filtered_file(opt$fastq2,
  file.path(opt$outdir, "filter_fastq"))

stopifnot(file.exists(fastq1), file.exists(fastq2))

message("using R1: ", fastq1)
message("using R2: ", fastq2)

# de-replicate the samples
# (depending on processing time, we may pre-apply this step)
derep_end1 <- dada2::derepFastq(fastq1)
derep_end2 <- dada2::derepFastq(fastq2)

### load error rates
message("using error rates")
message("* ", opt$learned_rates1)
message("* ", opt$learned_rates2)

error_rates_end1 <- readRDS(opt$learned_rates1)
error_rates_end2 <- readRDS(opt$learned_rates2)

## process both ends, and generate dada objects
message("Processing R1 end with dada2")
dada_end1 <- dada2::dada(derep_end1,
  err = error_rates_end1, multithread = TRUE)
message("Processing R2 end with dada2")
dada_end2 <- dada2::dada(derep_end2,
  err = error_rates_end2, multithread = TRUE)

## merge pairs
message("Merging paired dada2 output")


merged_pairs <- dada2::mergePairs(
  dada_end1, derep_end1, dada_end2, derep_endR2,
  minOverlap = microbiome.onglab::get_param_merge_pairs("minOverlap", params),
  maxMismatch = microbiome.onglab::get_param_merge_pairs("maxMismatch", params))

out_file <- file.path(opt$outdir, "merged_pairs",
  stringr::str_c(opt$sample_name, "_merged_pairs.rds"))

saveRDS(merged_pairs, file = out_file)

merge_summary <- tibble(
  sample = opt$sample_name,
  denoised_R1 = sum(dada_end1@.Data[[2]]$abundance),
  denoised_R2 = sum(dada_end2@.Data[[2]]$abundance),
  merged_pairs = sum(merged_pairs$abundance))

summary_file <- file.path(opt$outdir, "merged_pairs_summary",
  stringr::str_c(opt$sample_name, "_summary.rds"))
saveRDS(merge_summary, summary_file)

message("Done! merged paired table saved at ", out_file)
message("		and summary file at ", summary_file)
