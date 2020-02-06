#!/usr/bin/env Rscript

## dada2::dada and dada2::mergePairs wrap
##
## This script is a bit more elaborated, because instead of wrapping only one
## dada2 method, this wraps the following:
##
## * dada2::derepFastq
## * dada2::dada
## * dada2::mergePairs
##
## The first two methods barely require a parameter, while the third one
## requires a whole set of parameters that tune the behaviour of the matches reads


info=Sys.info();
message(paste0(names(info)," : ",info,"\n"))


library(optparse)

opt_list <- list(
  make_option("--sample_name", action = "store_true", type = "character",
              default = "dada2_sample_",
              help = "Name of this sample"),
  make_option("--fastq1", action = "store_true", type = "character",
              help = "Location of the filtered and trimmed R1 fastq.gz file"),
  make_option("--fastq2", action = "store_true", type = "character",
              help = "Location of the filtered and trimmed R2 fastq.gz file"),
  make_option("--learned_rates1", action = "store_true", type = "character",
              help = "Name of the error rates learned with the fastq R1 files"),
  make_option("--learned_rates2", action = "store_true", type = "character",
              help = "Name of the error rates learned with the fastq R2 files"),
  make_option("--param_file", action = "store_true", type = "character",
              help = "Location of the parameter file in json format, if not provided
        	the script will do dada2::mergePairs with all the default parameters"),
  make_option("--outdir", action = "store_true", type = "character",
              default = tempdir(),
              help = "Location of the output directory"),
  make_option("--cores", action = "store_true", type = "numeric",
              default = 4,
              help = "Number of parallel cpus to use")
)

# [1] "dust_2018_23_GCGTAGTA_TAATCTTA_L001" "190114_C76L2_162_S81_L001"
# [3] "190114_C76L2_151_S69_L001"           "190114_C75PR_238_S44_L001"
# [5] "dust_2018_37_TACGCTGC_GGCTCTGA_L001"
# >


# [1] "190114_C76L2_107_S22_L001"     "190114_C76L2_158_S76_L001"
# [3] "190114_C76L2_237_S117_L001"    "190114_C75PR_193_S10_L001"
# [5] "181221_C7HN5_241_16s_S78_L001"


opt <- parse_args(OptionParser(option_list = opt_list))


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


if(!file.exists(opt$param_file)){
  params <- data.frame()
}else{
  params <- fromJSON(opt$param_file, flatten = TRUE)
}

dada2_params <- c("minOverlap","maxMismatch")

stopifnot( all( names(params) %in% dada2_params))

fastq1 <- parse_filtered_file(opt$fastq1, file.path(opt$outdir,"filter_fastq"))
fastq2 <- parse_filtered_file(opt$fastq2, file.path(opt$outdir,"filter_fastq"))

stopifnot(file.exists(fastq1),file.exists(fastq2))

message("using R1: ", fastq1)
message("using R2: ", fastq2)

### de-replicate the samples (depending on processing time, we may pre-apply this step)
derep_R1 <- derepFastq(fastq1)
derep_R2 <- derepFastq(fastq2)

### load error rates
message("using error rates")
message("* ", opt$learned_rates1)
message("* ", opt$learned_rates2)

error_rates_R1 <- readRDS(opt$learned_rates1)
error_rates_R2 <- readRDS(opt$learned_rates2)

## process both ends, and generate dada objects
message("Processing R1 end with dada2")
dada_R1 <- dada(derep_R1, err = error_rates_R1, multithread = TRUE)
message("Processing R2 end with dada2")
dada_R2 <- dada(derep_R2, err = error_rates_R2, multithread = TRUE)

## merge pairs
message("Merging paired dada2 output")


merged_pairs <- mergePairs(
  dada_R1, derep_R1, dada_R2, derep_R2,
  minOverlap = get_param_merge_pairs("minOverlap",params),
  maxMismatch = get_param_merge_pairs("maxMismatch",params))


out_file <- file.path(opt$outdir,"merged_pairs",paste0(opt$sample_name,"_merged_pairs.rds"))

saveRDS(merged_pairs, file = out_file)

merge_summary <- tibble(
  sample = opt$sample_name,
  denoised_R1 = sum(dada_R1@.Data[[2]]$abundance),
  denoised_R2 = sum(dada_R2@.Data[[2]]$abundance),
  merged_pairs = sum(merged_pairs$abundance))

summary_file <- file.path(opt$outdir,"merged_pairs_summary",paste0(opt$sample_name,"_summary.rds"))
saveRDS(merge_summary,summary_file)

message("Done! merged paired table saved at ",out_file)
message("		and summary file at ",summary_file)
