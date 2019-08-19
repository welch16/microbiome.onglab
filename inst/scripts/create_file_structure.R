#!/usr/bin/bash

## this script will create the file structure to use this pipeline,
##

library(optparse)

opt_list <- list(
  make_option("--outdir", action = "store_true", type = "character",
              help = "Name of the output directory")
)

opt <- parse_args(OptionParser(option_list = opt_list))

message("creating directory ", opt$outdir)

dir.create(opt$outdir, showWarnings = FALSE)

### condor stuff

dir.create(file.path(opt$outdir,"err"), showWarnings = FALSE)
dir.create(file.path(opt$outdir,"log"), showWarnings = FALSE)
dir.create(file.path(opt$outdir,"out"), showWarnings = FALSE)

dir.create(file.path(opt$outdir,"summary"), showWarnings = FALSE)
dir.create(file.path(opt$outdir,"filter_fastq"), showWarnings = FALSE)
dir.create(file.path(opt$outdir,"filter_fastq_summary"), showWarnings = FALSE)
dir.create(file.path(opt$outdir,"error_rates"), showWarnings = FALSE)
dir.create(file.path(opt$outdir,"merged_pairs"), showWarnings = FALSE)
dir.create(file.path(opt$outdir,"merged_pairs_summary"), showWarnings = FALSE)
dir.create(file.path(opt$outdir,"ASV_tables"), showWarnings = FALSE)
dir.create(file.path(opt$outdir,"ASV_tables","idtaxa"), showWarnings = FALSE)
