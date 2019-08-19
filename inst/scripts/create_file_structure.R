#!/usr/bin/bash

## this script will create the file structure to use this pipeline,
##

library(optparse)
library(microbiome.onglab)

opt_list <- list(
  make_option("--outdir", action = "store_true", type = "character",
              help = "Name of the output directory")
)

opt <- parse_args(OptionParser(option_list = opt_list))

message("creating directory ", opt$outdir)

create_file_structure(opt$outdir)
