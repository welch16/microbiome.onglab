#!/usr/bin/env Rscript

## dada2::makeSequenceTable wrap
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

library(optparse)

opt_list <- list(
  make_option("--queue_file", action = "store_true", type = "character",
              help = "Name of a text file with the location of the
        		filtered and trimmed files used to learn the error rates"),
  make_option("--outprefix", action = "store_true", type = "character",
              default = "learned_error_rates",
              help = "Name of the output file with learned error rates, the full file name will
        	be {outdir}/{outprefix}_sequence_table.rds"),
  make_option("--outdir", action = "store_true", type = "character",
              default = tempdir(),
              help = "Location of the output directory"),
  make_option("--cores", action = "store_true", type = "numeric",
              default = 4,
              help = "Number of parallel cpus to use")
)

opt <- parse_args(OptionParser(option_list = opt_list))

options(mc.cores = opt$cores)

stopifnot(file.exists(opt$queue_file))

library(parallel)
library(magrittr)
library(tidyverse)
library(dada2)

all_files <- read_csv(opt$queue_file, col_names = FALSE)

all_files %<>%
  rename( sample_name = X1) %>%
  mutate(
    out_file = file.path(opt$outdir,"merged_pairs",paste0(sample_name,"_merged_pairs.rds"))
  ) %>%
  select(-X2,-X3)

all_files %<>%
  mutate(
    merged_pairs = mclapply(out_file,readRDS),
    nASVs = map_int(merged_pairs,nrow)
  )

sequence_table <- pull(all_files,merged_pairs)
sequence_table <- makeSequenceTable(sequence_table)
rownames(sequence_table) <- pull(all_files,sample_name)

out_file <- file.path(opt$outdir,"ASV_tables",paste0(opt$outprefix, "_sequence_table.rds"))


saveRDS(sequence_table, out_file)

message("Done! sequence table saved at ", out_file)
