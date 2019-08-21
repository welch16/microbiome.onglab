#!/usr/bin/env Rscript

library(optparse,quietly = TRUE)

optList <-  list(
  make_option("--reads_file", action = "store_true", type = "character",
              help = "File with the name of the sequence files"),
  make_option("--outprefix", action = "store_true", type = "character",
              default = "dada2",
              help = "Name of the output file with the labelled ASVs after
                    	removing the bimeras, the full file name will
        				be {outdir}/summary/{outprefix}_nreads.rds"),
  make_option("--outdir", action = "store_true", type = "character",
              default = tempdir(),
              help = "Location of the output directory"),
  make_option("--cores", action = "store_true", type = "numeric",
              default = 4,
              help = "Number of parallel cpus to use"))


opt <- parse_args(OptionParser(option_list = optList))

out_file = file.path(opt$outdir,"summary",paste0(opt$outprefix,"_nreads.rds"))

stopifnot(
  file.exists(opt$reads_file),
  dir.exists(opt$outdir))


library(magrittr)
library(tidyverse)
library(furrr)
library(microbiome.onglab)

input_files <- summarize_number_reads(opt$reads_file, optoutprefix,
                                      opt$outdir,opt$cores)

input_files %>% saveRDS(out_file)
