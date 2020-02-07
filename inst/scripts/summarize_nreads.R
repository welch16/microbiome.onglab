#!/usr/bin/env Rscript
# Summarize numbers of reads per step, makes some plots

info=Sys.info();
message(paste0(names(info)," : ",info,"\n"))

library(optparse,quietly = TRUE)

optList <-  list(
  make_option("--queue_file", action = "store_true", type = "character",
              help = "Three-column file with the names of samples and input fastqs"),
  make_option("--prefix", action = "store_true", type = "character",
              default = "dada2",
              help = "Prefix of the output ASV table sans chimeras. Will also be used in output file name
                be {outdir}/summary/{outprefix}_nreads.rds"),
  make_option("--outdir", action = "store_true", type = "character",
              default = tempdir(),
              help = "Top level of output directory (top of file structure for this run)"),
  make_option("--cores", action = "store_true", type = "numeric",
              default = 4,
              help = "Number of parallel cpus to use"))

opt <- parse_args(OptionParser(option_list = optList))

opt=list()
opt$queue_file="run/sample_table.csv"
opt$prefix="stool_vg_2018"
opt$outdir="./run"
opt$cores=4

out_file = file.path(opt$outdir,"summary",paste0(opt$prefix,"_nreads.rds"))

# condor job will hold if we stop
stopifnot(
  file.exists(opt$queue_file),
  dir.exists(opt$outdir))

library(magrittr)
library(tidyverse)
library(furrr)
library(microbiome.onglab)


input_files <- summarize_number_reads(opt$queue_file, opt$prefix,
                                      opt$outdir,opt$cores)
input_files %>% saveRDS(out_file)

summary_absolute <- microbiome.onglab::plot_abundance_per_step(input_files)
summary_relative <- microbiome.onglab::plot_abundance_per_step(input_files, relative = TRUE)


figdir=file.path(opt$outdir,"figs")
if (!dir.exists(figdir)) dir.create(figdir)

abs_fname=file.path( figdir, sprintf("%s_abu_per_step.png", opt$prefix))
rel_fname=file.path( figdir, sprintf("%s_rel_abu_per_step.png", opt$prefix))
ggsave(
  filename = abs_fname,
  plot = summary_absolute,
  width = 6,
  height = 4,
  units = "in")
ggsave(
  filename = rel_fname,
  plot = summary_relative,
  width = 6,
  height = 4,
  units = "in")

message("Saved summary to ", out_file, ", ", abs_fname, ", ", rel_fname)