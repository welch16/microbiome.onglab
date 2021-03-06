#!/usr/bin/env Rscript

#' Extract fasta file from an ASV table file (post chimera removal).
#' (We need to provide kraken2 the sequences in fasta format.)
#' Puts the resulting *fna file in {outdir}/fasta/

library(optparse)

info <- Sys.info();
message(stringr::str_c(names(info), " : ", info, "\n"))

opt_list <- list(
  optparse::make_option("--asv_file", action = "store_true", type = "character",
    help = "Name of input ASV sequence table (after chimera filter)."),
  optparse::make_option("--outdir", action = "store_true", type = "character",
    help = "Location of the output directory (top level)"),
  optparse::make_option("--prefix", action = "store_true", type = "character",
    help = "Prefix for output filename (will be {outdir}/fasta/{prefix}.fna")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = opt_list))

### get asv tables after removing chimeras
library(magrittr)
library(tidyverse)
library(seqinr)

base_dr <- opt$outdir
my_prefix <- opt$prefix
asv_file <- opt$asv_file

# check file, directory existence
fasta_dir <- file.path(base_dr, "fasta")
out_file <- file.path(fasta_dir, sprintf("%s.fna", my_prefix))

if (!file.exists(out_file)) {
  message("Output fasta file ", out_file,
    " already exists. Remove if you want to run again.")
} else if (!file.exists(asv_file)) {
  message("Input sequence table ", asv_file, " does not exist.")
} else {
  my_asv <- asv_file %>% readRDS()

  # get sequences from the asv table
  get_sequences <- function(asv) {
    seqs <- colnames(asv)
    seqs <- stringr::str_split(seqs, pattern = "")
    names(seqs) <- stringr::str_c("asv", seq_along(seqs), sep = "_")
    seqs
  }
  sequences <- get_sequences(my_asv)

  # make the fasta dir and write output if we got this far
  dir.create(fasta_dir, showWarnings = FALSE)
  seqinr::write.fasta(
    sequences,
    names(sequences),
    out_file, nbchar = 500)

  message("Done! Wrote fasta file ", out_file)
}