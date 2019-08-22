#!/usr/bin/env Rscript

library(optparse,quietly = TRUE)

optList <-  list(
  make_option("--asv_file", action = "store_true", type = "character",
              help = "A rds file with the ASV table previously computed"),
  make_option("--taxa_model", action = "store_true", type = "character",
              help = "Name of the learned model used to label the ASV sequences"),
  make_option("--outprefix", action = "store_true", type = "character",
              default = "dada2",
              help = "Name of the output file with the labelled ASVs after
                    	removing the bimeras, the full file name will
        				be {outdir}/ASV_tables/idtaxa/{outprefix}_labels.rds"),
  make_option("--outdir", action = "store_true", type = "character",
              default = tempdir(),
              help = "Location of the output directory"),
  make_option("--cores", action = "store_true", type = "numeric",
              default = 4,
              help = "Number of parallel cpus to use"))

opt <- parse_args(OptionParser(option_list = optList))

options(mc.cores = opt$cores)

out_file = file.path(opt$outdir,"ASV_tables","idtaxa",
                     paste0(opt$outprefix,"_labels.rds"))


stopifnot(file.exists(opt$asv_file))

library(magrittr)
library(tidyverse)
library(dada2)
library(microbiome.onglab)
library(DECIPHER)
library(jsonlite)

message("Input files in: ", opt$dada_file)
message("Output file: ", out_file)

if(!file.exists(out_file)){

  asv_table <- readRDS(opt$asv_file)
  dna <- DNAStringSet(getSequences(asv_table)) # Create a DNAStringSet from the ASVs
  names(dna) <- paste0("seq-",seq_along(dna))

  ## label taxa
  taxa_id <- get_tax_ids(dna, opt$taxa_model,opt$cores)

  saveRDS(taxa_id,out_file)

}


