# dada2::makeSequenceTable wrap
#
# This script is a bit more elaborated, because instead of wrapping only one
# dada2 method, this wraps the following:
#
# * dada2::derepFastq
# * dada2::dada
# * dada2::mergePairs
#
# The first two methods barely require a parameter, while the third one
# requires a whole set of parameters that tune the behaviour of the
# matched reads

info <- Sys.info();
message(stringr::str_c(names(info), " : ", info, "\n"))

library(optparse)

opt_list <- list(
  optparse::make_option("--queue_file", action = "store_true",
    type = "character",
    help = "Name of the master queue file. We will use the sample names to
      locate {outdir}/merged_pairs/{sample_name}_merged_pairs.rds"),
  optparse::make_option("--outprefix", action = "store_true",
    type = "character", default = "learned_error_rates",
    help = "Name of the output file with learned error rates, the full file 
      name will be {outdir}/{outprefix}_sequence_table.rds"),
  optparse::make_option("--outdir", action = "store_true", type = "character",
    default = tempdir(),
    help = "Location of the output directory"),
  optparse::make_option("--cores", action = "store_true", type = "numeric",
    default = 4,
    help = "Number of parallel cpus to use")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = opt_list))

options(mc.cores = opt$cores)

stopifnot(file.exists(opt$queue_file))

library(parallel)
library(magrittr)
library(tidyverse)
library(dada2)
library(microbiome.onglab)

all_files <- readr::read_csv(opt$queue_file, col_names = FALSE)

all_files %<>%
  dplyr::rename(
    sample_name = X1) %>%
  dplyr::mutate(
    out_file = file.path(opt$outdir,
      "merged_pairs", stringr::str_c(sample_name, "_merged_pairs.rds"))
  ) %>%
  dplyr::select(-X2, -X3)

# We will remove rows where file does not exist
all_files %<>%
  dplyr::mutate(has_pairs = file.exists(out_file))

# report which samples we're dropping
all_files %>%
    dplyr::filter(!has_pairs) %>%
    dplyr::select(sample_name, out_file) %>%
    tibble::deframe() %>%
    purrr::walk2(names(.), .f = ~ message("Dropping sample ", .y,
      ". No merged pairs file ", .x))

all_files %<>%
  dplyr::filter(has_pairs) %>%
  dplyr::mutate(
    merged_pairs = mclapply(out_file, readRDS),
    n_ASVs = map_int(merged_pairs, nrow)
  )

sequence_table <- dplyr::pull(all_files, merged_pairs)
sequence_table <- dada2::makeSequenceTable(sequence_table)
rownames(sequence_table) <- dplyr::pull(all_files, sample_name)

out_file <- file.path(opt$outdir, "ASV_tables",
  stringr::str_c(opt$outprefix, "_sequence_table.rds"))

saveRDS(sequence_table, out_file)
message("Done! sequence table saved at ", out_file)
