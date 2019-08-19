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

plan(multiprocess,workers = opt$cores)

input_files = read_csv(opt$reads_file,col_names = FALSE) %>%
  set_names(c("name","R1","R2"))

input_files %<>%
  mutate(
    summ = file.path(opt$outdir,"filter_fastq_summary",paste0(name,"_trim_summary.rds")) %>%
      future_map(readRDS)
  ) %>%
  unnest() %>%
  select(-fold_change)

input_files %<>%
  select(-R1,-R2) %>%
  mutate(
    reads.merged_pairs = file.path(opt$outdir,"merged_pairs",
                                   paste0(name,"_merged_pairs.rds")) %>%
      future_map(readRDS) %>%
      future_map_dbl( ~ sum(.$abundance)))

asv_table <- readRDS(file.path(opt$outdir,"ASV_tables",paste0(opt$outprefix,"_sequence_table.rds")))

asv_table %<>%
  {
    rs = rowSums(.)
    tibble(
      name = names(rs),
      reads.asv_table = rs)}

input_files %<>% inner_join(asv_table,by = "name")

input_files %<>%
  mutate(
    perc.out = reads.out / reads.in,
    perc.merged_pairs = reads.merged_pairs / reads.in,
    perc.asv_table = reads.asv_table / reads.in
  )

input_files %>% saveRDS(out_file)
