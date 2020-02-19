#!/usr/bin/env Rscript

info <- Sys.info();
message(stringr::str_c(names(info), " : ", info, "\n"))

library(optparse, quietly = TRUE)

opt_list <-  list(
  optparse::make_option("--asv_file", action = "store_true", type = "character",
    help = "A rds file with the ASV table previously computed"),
  optparse::make_option("--outdir", action = "store_true", type = "character",
    default = tempdir(),
    help = "Location of the output directory"),
  optparse::make_option("--prefix", action = "store_true",
    type = "character", default = "dada2",
    help = "Name of the output file with the labelled ASVs after removing the
       bimeras, the full file name will be
      {outdir}/kraken/{outprefix}_labels.rds"),
  optparse::make_option("--krakendir", action = "store_true",
    default = "/z/Comp/onglab/programs/kraken",
    help = "Directory where kraken is installed, by default we use the full
      path of '/z/Comp/onglab/programs/kraken"),
  optparse::make_option("--krakenmodel", action = "store_true",
    type = "character",
    help = "Name of the file with the learned model used by kraken to
      classify the ASV sequences"),
  optparse::make_option("--krakenprefix", action = "store_true",
    type = "character", default = "kraken",
    help = "Prefix used to name the kraken results according to the db"),
  optparse::make_option("--cores", action = "store_true", type = "numeric",
    default = 4,
    help = "Number of parallel cpus to use"))

opt <- optparse::parse_args(optparse::OptionParser(option_list = opt_list))
options(mc.cores = opt$cores)

stopifnot(file.exists(opt$asv_file))

library(magrittr)
library(tidyverse)
library(taxizedb)
library(furrr)
library(future)

future::plan(multiprocess, workers = opt$cores)

message("Input files in: ", opt$asv_file)


kraken_base <- file.path(opt$outdir, "kraken",
  stringr::str_c(opt$prefix, opt$kraken_prefix, sep = "_"))
kraken_out <- stringr::str_c(kraken_base, ".out")
kraken_summary <- stringr::str_c(kraken_base, ".summary")
kraken_classified <- stringr::str_c(kraken_base, ".classified")
kraken_unclassified <- stringr::str_c(kraken_base, ".unclassified")
out_file <- stringr::str_c(kraken_base, "taxa.rds", sep = "_")

message("Output file: ", out_file)

## generate fasta file
fasta_script <- system.file("scripts", "extract_fasta_before_kraken.R",
  package = "microbiome.onglab")

fasta_file <- file.path(opt$outdir, "fasta",
  stringr::str_c(opt$outprefix, ".fna"))

message("generating fasta file from ASV matrix")
system(
  stringr::str_c(
    fasta_script,
    "--asv_file", opt$asv_file,
    "--outdir", opt$outdir,
    "--prefix", opt$prefix, sep = " ")
)


asv_matrix <- readRDS(opt$asv_file)
asvs <- 
  tibble::tibble(asv = stringr::str_c("asv", seq_len(ncol(asv_matrix)),
    sep = "_"))
ntaxa <- nrow(asvs)


ncbi_db <- taxizedb::db_download_ncbi()

taxa <- c("phylum", "class", "order", "family", "genus", "species")

kraken_output <- readr::read_tsv(kraken_out)
names(kraken_output) <- c("rank", "asv", "id", "seq_length", "id_bp")

kraken_output %<>%
  parse_kraken_taxa_from_id()


# kraken_results %<>%
#   dplyr::mutate(
#     labels = purrr::map( labels, set_names,
#       c("rank", "asv", "id", "seq_length", "id_bp")),
#     labels = furrr::future_map( labels, get_taxa_from_id))
   
# kraken_results %<>%
#   dplyr::mutate(
#     unlabeled = purrr::map( labels,
#       ~ dplyr::mutate(.,
#         labelled = purrr::map_lgl(id_taxa, is.data.frame))) %>%
#           purrr::map( ~ dplyr::filter(., !labelled)) %>%
#           purrr::map( dplyr::select, -id_taxa, -labelled))

# merge_results <- function(label, unlabel) {

#   is_na <- NULL

#   `%<>%` <- magrittr::`%<>%`
#   `%>%` <- dplyr::`%>%`

# 	label %<>%
#     dplyr::filter( ! is_na) %>%
#     dplyr::select(-id_bp) %>%
#     dplyr::mutate(
#       id_taxa = purrr::map(id_taxa, as_tibble))

#   label %<>% dplyr::select(-is_na)

#   label %>%
#     dplyr::mutate( original = TRUE) %>%
#    	dplyr::arrange( as.numeric(stringr::str_remove(asv,"asv_")))
# }

# kraken_results %<>%
# 	dplyr::mutate(
#     out = purrr::map2(labels, unlabeled, merge_results))
#     # ,
#     # outfile = file.path(basedr, "dust_comp", "kraken",
#     # 	stringr::str_c( "2020_02_17_dustcomp_",db, "_kraken_labels.rds")))

# # kraken_results %>%
# # 	dplyr::mutate(
# #     out = purrr::walk2( out, outfile, saveRDS))
  


# clean_id_taxa <- function(id_taxa, taxa) {

#   `%<>%` <- magrittr::`%<>%`
#   name <- NULL

#   if (is.data.frame(id_taxa)) {
#     id_taxa %<>%
#       dplyr::filter(rank %in% taxa) %>%
#       dplyr::select(-id)
#     id_taxa %<>% tidyr::spread(rank, name)
# 	}
#   id_taxa
# }


# clean_id_taxa_wrap <- function(labels, taxa, asvs) {

#   labels %>%
#     dplyr::mutate(
#       id_taxa_clean = furrr::future_map(id_taxa, clean_id_taxa, taxa)) %>%
#     dplyr::select(asv, id_taxa_clean) %>%
#     tidyr::unnest(cols = c(id_taxa_clean)) %>%
#     dplyr::select(asv, tidyselect::one_of(taxa)) %>%
#     dplyr::right_join(asvs, by = "asv")

# }

# kraken_results %<>%
#   dplyr::mutate( out = map(out, clean_id_taxa_wrap, taxa, asvs))

# kraken_results %<>%
#   mutate(
#     total = map_int(out, nrow),
#     labelled = map(out, na.omit) %>% map_int(nrow),
#     plabelled = labelled / total)

# dir.create(file.path(basedr, "rds", "kraken"), showWarnings = FALSE)

# kraken_results %>%
#   dplyr::select_if(negate(is.list)) %>%
# 	saveRDS(file.path(basedr, "dust_comp", "kraken",
#     "2020_02_17_kraken_summary.rds"))


# outdr <- file.path(basedr, "dust_comp", "kraken")

# outfiles <- str_c("dustcomp_", pull(kraken_results,db), "_taxa.rds")
# outfiles %<>% file.path(outdr, .)

# walk2( pull(kraken_results, out), outfiles, saveRDS)
