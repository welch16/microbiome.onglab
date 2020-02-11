#!/usr/bin/env Rscript

info <- Sys.info();
message(stringr::str_c(names(info), " : ", info, "\n"))

# wrapper for running sourcetracker as QC
#

library(optparse)

opt_list <- list(
  optparse::make_option("--asv_file", action = "store_true",
    type = "character",
    help = "Name of the dada2 asv matrix. This is the matrix after
      removing the chimeras."),
  optparse::make_option("--neg_control_file", action = "store_true",
    type = "character",
    help = "Name of an rds file with the relationship between samples and
      controls. It is expected a tibble with two columns: sample | neg_control,
      where the elements of the second column are character vectors with the
      negative controls that correspond to the sample"),
  optparse::make_option("--param_file", action = "store_true",
    type = "character",
    help = "Location of the parameter file in json format, if not provided
      the script will use sourcetracker default parameters"),
  optparse::make_option("--top_asvs", action = "store_true", type = "numeric",
    default = 50,
    help = "Number of the top ASVs to use"),
  optparse::make_option("--max_split_size", action = "store_true",
    type = "numeric", default = 30,
    help = "Max number of samples that are going to be processed at once"),
  optparse::make_option("--rarefaction_depth", action = "store_true",
    type = "character", default = "inf",
    help = "Depth to rarefy the samples [Numeric or 'inf']"),
  optparse::make_option("--outfile", action = "store_true", type = "character",
    help = "Name of the file where the results are saved"),
  optparse::make_option("--cores", action = "store_true", type = "numeric",
    default = 4,
    help = "Number of cores to be used for parallel processing"))


opt <- optparse::parse_args(optparse::OptionParser(option_list = opt_list))

out_file <- opt$outfile ## may change if we change the directory structure

stopifnot(file.exists(opt$asv_file),
          file.exists(opt$neg_control_file))

if (tolower(opt$rarefaction_depth) == "inf") {
  opt$rarefaction_depth <- NULL
}else{
  opt$rarefaction_depth <- as.numeric(opt$rarefaction_depth)
}

library(magrittr)
library(tidyverse)
library(microbiome.onglab)
library(jsonlite)
library(BiocParallel)


asv_table <- readRDS(opt$asv_file)
neg_controls <- readRDS(opt$neg_control_file)

if (!file.exists(opt$param_file)) {
  params <- data.frame()
}else{
  params <- jsonlite::fromJSON(opt$param_file, flatten = TRUE)
}

dada2_params <- c("alpha1", "alpha2", "beta")

stopifnot(all(names(params) %in% dada2_params))

neg_controls_groups <- neg_controls %>%
  dplyr::mutate(
    neg_control_full = purrr::map_chr(neg_control, stringr::str_c,
      collapse = "_")) %>%
  dplyr::count(neg_control_full) %>%
  dplyr::mutate(
    id = str_c("G", dplyr::row_number()))

neg_controls %<>%
  dplyr::mutate(
    neg_control_full = purrr::map_chr(neg_control, stringr::str_c,
      collapse = "_")) %>%
  dplyr::left_join(neg_controls_groups, by = "neg_control_full") %>%
  dplyr::select(-n, -neg_control_full)

samples <- rownames(asv_table)

all_negative_controls <- neg_controls %>%
  dplyr::pull(neg_control) %>%
  unlist() %>%
  unique()

## only use the top_asvs in order to make the code run faster
asv_table <- asv_table[, seq_len(opt$top_asvs)]


microbiome.onglab::load_sourcetracker()
alpha1 <- microbiome.onglab::get_param_sourcetracker("alpha1", params)
alpha2 <- microbiome.onglab::get_param_sourcetracker("alpha2", params)
beta <- microbiome.onglab::get_param_sourcetracker("beta", params)

neg_controls %<>% tidyr::nest(-id, .key = "group")

create_sourcetracker <- function(group, asv, rarefaction_depth) {

  samples <- group %>% dplyr::pull(sample)
  neg_controls <- dplyr::pull(group, neg_control) %>%
    unlist() %>%
    unique()

  mat <- microbiome.onglab::aux_get_matrix(asv, samples)
  # sourcetracker is a function obtained from evaluation the source code
  # obtained from RCurl
  sourcetracker(mat, neg_controls, rarefaction_depth = rarefaction_depth)
}

# train SourceTracker object on training data
message("creating sourcetracker objects")
neg_controls %<>%
  dplyr::mutate(
    stracker = purrr::map(group, create_sourcetracker,
      asv_table, opt$rarefaction_depth))

message("separating large samples")
neg_controls %<>%
    dplyr::mutate(
      st_pred = purrr::map(group,
        microbiome.onglab::aux_split_run, asv_table, opt$max_split_size))

results <- neg_controls %>%
  dplyr::select(id, st_pred) %>%
  tidyr::unnest(cols = c(st_pred)) %>%
  dplyr::left_join(
    dplyr::select(neg_controls, id, stracker), by = "id")

# fix st pred, the issue was that when splitting into groups, there was one
# group with one sample, which caused problems when merging the results
check_if_matrix <- function(pred) {

  if (class(pred) != "matrix") {
    pred %<>% as.matrix() %>% t()
  }
  pred
}

# for the few cases where the submatrix is actually a vector, needs to convert
# to matrix again
results %<>%
  dplyr::mutate(
    st_pred = purrr::map(st_pred, check_if_matrix))

message("processing samples with sourcetracker mixture model")
results %<>%
  dplyr::mutate(
    pred = BiocParallel::bpmapply(function(x, y)
      predict(x, y, alpha1 = alpha1, alpha2 = alpha2, beta = beta,
      verbosity = TRUE), stracker, st_pred,
      BPPARAM = BiocParallel::MulticoreParam(workers = opt$cores),
      SIMPLIFY = FALSE))

nest <- tidyr::nest_legacy

# merge parallel results and clean them to tibble / data.frame
message("merging results")
results %<>%
  dplyr::select(id, pred) %>%
  nest(pred, .key = "preds") %>%
  dplyr::mutate(
    preds = purrr::map(preds, pull, pred),
    results = purrr::map(preds, microbiome.onglab::merge_sourcetracker_results),
    draws = purrr::map(results, "draws"),
    proportions = purrr::map(results, "proportions"),
    proportions_sd = purrr::map(results, "proportions_sd"),
    train_envs = purrr::map(results, "train.envs"),
    samplenames = purrr::map(results, "samplenames"))

clean_proportions <- function(proportions, proportions_sd) {
  proportions %<>% as.data.frame() %>%
    tibble::rownames_to_column("sample") %>%
    tibble::as_tibble()
  proportions_sd %<>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample") %>%
    tibble::as_tibble()
  proportions %<>%
    tidyr::gather(control, value, -sample)
  proportions_sd %<>%
    tidyr::gather(control, value, -sample)

  dplyr::bind_rows(
    proportions %>%
      dplyr::mutate(
        type = "mean"),
    proportions_sd %>%
      dplyr::mutate(
        type = "sd"))

}

clean_draws <- function(draw, train_envs, samplenames) {

  draw <- aperm(draw, perm = c(3, 2, 1))
  train_envs %<>% stringr::str_replace_all("\\-", "\\_")
  samplenames %<>% stringr::str_replace_all("\\-", "\\_")

  dimnames(draw) <- list(samplenames,
                        train_envs,
                        NULL)

  flat_draw <- draw %>%
    as.data.frame() %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      sample = samplenames) %>%
    tidyr::gather(control_draw, proportion, -sample)

  flat_draw %>%
    tidyr::separate(control_draw,
      into = c("control", "draw"), sep = "\\.")

}

results %<>%
  dplyr::left_join(
    dplyr::select(neg_controls, id, group),
    by = "id") %>%  ## brief correction for groups of length 1
  dplyr::mutate(
    samplenames2 = purrr::map(group, pull, sample),
    samplenames = dplyr::if_else(
      purrr::map_lgl(samplenames, is.null),
        samplenames2, samplenames)) %>%
  dplyr::select(-samplenames2) %>%
  dplyr::mutate(
    draws = purrr::pmap(list(draws, train_envs, samplenames), clean_draws),
    proportions = purrr::map2(proportions, proportions_sd, clean_proportions))

results %<>%
  dplyr::select(
    -proportions_sd,
    -train.envs,
    -samplenames,
    -preds,
    -results) %>%
  dplyr::select(id, group, draws, proportions)

results %>% saveRDS(out_file)
