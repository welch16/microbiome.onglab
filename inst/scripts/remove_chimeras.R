#!/usr/bin/env Rscript

library(optparse,quietly = TRUE)

optList <-  list(
  make_option("--dada_file", action = "store_true", type = "character",
              help = "A file with the name of all the dada2 tables to merge"),
  make_option("--param_file", action = "store_true",type = "character",
              help = "json file with the parameters used to remove bimeras"),
  make_option("--outprefix", action = "store_true", type = "character",
              default = "dada2",
              help = "Name of the output file with the labelled ASVs after
                    	removing the bimeras, the full file name will
        				be {outdir}/ASV_tables/{outprefix}_asv_wo_bimeras.rds"),
  make_option("--outdir", action = "store_true", type = "character",
              default = tempdir(),
              help = "Location of the output directory"),
  make_option("--cores", action = "store_true", type = "numeric",
              default = 4,
              help = "Number of parallel cpus to use"))

opt <- parse_args(OptionParser(option_list = optList))

options(mc.cores = opt$cores)

out_file = file.path(opt$outdir, "ASV_tables",paste0(opt$outprefix,"_asv_wo_bimeras.rds"))

# opt$dada_file = "./condor_runs/2019_06_25_full_dust_ASV_tables.txt"
# opt$param_file = "./dada2_pipeline_params/2019_06_21_dada2_bimera_parameters.json"
# opt$taxa_model = "../idtaxa_trained_models/GTDB_r86-mod_September2018.RData"
#


stopifnot(file.exists(opt$dada_file))

library(magrittr)
library(tidyverse)
library(dada2)
library(microbiome.onglab)
library(jsonlite)

if(!file.exists(opt$param_file)){
  params <- data.frame()
}else{
  params <- fromJSON(opt$param_file, flatten = TRUE)
}

dada2_params <- c("minSampleFraction","ignoreNNegatives",
                  "minFoldParentOverAbundance","allowOneOf","minOneOffParentDistance","maxShift")

stopifnot( all( names(params) %in% dada2_params))

message("Removing bimeras")
message("Input files in: ", opt$dada_file)
message("Output file: ", out_file)

if(!file.exists(out_file)){

  asv_tables <- read_delim(opt$dada_file, " ",col_names = FALSE) %>% pull(X1)
  asv_tables %<>% map(readRDS)

  if(length(asv_tables) > 1){
    asv_table <- mergeSequenceTables(tables = asv_tables)
  }else{
    asv_table <- asv_tables[[1]]
  }

  asv_table <- removeBimeraDenovo(
    asv_table, method = "consensus",
    minSampleFraction = get_param_chimeras("minSampleFraction",params),
    ignoreNNegatives = get_param_chimeras("ignoreNNegatives",params),
    minFoldParentOverAbundance = get_param_chimeras("minFoldParentOverAbundance",params),
    allowOneOf = get_param_chimeras("allowOneOf",params),
    minOneOffParentDistance = get_param_chimeras("minOneOffParentDistance",params),
    maxShift = get_param_chimeras("maxShift",params),
    multithread = TRUE)

  saveRDS(asv_table,out_file)


}


