
#' auxiliar function to label the taxa
#'
#' @param dna \code{DNAStringSet} with the ASV sequences
#' @param rdata_file idtaxa model learned saved as a RData file
#' @param confidence_thr confidence threshold, higher levels
#'  indicates higher accuracy
#' @param cores number of cpus to use
#' @return a \code{data.frame} with the labelled taxa
#' @export
get_tax_ids <- function(dna, rdata_file, confidence_thr = 60, cores = NULL) {
  training_set <- NULL
  load(rdata_file)
  ids <- IdTaxa(dna, training_set, strand = "both",
                threshold = confidence_thr,
                processors = cores, verbose = FALSE) # use all processors
  ranks <- c("domain", "phylum", "class", "order",
    "family", "genus", "species") # ranks of interest
  # Convert the output object of class "Taxa" to a matrix
  # analogous to the output from assignTaxonomy
  taxid <- t(sapply(ids, function(x) {
    m <- match(ranks, x$rank)
    taxa <- x$taxon[m]
    taxa[startsWith(taxa, "unclassified_")] <- NA
    taxa
  }))
  colnames(taxid) <- ranks
  rownames(taxid) <- as.character(dna)
  taxid
}

#' generate the input file to label taxa with condor
#'
#' @param asv_file rds file with the ASV table without chimeras
#' @param taxa_model_file RData file with ta learned taxa model to use
#'  with \code{DECIPHER}
#' @param taxa_prefix name of the output file
#' @param taxa_outdir directory where the output is going to be saved
#' @param condor_file name of the file with condor instructions
#' @param batch_name string with the name of the batch
#' @param confidence_thr confidence threshold, higher levels
#'  indicates higher accuracy
#' @param request_cores number of cpus per machine
#' @param request_mem number of GB required as memory
#' @export
condor_label_taxa <- function(
  asv_file = "",
  taxa_model_file = "",
  taxa_prefix = "taxa",
  taxa_outdir = ".",
  condor_file = "condor_label_taxa",
  batch_name = "dada2_label_taxa",
  confidence_thr = 60,
  request_cores = 4,
  request_mem = "4 GB"
) {

  if (!file.exists(taxa_model_file)) {
    warning("need to define taxa_model_file")
  }

  warning_file(asv_file, "Need to define the remove chimeras file")

  file_connection <- file(condor_file)
  rscript <- system("which Rscript", intern = TRUE)

  writeLines(
    c(
      "universe         = vanilla",
      stringr::str_c("batch_name       = ", batch_name),
      stringr::str_c("executable       = ", rscript),
      stringr::str_c("args             = $(script_r) --asv_file $(infile)",
        "--taxa_model $(taxamodel) --thr $(thr) --outprefix $(outprefix)",
        "--outdir $(outdir) --cores ", request_cores, sep = " "),
      stringr::str_c("request_cpus     = ", request_cores),
      stringr::str_c("request_memory   = ", request_mem),
      "on_exit_hold     = (ExitBySignal == True) || (ExitCode != 0)",
      stringr::str_c("periodic_release = (NumJobStarts < 5) &&",
        "((CurrentTime - EnteredCurrentStatus) > 180)", sep = " "),
      stringr::str_c("script_r         = ",
        system.file("scripts/label_taxa.R", package = "microbiome.onglab")),
      stringr::str_c("infile           = ", asv_file),
      stringr::str_c("outdir           = ", taxa_outdir),
      stringr::str_c("output           = $(outdir)/out/dada2_label_taxa_",
        "$(outprefix).$(cluster).$(process).out"),
      stringr::str_c("error            = $(outdir)/err/dada2_label_taxa_",
        "$(outprefix).$(cluster).$(process).err"),
      stringr::str_c("log              = $(outdir)/log/dada2_label_taxa_",
        "$(outprefix).$(cluster).$(process).log"),
      stringr::str_c("thr              = ", confidence_thr),
      stringr::str_c("taxamodel        = ", taxa_model_file),
      stringr::str_c("outprefix        = ", taxa_prefix),
      "queue 1"), file_connection)

  close(file_connection)
}
