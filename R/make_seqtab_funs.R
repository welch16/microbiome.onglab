
##' generate the input to make the sequence table files with condor
##'
##' @param make_queue_file name of the file with the queue, it has 3 columns sample_name | R1.fastq | R2.fastq
##' @param make_prefix prefix name for the sequence table output file
##' @param make_outdir directory where the output is going to be saved
##' @param condor_file name of the file with condor instructions
##' @param batch_name string with the name of the batch
##' @param request_cores number of cpus per machine
##' @param request_mem number of GB required as memory
##' @export
condor_make_sequence_table <- function(
  make_queue_file = "",
  make_prefix = "seqtab",
  make_outdir = ".",
  condor_file = "./condor_make_seqtab",
  batch_name = "dada2_make_seqtab",
  request_cores = 4,
  request_mem = "4 GB"
)
{

  str_c <- stringr::str_c

  warning_file(make_queue_file, "Need to define make sequence tab queue file")

  file_connection <- file(condor_file)

  rscript <- system("which Rscript", intern = TRUE)

  writeLines(
    c(
    "universe         = vanilla",
    str_c("batch_name       = ", batch_name),
    str_c("executable       = ", rscript),
    str_c("args             = $(script_r) --queue_file $(queue_file) --outprefix $(outprefix) --outdir $(outdir) --cores ", request_cores),
    str_c("request_cpus     = ", request_cores),
    str_c("request_memory   = ", request_mem),
    "on_exit_hold     = (ExitBySignal == True) || (ExitCode != 0)",
    "periodic_release = (NumJobStarts < 5) && ((CurrentTime - EnteredCurrentStatus) > 180)",
    "basedr           = .",
    str_c("script_r         = ", system.file("scripts/make_sequence_table.R", package = "microbiome.onglab")),
    str_c("outdir           = ", make_outdir),
    "output           = $(outdir)/out/dada2_seqtab_$(outprefix).$(cluster).$(process).out",
    "error            = $(outdir)/err/dada2_seqtab_$(outprefix).$(cluster).$(process).err",
    "log              = $(outdir)/log/dada2_seqtab_$(outprefix).$(cluster).$(process).log",
    str_c("queue_file        = ", make_queue_file),
    str_c("outprefix        = ", make_prefix),
    "queue 1"), file_connection)


  close(file_connection)
}

##' function to add default parameters to remove_chimeras script
##'
##' @param param_name name of the parameter to look for
##' @param param_frame data.frame with the chosen parameters
##' @return the option in param_frame of the default value for
##'     \code{param_name}
##' @export
get_param_chimeras <- function(param_name, param_frame){

  out <- param_frame[[param_name]]

  if(is.null(out)){

    if(param_name %in% c("minSampleFraction")) out <- 0.9
    if(param_name %in% c("ignoreNNegatives")) out <- 1
    if(param_name %in% c("minFoldParentOverAbundance")) out <- 1.5
    if(param_name %in% c("allowOneOf")) out <- FALSE
    if(param_name %in% c("minOneOffParentDistance")) out <- 4
    if(param_name %in% c("maxShift")) out <- 16

  }else{
    if(class(out) == "list"){
      out <- out[[1]]
    }else{
      if(out == "Inf") out <- Inf
      if(out == "FALSE") out <- FALSE
    }
  }
  out
}


##' generate the input file to remove chimeras from the sequence table with condor
##'
##' @param sequence_table_file name of the file with the sequence tables to merge (one p/line)
##' @param chim_param_file name of the remove chimeras parameter file
##' @param chim_prefix prefix name for the sequence table output file
##' @param chim_outdir directory where the output is going to be saved
##' @param condor_file name of the file with condor instructions
##' @param batch_name string with the name of the batch
##' @param request_cores number of cpus per machine
##' @param request_mem number of GB required as memory
##' @export
condor_remove_chimeras <- function(
  sequence_table_file,
  chim_param_file,
  chim_prefix,
  chim_outdir,
  condor_file,
  batch_name,
  request_cores,
  request_mem
)
{

  str_c <- stringr::str_c
  warning_file(sequence_table_file, "Need to define the remove chimeras file")
  warning_file(chim_param_file, "Need to define the parameter file for removing chimeras")

  file_connection <- file(condor_file)

  rscript <- system("which Rscript", intern = TRUE)

  writeLines(
    c(
      "universe         = vanilla",
      str_c("batch_name       = ", batch_name),
      str_c("executable       = ", rscript),
      str_c("args             = $(script_r) --dada_file $(infile) --param_file $(param) --outprefix $(outprefix) --outdir $(outdir) --cores ", request_cores),
      str_c("request_cpus     = ", request_cores),
      str_c("request_memory   = ", request_mem),
      "on_exit_hold     = (ExitBySignal == True) || (ExitCode != 0)",
      "periodic_release = (NumJobStarts < 5) && ((CurrentTime - EnteredCurrentStatus) > 180)",
      str_c("script_r         = ", system.file("scripts/remove_chimeras.R", package = "microbiome.onglab")),
      str_c("infile           = ", sequence_table_file),
      str_c("param            = ", chim_param_file),
      str_c("outdir           = ", chim_outdir),
      str_c("outprefix        = ", chim_prefix),
      "output           = $(outdir)/out/dada2_remove_chimera_$(outprefix).$(cluster).$(process).out",
      "error            = $(outdir)/err/dada2_remove_chimera_$(outprefix).$(cluster).$(process).err",
      "log              = $(outdir)/log/dada2_remove_chimera_$(outprefix).$(cluster).$(process).log",
      "queue 1"), file_connection)

  close(file_connection)

}

get_tax_ids <- function(dna,rdata_file)
{
  trainingSet <- NULL
  load(rdata_file)
  ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
  ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
  # Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
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

##' generate the input file to label taxa with condor
##'
##' @param asv_file rds file with the ASV table without chimeras
##' @param taxa_model_file RData file with ta learned taxa model to use with \code{DECIPHER}
##' @param taxa_prefix name of the output file
##' @param taxa_outdir directory where the output is going to be saved
##' @param condor_file name of the file with condor instructions
##' @param batch_name string with the name of the batch
##' @param request_cores number of cpus per machine
##' @param request_mem number of GB required as memory
##' @export
condor_label_taxa <- function(
  asv_file,
  taxa_model_file,
  taxa_prefix,
  taxa_outdir,
  condor_file,
  batch_name,
  request_cores,
  request_mem = "4 GB"
)
{
  str_c <- stringr::str_c

  stopifnot(file.exists(taxa_model_file))

  warning_file(asv_file, "Need to define the remove chimeras file")

  if(!dir.exists(taxa_outdir))dir.create(taxa_outdir,showWarnings = FALSE)

  dir.create(file.path(taxa_outdir,"err"),showWarnings = FALSE)
  dir.create(file.path(taxa_outdir,"log"),showWarnings = FALSE)
  dir.create(file.path(taxa_outdir,"out"),showWarnings = FALSE)

  file_connection <- file(condor_file)

  rscript <- system("which Rscript", intern = TRUE)

  writeLines(
    c(
      "universe         = vanilla",
      str_c("batch_name       = ", batch_name),
      str_c("executable       = ", rscript),
      str_c("args             = $(script_r) --asv_file $(infile) --taxa_model $(taxamodel) --outprefix $(outprefix) --outdir $(outdir) --cores ", request_cores),
      str_c("request_cpus     = ", request_cores),
      str_c("request_memory   = ", request_mem),
      "on_exit_hold     = (ExitBySignal == True) || (ExitCode != 0)",
      "periodic_release = (NumJobStarts < 5) && ((CurrentTime - EnteredCurrentStatus) > 180)",
      str_c("script_r         = ", system.file("scripts/label_taxa.R",package = "microbiome.onglab")),
      str_c("infile           = ", asv_file),
      str_c("outdir           = ", taxa_outdir),
      "output           = $(outdir)/out/dada2_label_taxa_$(outprefix).$(cluster).$(process).out",
      "error            = $(outdir)/err/dada2_label_taxa_$(outprefix).$(cluster).$(process).err",
      "log              = $(outdir)/log/dada2_label_taxa_$(outprefix).$(cluster).$(process).log",
      str_c("taxamodel        = ", taxa_model_file),
      str_c("outprefix        = ", taxa_prefix),
      "queue 1"), file_connection)

  close(file_connection)

}

