##' @import furrr
##' @import future
##' @importFrom tidyr unnest
NULL

##' summarizes the number of reads resulted after each step into one table
##'
##' @param reads_file, a 3 column csv file with columns: sample | R1.fastq | R2.fastq
##' @param outprefix prefix used in the final ASV table
##' @param outdir names of the output directory where the results are saved
##' @param cores number of cpus used to run the process
##' @export
summarize_number_reads <- function(reads_file,outprefix,outdir,cores)
{
  `.` <- R1 <- R2 <- NULL
  name <- fold_change <- reads.out <- reads.in <- reads.merged_pairs <- NULL
  reads.asv_table <- NULL

  future::plan(future::multiprocess,workers = cores)

  stopifnot(
    file.exists(reads_file),
    dir.exists(outdir))

  ### read input file, i.e. the 3 column file
  input_files = read_csv(reads_file,col_names = FALSE) %>% set_names(c("name","R1","R2"))

  input_files %<>%
    mutate(
      summ = file.path(outdir,"filter_fastq_summary",paste0(name,"_trim_summary.rds")) %>%
        furrr::future_map(readRDS)) %>%
    tidyr::unnest() %>%
    select(-fold_change)

  input_files %<>% select(-R1,-R2) %>%
    mutate(
      reads.merged_pairs = file.path(outdir,"merged_pairs",
                                   paste0(name,"_merged_pairs.rds")) %>%
        furrr::future_map(readRDS) %>%
        furrr::future_map_dbl( ~ sum(.$abundance)))

  asv_table <- readRDS(file.path(outdir,"ASV_tables",paste0(outprefix,"_sequence_table.rds")))

  asv_table %<>%
    {
      rs = rowSums(.)
      tibble(
        name = names(rs),
        reads.asv_table = rs)}

  input_files %<>% inner_join(asv_table,by = "name")
  input_files %>%
    mutate(
      perc.out = reads.out / reads.in,
      perc.merged_pairs = reads.merged_pairs / reads.in,
      perc.asv_table = reads.asv_table / reads.in  )

}
