#' @import future
#' @import furrr
#' @import ggplot2
#' @importFrom tidyr unnest
#' @importFrom tidyr gather
#' @importFrom forcats fct_relevel
#' @importFrom stats median
NULL

#' summarizes the number of reads resulted after each step into one table
#'
#' @param reads_file, a 3 column csv file with columns:
#'   sample | R1.fastq | R2.fastq
#' @param outprefix prefix used in the final ASV table
#' @param outdir names of the output directory where the results are saved
#' @param cores number of cpus used to run the process
#' @export
summarize_number_reads <- function(reads_file, outprefix, outdir, cores) {
  `.` <- end1 <- end2 <- NULL
  name <- fold_change <- reads.out <- NULL
  reads.in <- reads.merged_pairs <- NULL
  reads.asv_table <- NULL

  future::plan(future::multiprocess, workers = cores)

  stopifnot(
    file.exists(reads_file),
    dir.exists(outdir))

  ### read input file, i.e. the 3 column file
  input_files <- read_csv(reads_file, col_names = FALSE) %>%
    set_names(c("name", "end1", "end2"))

  input_files %<>%
    dplyr::mutate(
      summ = file.path( outdir, "filter_fastq_summary",
        stringr::str_c(name, "_trim_summary.rds")) %>%
        furrr::future_map(readRDS)) %>%
    tidyr::unnest() %>%
    dplyr::select(-fold_change)

  input_files %<>% select(-R1,-R2) %>%
    mutate(
      reads.merged_pairs = file.path(outdir,"merged_pairs",
                                   paste0(name,"_merged_pairs.rds")) %>%
        furrr::future_map(readRDS) %>%
        furrr::future_map_dbl( ~ sum(.$abundance)))

  asv_table <- readRDS(file.path(outdir,"ASV_tables",paste0(outprefix,"_sequence_table.rds")))

  asv_table %<>% {
    rs <- rowSums(.)
    tibble::tibble(
      name = names(rs),
      reads.asv_table = rs)}

  input_files %<>% inner_join(asv_table,by = "name")
  input_files %>%
    mutate(
      perc.out = reads.out / reads.in,
      perc.merged_pairs = reads.merged_pairs / reads.in,
      perc.asv_table = reads.asv_table / reads.in )

}


##' compares the (relative) abundance for each step
##' @param nreads_summary output of the \code{summarize_number_reads} function
##' @param summary_fun function used to summarize the trend, the default function is the \code{median}
##' @param relative boolean indicator determining wheter comparing the relative abundance
##' @return a \code{ggplot} object
##' @export
plot_abundance_per_step <- function(nreads_summary, summary_fun = median, relative = FALSE)
{
  name <- step <- sample <- NULL
  if(relative){

    nreads_summary <- select(nreads_summary,name,sample,contains("perc"))
    nreads_summary <- mutate(nreads_summary,perc.in = 1)

  }else{

    nreads_summary <- select(nreads_summary,name,sample,contains("reads"))

  }

  nreads_summary <- tidyr::gather(nreads_summary,step,value,-name,-sample)
  nreads_summary <- mutate(nreads_summary,
      step = stringr::str_remove(step,"reads."),
      step = stringr::str_remove(step,"perc."),
      step = factor(step) %>% fct_relevel("in","out","merged_pairs"))

  fun_summary <- group_by(nreads_summary,step)
  fun_summary <- summarize(fun_summary, value = summary_fun(value))

  base_plot <- nreads_summary %>%
    ggplot(aes(step,value))+
    geom_boxplot()+
    geom_point(alpha = 1/4,shape = 21)+
    geom_line(aes(group = name),alpha = 1/4)+
    geom_line(data = fun_summary,size = 1.5,aes(group = 1), linetype = 2,
              colour = "red")+
    theme(
      legend.position = "none",
      strip.text.y = element_text(angle = -90, size =10),
      strip.background = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank()
    )+
    labs(
      x = "processing step",
      y = str_c(if_else(relative,"relative",""), "abundance after step",sep =" "))

  if(relative){
    base_plot +
      scale_y_continuous(labels = scales::percent_format())
  }else{
    base_plot +
      scale_y_continuous(labels = scales::comma_format(scale = 1e-3,suffix = "K"))
  }




}
