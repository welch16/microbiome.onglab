#' @import future
#' @importFrom magrittr `%<>%`
#' @importFrom magrittr `%>%`
#' @importFrom cowplot plot_grid
#' @import ggplot2
#' @importFrom purrr map
#' @importFrom purrr map2
#' @import furrr
#' @import tidyselect
#' @importFrom tidyr unnest
#' @importFrom tidyr gather
#' @importFrom forcats fct_relevel
#' @importFrom stats median
#' @importFrom grDevices dev.off
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
  name <- fold_change <- NULL

  future::plan(future::multiprocess, workers = cores)

  stopifnot(file.exists(reads_file), dir.exists(outdir))
  ### read input file, i.e. the 3 column file
  input_files <- read_csv(reads_file, col_names = FALSE) %>%
    set_names(c("name", "end1", "end2"))

  input_files %<>%
    dplyr::mutate(
      summ = file.path(outdir, "filter_fastq_summary",
        stringr::str_c(name, "_trim_summary.rds")) %>%
        furrr::future_map(readRDS)) %>%
    tidyr::unnest(cols = c(summ)) %>%
    dplyr::select(-fold_change)

  input_files %<>%
    dplyr::select(-end1, -end2) %>%
    dplyr::mutate(
      reads.merged_pairs = file.path(outdir, "merged_pairs",
                                   paste0(name, "_merged_pairs.rds")) %>%
        furrr::future_map(readRDS) %>%
        furrr::future_map_dbl(~ sum(.$abundance)))

  asv_table <- readRDS(
    file.path(outdir, "ASV_tables",
      stringr::str_c(outprefix, "_sequence_table.rds")))

  asv_table %<>% {
    rs <- rowSums(.)
    tibble::tibble(
      name = names(rs),
      reads.asv_table = rs)
  }

  input_files %<>% dplyr::inner_join(asv_table, by = "name")
  input_files %>%
    dplyr::mutate(
      perc.out = reads.out / reads.in,
      perc.merged_pairs = reads.merged_pairs / reads.in,
      perc.asv_table = reads.asv_table / reads.in)

}


#' compares the (relative) abundance for each step
#' @param nreads_summary output of the \code{summarize_number_reads} function
#' @param summary_fun function used to summarize the trend,
#'  the default function is the \code{median}
#' @param relative boolean indicator determining wheter comparing the
#'  relative abundance
#' @return a \code{ggplot} object
#' @export
plot_abundance_per_step <- function(
  nreads_summary,
  summary_fun = median,
  relative = FALSE) {
  name <- step <- sample <- NULL

  if (relative) {

    nreads_summary <- dplyr::select(nreads_summary,
      name, sample, tidyselect::contains("perc"))
    nreads_summary <- dplyr::mutate(nreads_summary, perc.in = 1)

  }else{
    nreads_summary <- dplyr::select(nreads_summary,
      name, sample, tidyselect::contains("reads"))
  }

  nreads_summary <- tidyr::gather(nreads_summary,
    step, value, -name, -sample)
  nreads_summary <- mutate(nreads_summary,
      step = stringr::str_remove(step, "reads."),
      step = stringr::str_remove(step, "perc."),
      step = factor(step) %>% fct_relevel("in", "out", "merged_pairs"))

  fun_summary <- dplyr::group_by(nreads_summary, step)
  fun_summary <- dplyr::summarize(fun_summary, value = summary_fun(value))

  base_plot <- nreads_summary %>%
    ggplot2::ggplot(aes(step, value)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_point(alpha = 0.25, shape = 21) +
    ggplot2::geom_line(aes(group = name), alpha = 0.25) +
    ggplot2::geom_line(data = fun_summary, size = 1.5,
      ggplot2::aes(group = 1), linetype = 2, colour = "red") +
    ggplot2::theme(
      legend.position = "none",
      strip.text.y = ggplot2::element_text(angle = -90, size = 10),
      strip.background = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank()) +
    ggplot2::labs(
      x = "processing step",
      y = stringr::str_c(
        dplyr::if_else(relative, "relative", ""),
          "abundance after step", sep = " "))

  if (relative) {
    base_plot +
      ggplot2::scale_y_continuous(labels = scales::percent_format())
  }else{
    base_plot +
      ggplot2::scale_y_continuous(
        labels = scales::comma_format(scale = 1e-3, suffix = "K"))
  }
}


#' plots the quality profiles generated by dada2's
#' \code{plotQualityProfile} function, both ends per page, and all samples in a
#' single pdf file
#' @param sample_tibble a tibble / data.frame with 3 columns:
#'   sample_name | end 1 file | end 2 file
#' @param pdf_filename a string with the name of the file where the plots are
#'   saved
#' @export
plot_quality_profiles <- function(sample_tibble, pdf_filename) {

  if (ncol(sample_tibble) > 3) {
    warning("sample_tibble has more than 3 columns, using the first 3")
  }
  vars <- names(sample_tibble)
  sample_tibble %<>% dplyr::select(tidyselect::one_of(vars[1:3])) %>%
    magrittr::set_names(c("name", "end1", "end2"))

  profile_plot <- function(file) {
    dada2::plotQualityProfile(file) +
      ggplot2::theme(strip.background = ggplot2::element_blank())
  }

  sample_tibble %<>%
    dplyr::mutate(
      end1_plot = purrr::map(end1, profile_plot),
      end2_plot = purrr::map(end2, profile_plot))

  plot_all <- function(end1, end2, filename) {
    grDevices::pdf(file = filename, width = 9, height = 4)
    u <- purrr::map2(end1, end2, cowplot::plot_grid, nrow = 1)
    u <- purrr::map(u, print)
    grDevices::dev.off()
  }

  plot_all(sample_tibble$end1_plot, sample_tibble$end2_plot,
    pdf_filename)

}