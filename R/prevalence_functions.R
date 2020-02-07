#' @importFrom forcats fct_explicit_na
#' @importFrom tibble as_tibble
#' @importFrom tibble tibble
NULL

#' function that checks if any element is na in a vector
#' @param x a vector
#' @return logical value if there is at least one \code{NA} value
#' @export
any_na <- function(x)any(is.na(x))


#' computes the prevalence statistics for each ASV
#'
#' @param asv_table matrix with samples as rows and ASVs as columns
#' @param idtaxa matrix with taxa classification for each ASV
#' @return a \code{tibble} with the prevalence statistics per ASV
#' @export
compute_prevalence <- function(asv_table, idtaxa) {

  idtaxa <- tibble::as_tibble(idtaxa)

  prevalence <- tibble::tibble(
    id = stringr::str_c("asv", seq_len(ncol(asv_table)), sep = "_"),
    seq = colnames(asv_table),
    prev = apply(asv_table, 2, FUN = function(x)sum(x  > 0)),
    nreads = colSums(asv_table),
    ave_reads = colMeans(asv_table))

  idtaxa <- dplyr::mutate(idtaxa,
    id = stringr::str_c("asv", seq_len(ncol(asv_table)), sep = "_"))
  idtaxa <- dplyr::mutate_if(idtaxa, any_na,
    list(~ forcats::fct_explicit_na(., "unlabelled")))

  dplyr::inner_join(prevalence, idtaxa, by = "id")
}

#' plots the prevalence statistic
#'
#' @param prevalence tibble with prevalence statistic
#'  computed by \code{compute_prevalence}
#' @param nsamples integer with the number of samples in the data
#' @return a \code{ggplot} object with the plot
#' @export
plot_prevalence <- function(prevalence, nsamples) {
  nreads <- NULL
  prev <- NULL
  to_check <- NULL

  if (!any(names(prevalence) == "to_check")) {
    prevalence <- dplyr::mutate(prevalence = "no")
  }

  prev_plot <- ggplot2::ggplot(prevalence) +
    ggplot2::geom_point(
      ggplot2::aes(nreads, prev / nsamples, colour = to_check), alpha = .7) +
    ggplot2::scale_x_log10(labels = scales::comma_format(accuracy = 1)) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(
      colour = "reads in \n single asv",
      x = "total abundance",
      y = "prevalence") +
    ggplot2::facet_wrap(. ~ phylum) +
    ggplot2::scale_color_manual(
      values = c(`yes` = "red", `no` = "navyblue")) +
    ggplot2::theme(
      legend.position = "top",
      strip.background = ggplot2::element_blank())

}
