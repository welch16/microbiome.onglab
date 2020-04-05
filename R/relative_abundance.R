#' @importFrom tidyr gather
#' @importFrom dplyr group_by summarize ungroup mutate select filter
#' @importFrom Polychrome green.armytage.colors dark.colors
#' @importFrom ggplot2 ggplot aes geom_col facet_grid
#' @importFrom ggplot2 scale_y_continuous theme scale_fill_manual
#' @importFrom ggplot2 guides guide_legend
#' @importFrom grDevices colorRampPalette
NULL

#' summarizes the asv table (in tibble format) grouping by a taxonomy level
#' @param asv_tibble a tibble with asv information merged with taxonomy
#'  annotations
#' @param taxa_rank a string with the taxonomic rank for which the data
#'  is summarized
#' @return a tibble with the abundance and prevalence for each taxa
#' @export
summarize_per_taxa <- function(asv_tibble, taxa_rank) {

  abu <- abundance <- prevalence <- NULL

  samples <- names(asv_tibble)

  asv_tibble %>%
    tidyr::gather(sample, abu, all_of(samples)) %>%
    dplyr::group_by(!! rlang::sym(taxa_rank), sample) %>%
    dplyr::summarize(
      abundance = sum(abu),
      prevalence = sum(abu > 0)) %>%
    dplyr::ungroup()
}

#' plots a relative abundance heatmap with samples as columns
#' @param taxa_summary  a tibble with the abundance and prevalence for each taxa
#' @param taxa_rank a string with the taxonomic rank for which the data
#'  is summarized
#' @param meta a tibble with the meta information, must contains a key column
#'  with the same value than taxa_summary
#' @param meta_key name of that key, by default uses the value "sample"
#' @param min_prop the minimum frequency required to group taxa into "Other"
#' @export
relative_abundance_heatmap <- function(taxa_summary,
  taxa_rank, meta, meta_key = "sample", min_prop = 0.01) {

  plot_data <- taxa_summary %>%
    dplyr::mutate(
      taxa = factor(.data[[taxa_rank]]),
      taxa = fct_explicit_na(taxa, "Unlabel")) %>%
    dplyr::group_by(sample) %>%
    dplyr::select(-prevalence) %>%
    dplyr::mutate(rel = abundance / sum(abundance)) %>%
    dplyr::ungroup()

  plot_data %<>%
    dplyr::mutate(
      taxa = fct_lump_prop(taxa, min_prop, w = rel),
      taxa = factor(taxa))

  ncolors <- plot_data %>%
    dplyr::pull(taxa) %>%
    nlevels()

  if (ncolors >= 26) {
    pal1 <- Polychrome::green.armytage.colors(26)
    pal2 <- Polychrome::dark.colors(min(24, ncolors - 26))
    pal <- c(pal1, pal2)
    pal <- grDevices::colorRampPalette(pal)
    pal <- pal(ncolors)
  }else{
   pal <- Polychrome::green.armytage.colors(ncolors)
  }

  names(pal) <- plot_data %>%
    dplyr::pull(taxa) %>%
    levels()
  pal["Other"] <- "darkgrey"
  pal["Unlabel"] <- "lightblue"

  plot_data %<>%
    dplyr::left_join(meta, by = c(sample = meta_key))

  plot_data %>%
    dplyr::filter(!is.na(farm)) %>%
    ggplot2::ggplot(aes(sample, rel, fill = taxa)) +
    ggplot2::geom_col(width = .8) +
    ggplot2::facet_grid(. ~ interaction(farm, location), 
      scales = "free_x", space = "free_x") +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(1), expand = c(0, 0)) +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 8))

}
