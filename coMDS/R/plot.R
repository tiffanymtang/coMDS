#' Plotting function for coMDS and locoMDS objects
#'
#' @param obj A `coMDS` object, `locoMDS` object, or list of `locoMDS` objects
#'   (e.g., output of [locoMDS()] with multiple hyperparameter combinations).
#' @param color Optional vector of colors for points in the plot. Only used if
#'   type is "scores". If `NULL`, no color will be applied.
#' @param type Character string indicating the type of visualization. Can be
#'   either "scores" or "relative_errors" to show the consensus component scores
#'   or the relative errors from each input source, respectively. Default is
#'   "scores".
#' @param ... Additional arguments passed to the ggplot2::geom_point() if
#'   type is "scores" or ggplot2::geom_bar() if type is "relative_errors".
#'
#' @returns A ggplot object visualizing the consensus components of the input
#'   object.
#'
#' @examples
#' data(iris)
#' # remove duplicates so that tSNE can run
#' iris <- dplyr::distinct(iris)
#' X <- iris[, 1:4]
#' species <- iris$Species
#'
#' # fit various dimension reduction methods
#' pca_scores <- prcomp(X, center = TRUE, scale = TRUE)$x
#' tsne_scores <- Rtsne::Rtsne(X, dims = 2, perplexity = 30, verbose = FALSE)$Y
#' umap_scores <- umap::umap(X, n_components = 2, verbose = FALSE)$layout
#' dr_list <- list(
#'   pca = pca_scores,
#'   tsne = tsne_scores,
#'   umap = umap_scores
#' )
#'
#' # fit coMDS using dimension reduction embeddings directly as input
#' comds_out <- coMDS(embed_list = dr_list, ndim = 2)
#' # plot coMDS scores
#' plot_coMDS(comds_out, color = species, type = "scores")
#' # plot coMDS relative errors
#' plot_coMDS(comds_out, type = "relative_errors")
#'
#' # fit LoCoMDS with specific hyperparameters
#' locomds_out <- locoMDS(
#'   embed_list = dr_list, ndim = 2, tau = 0.1, percentile = 0.5
#' )
#' # plot LoCoMDS scores
#' plot_coMDS(locomds_out, color = species, type = "scores")
#' # plot LoCoMDS relative errors
#' plot_coMDS(locomds_out, type = "relative_errors")
#'
#' # fit LoCoMDS with multiple possible hyperparameters
#' locomds_multi_out <- locoMDS(
#'   embed_list = dr_list, ndim = 2, tau = c(0.01, 0.1), percentile = c(0.5, 0.8)
#' )
#' # plot LoCoMDS scores for multiple hyperparameters
#' plot_coMDS(locomds_multi_out, color = species, type = "scores")
#' # plot LoCoMDS relative errors for multiple hyperparameters
#' plot_coMDS(locomds_multi_out, type = "relative_errors")
#'
#' @export
plot_coMDS <- function(obj, color = NULL,
                       type = c("scores", "relative_errors"),
                       ...) {
  type <- match.arg(type)

  # check input
  input_mode <- "single"
  if (!inherits(obj, "coMDS") && !inherits(obj, "locoMDS")) {
    input_mode <- "multi"
    if (!isTRUE(all(sapply(obj, function(x) inherits(x, "locoMDS"))))) {
      stop("Input `obj` should either be a `coMDS` or `locoMDS` object.")
    }
  }
  if (type == "scores") {
    # create plotting data frame
    if (input_mode == "single") {
      plt_df <- as.data.frame(obj$gspace) |>
        setNames(paste0("Consensus Component ", 1:ncol(obj$gspace))) |>
        dplyr::mutate(.color = color)
    } else {
      plt_df <- purrr::map(
        obj,
        ~ as.data.frame(.x$gspace) |>
          setNames(paste0("Consensus Component ", 1:ncol(.x$gspace))) |>
          dplyr::mutate(
            # rescale to be between 0 and 1
            dplyr::across(
              tidyselect::everything(),
              ~ (.x - min(.x)) / (max(.x) - min(.x))
            ),
            tau = .x$tau,
            percentile = .x$percentile,
            .color = color
          )
      ) |>
        dplyr::bind_rows()
    }
    # make plot
    if (!is.null(color)) {
      plt <- ggplot2::ggplot(plt_df) +
        ggplot2::aes(
          x = `Consensus Component 1`,
          y = `Consensus Component 2`,
          color = .color
        )
    } else {
      plt <- ggplot2::ggplot(plt_df) +
        ggplot2::aes(
          x = `Consensus Component 1`,
          y = `Consensus Component 2`
        )
    }
    plt <- plt +
      ggplot2::geom_point(...) +
      ggplot2::labs(
        x = "Consensus Component 1",
        y = "Consensus Component 2",
        color = ""
      ) +
      ggplot2::theme_minimal()
  } else if (type == "relative_errors") {
    if (input_mode == "single") {
      plt_df <- as.data.frame(obj$sps_norm) |>
        setNames("Relative Error") |>
        tibble::rownames_to_column("Source") |>
        dplyr::arrange(`Relative Error`) |>
        dplyr::mutate(
          start = 0,
          Source = forcats::fct_inorder(Source)
        )
    } else {
      plt_df <- purrr::map(
        obj,
        ~ as.data.frame(.x$sps_norm) |>
          setNames("Relative Error") |>
          tibble::rownames_to_column("Source") |>
          dplyr::arrange(`Relative Error`) |>
          dplyr::mutate(
            start = 0,
            tau = .x$tau,
            percentile = .x$percentile,
            Source = forcats::fct_inorder(Source)
          )
      ) |>
        dplyr::bind_rows()
    }
    plt <- plt_df |>
      ggplot2::ggplot() +
      ggplot2::geom_bar(
        ggplot2::aes(
          x = Source,
          y = `Relative Error`
        ),
        fill = "#6FBBE3",
        stat = "identity",
        ...
      ) +
      ggplot2::geom_hline(
        yintercept = 100 / length(unique(plt_df$Source)),
        linetype = "dashed",
        color = "black"
      ) +
      ggplot2::theme_minimal()
  }
  if (input_mode == "multi") {
    plt <- plt +
      ggplot2::facet_grid(tau ~ percentile, labeller = ggplot2::label_both)
  }
  return(plt)
}

