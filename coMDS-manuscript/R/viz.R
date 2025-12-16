#' Plot 2D dimension reduction visualization
plot_dr <- function(df, color = NULL, point_size = 0.5, blank_axes = FALSE) {
  plt_df <- as.data.frame(df) |>
    dplyr::mutate(.color = color)
  if (!is.null(color)) {
    plt <- ggplot2::ggplot(plt_df) +
      ggplot2::aes(x = `Component 1`, y = `Component 2`, color = .color) +
      ggplot2::labs(color = "")
  } else {
    plt <- ggplot2::ggplot(plt_df) +
      ggplot2::aes(x = `Component 1`, y = `Component 2`)
  }
  plt <- plt +
    ggplot2::geom_point(
      size = point_size
    ) +
    ggplot2::theme_minimal()
  if (blank_axes) {
    plt <- plt +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank()
      )
  }
  return(plt)
}


#' Plot full distribution of evaluation metrics
plot_evaluation_distribution <- function(data,
                                         metrics = c("triplet", "spearman", "rf")) {
  plt_df <- data |>
    dplyr::select(.method_name, tidyselect::all_of(paste0(metrics, "_raw"))) |>
    dplyr::rename_with(~ stringr::str_remove(.x, "_raw")) |>
    dplyr::mutate(
      .color = dplyr::case_when(
        .method_name %in% c("CoMDS", "LoCoMDS") ~ "CoMDS/LoCoMDS",
        stringr::str_detect(.method_name, "Meta-Spec") ~ "Meta-Spec",
        stringr::str_detect(.method_name, "Multi-SNE") ~ "Multi-SNE",
        TRUE ~ "Other"
      ),
      .method_name = stringr::str_replace(.method_name, " \\(", "\n\\(") |>
        stringr::str_replace("\\# neighbors", "\\#nbhrs") |>
        stringr::str_replace("perplexity", "perp.")
    )
  if ("rf" %in% metrics) {
    plt_df <- plt_df |>
      tidyr::unnest_wider("rf", names_sep = "_")
  }
  plt_df <- plt_df |>
    tidyr::pivot_longer(
      cols = -c(.method_name, .color),
      names_to = "metric",
      values_to = "values"
    ) |>
    tidyr::unnest_longer("values") |>
    dplyr::mutate(
      metric = factor(
        metric,
        levels = c(
          "triplet", "spearman",
          "rf_pred_accuracy", "rf_macro_roc", "rf_multi_roc",
          "rf_match_prob", "rf_diff_prob"
        )
      ),
      metric_name = as.character(metric)
    )
  plt_ls <- plt_df |>
    dplyr::group_by(metric) |>
    dplyr::group_map(
      function(.x, .y) {
        metric_name <- .x$metric_name[[1]]
        if (metric_name %in% c("rf_diff_prob")) {
          group_plt <- ggplot2::ggplot(.x) +
            ggplot2::aes(
              # order in decreasing median value
              x = reorder(.method_name, -values, FUN = median),
              y = values,
              fill = .color
            )
        } else {
          group_plt <- ggplot2::ggplot(.x) +
            ggplot2::aes(
              # order in increasing median value
              x = reorder(.method_name, values, FUN = median),
              y = values,
              fill = .color
            )
        }
        group_plt <- group_plt +
          ggplot2::geom_boxplot(outlier.size = 0.5) +
          ggplot2::scale_fill_manual(
            values = c("#00BFC4", "#F8766D", "#c77cff", "white")
          ) +
          ggplot2::labs(x = "Method", y = metric_name, fill = "Method Type") +
          vthemes::theme_vmodern()
      }
    )
  plt <- patchwork::wrap_plots(plt_ls, ncol = 1, guides = "collect") +
    patchwork::plot_layout(axis_titles = "collect")
  return(plt)
}


#' Plot mean evaluation metrics
plot_evaluation_means <- function(data,
                                  metrics = c("triplet", "spearman", "rf", "lcmc_local")) {
  plt_df <- data |>
    dplyr::select(.method_name, tidyselect::all_of(metrics)) |>
    dplyr::mutate(
      .color = dplyr::case_when(
        .method_name %in% c("CoMDS", "LoCoMDS") ~ "CoMDS/LoCoMDS",
        stringr::str_detect(.method_name, "Meta-Spec") ~ "Meta-Spec",
        stringr::str_detect(.method_name, "Multi-SNE") ~ "Multi-SNE",
        TRUE ~ "Other"
      ),
      .method_name = stringr::str_replace(.method_name, " \\(", "\n\\(") |>
        stringr::str_replace("\\# neighbors", "\\#nbhrs") |>
        stringr::str_replace("perplexity", "perp.")
    )
  if ("rf" %in% metrics) {
    plt_df <- plt_df |>
      tidyr::unnest_wider("rf", names_sep = "_")
  }
  plt_df <- plt_df |>
    tidyr::pivot_longer(
      cols = -c(.method_name, .color),
      names_to = "metric",
      values_to = "values"
    ) |>
    dplyr::mutate(
      metric = factor(
        metric,
        levels = c(
          "triplet", "spearman", "lcmc_local",
          "rf_pred_accuracy", "rf_macro_roc", "rf_multi_roc",
          "rf_match_prob", "rf_diff_prob"
        )
      ),
      metric_name = as.character(metric)
    )
  plt_ls <- plt_df |>
    dplyr::group_by(metric) |>
    dplyr::group_map(
      function(.x, .y) {
        metric_name <- .x$metric_name[[1]]
        if (metric_name %in% c("rf_diff_prob")) {
          group_plt <- ggplot2::ggplot(.x) +
            ggplot2::aes(
              # order in decreasing median value
              x = reorder(.method_name, -values, FUN = median),
              y = values,
              fill = .color
            )
        } else {
          group_plt <- ggplot2::ggplot(.x) +
            ggplot2::aes(
              # order in increasing median value
              x = reorder(.method_name, values, FUN = median),
              y = values,
              fill = .color
            )
        }
        group_plt <- group_plt +
          ggplot2::geom_bar(color = "black", stat = "identity") +
          ggplot2::scale_fill_manual(
            values = c("#00BFC4", "#F8766D", "#c77cff", "white")
          ) +
          ggplot2::labs(x = "Method", y = metric_name, fill = "Method Type") +
          vthemes::theme_vmodern()
      }
    )
  plt <- patchwork::wrap_plots(plt_ls, ncol = 1, guides = "collect") +
    patchwork::plot_layout(axis_titles = "collect")
  return(plt)
}


#' Plot LCMC curves
plot_lcmc_curves <- function(data) {
  plt_df <- data |>
    dplyr::select(.method_name, lcmc_local_raw) |>
    dplyr::mutate(
      .color = dplyr::case_when(
        .method_name %in% c("CoMDS", "LoCoMDS") ~ .method_name,
        stringr::str_detect(.method_name, "Meta-Spec") ~ .method_name,
        stringr::str_detect(.method_name, "Multi-SNE") ~ "Multi-SNE",
        TRUE ~ "Other"
      ),
      .method_name = stringr::str_replace(
        .method_name, "\\# neighbors", "\\# nbhrs"
      ) |>
        stringr::str_replace("perplexity", "perp.")
    ) |>
    dplyr::rename_with(~ stringr::str_remove(.x, "_raw")) |>
    tidyr::pivot_longer(
      cols = -c(.method_name, .color),
      names_to = "type",
      values_to = "lcmc"
    ) |>
    tidyr::unnest(lcmc) |>
    dplyr::mutate(
      type = R.utils::capitalize(stringr::str_remove(type, "lcmc_")) |>
        factor(levels = c("Local"))
    )
  plt <- ggplot2::ggplot(plt_df) +
    ggplot2::aes(
      x = k, y = lcmc, group = .method_name, color = .color
    ) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::facet_grid(~ type, scales = "free_x") +
    ggplot2::labs(
      x = "Neighborhood size (k value)",
      y = "Average LCMC",
      color = "Method"
    ) +
    vthemes::theme_vmodern(size_preset = "medium")
  return(plt)
}


#' Plot stability experiment results
plot_stability <- function(data, method = "pearson"){
  upper_fn <- function(data, mapping, ...) {
    x <- as.numeric(data[[rlang::as_name(mapping$x)]])
    y <- as.numeric(data[[rlang::as_name(mapping$y)]])
    corr <- stats::cor(x, y, method = method, use = "complete.obs")
  
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5,
                        label = paste0("correlation = \n", round(corr, 3)),
                        size = 4) +
      ggplot2::theme_void()
  }

  plt <- GGally::ggpairs(
    data,
    columns = which(!names(data) %in% "label"),
    mapping = ggplot2::aes(color = label),
    lower = list(
      continuous = GGally::wrap(
        "points",
        size = 0.25
      )
    ),
    diag  = list(
      continuous = GGally::wrap(
        "densityDiag",
        size = 0.4
      )
    ),
    upper = list(
      continuous = upper_fn  # custom upper function
    )
  )+
  ggplot2::theme_bw()
  return(plt)
}

#' Plot outlier points in 2D embeddings
plot_outlier <- function(df, point_size = 0.5, transparency = 0.2, blank_axes = FALSE, 
                label_size = NULL, max_overlaps = NULL, color = NULL){
  df$ID <- 1:nrow(df)
  plt_df <- as.data.frame(df)
  plt_df$.outlier <- plt_df$outlier
  plt_df$.color <- color

  if (is.null(label_size)) {
    label_size <- point_size * 2
  }

  if (is.null(max_overlaps)) {
    max_overlaps <- sum(plt_df$.outlier == 1)*3
  }

  this_plot <- ggplot2::ggplot(plt_df, ggplot2::aes(x = `Component 1`, y = `Component 2`)) +
    ggplot2::geom_point(
      ggplot2::aes_string(color = ".color", alpha = ".outlier"),
      size = point_size
    ) +
    ggplot2::scale_alpha_continuous(range = c(transparency, 1)) +
    ggrepel::geom_text_repel(
      data = plt_df[plt_df$.outlier == 1, ],
      ggplot2::aes_string(color = ".color", label = "ID"),
      size = label_size + 1,
      max.overlaps = max_overlaps,
      show.legend = FALSE
    ) +
    ggplot2::theme_minimal() +
    ggplot2::guides(alpha = "none")

  if (blank_axes) {
    this_plot <- this_plot +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )
  }
  return(this_plot)
}

