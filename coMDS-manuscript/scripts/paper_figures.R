rm(list = ls())

set.seed(331)

data_dir <- here::here("data")
results_dir <- here::here("results")
figures_dir <- here::here("results", "figures")

datasets <- list(
  "HIV" = "hiv",
  "4EQ" = "4eq",
  "8EQ" = "8eq",
  "Trajectory" = "traj",
  "Cycle" = "cycle",
  "Olive Oil" = "olive",
  "Wheat" = "wheat",
  "Wholesale" = "wholesale",
  "Star" = "star_rs",
  "Gaussian Simulation" = "gaussian",
  "Swiss Roll Simulation" = "swiss_roll"
)
method_colors <- c(
  "CoMDS" = "#295B69",
  "LoCoMDS" = "#62BCD3",
  "Meta-Spec (UMAP)" = "#E18235",
  "Meta-Spec (kPCA)" = "#FFB45F",
  "Multi-SNE" = "#84706C",
  "Multi-SNE (perplexity = 30)" = "#84706C"
)
method_names <- c(
  "PCA",
  "MDS",
  "Non-metric MDS",
  "Sammon",
  "kPCA (sigma = 0.01)",
  "kPCA (sigma = 0.001)",
  "LLE",
  "HLLE",
  "Isomap",
  "LEIM",
  "PHATE (# neighbors = 30)",
  "PHATE (# neighbors = 50)",
  "tSNE (perplexity = 30)",
  "tSNE (perplexity = 50)",
  "UMAP (# neighbors = 30)",
  "UMAP (# neighbors = 50)",
  "CoMDS",
  "LoCoMDS",
  "Meta-Spec (kPCA)",
  "Meta-Spec (UMAP)",
  "Multi-SNE (perplexity = 30)"
)
metric_names <- c(
  triplet = "Random Triplets",
  spearman = "Spearman Correlation",
  lcmc_local = "Local LCMC",
  multi_roc = "RF AUROC (Multi-class)",
  diff_prob = "RF Probability\nof Wrong Class"
)

cluster_colors <- c(
  "#25B0D6", "#E9BD52", "#DA6DA9", "#85B969", "#948EC9", 
  "#4E443A", "#EE8D2B", "#668189", "#984848"
)
eq_cluster_colors <- c(
  "#25B0D6", "#E9BD52", "#948EC9", "#EE8D2B",
  "#4E443A", "#DA6DA9", "#668189", "#85B969"
)
ordered_cluster_colors <- c(
  "#DA6DA9", "#E9BD52", "#85B969", "#25B0D6", "#948EC9"
)

for (fname in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) {
  source(fname)
}

clean_method_names <- function(x, multiline = FALSE) {
  if (stringr::str_detect(x, "kPCA") & stringr::str_detect(x, "sigma")) {
    sigma_val <- stringr::str_extract(x, "(?<=sigma = )\\d+\\.\\d+")
    if (multiline) {
      x <- bquote(bold(atop("kPCA", (sigma == .(sigma_val)))))
    } else {
      x <- bquote(bold(kPCA~(sigma == .(sigma_val))))
    }
  } else if (stringr::str_detect(x, "neighbors")) {
    n_neighbors_val <- stringr::str_extract(x, "(?<=# neighbors = )\\d+")
    method_name <- stringr::str_extract(x, "^[^ ]+")
    if (multiline) {
      x <- bquote(bold(atop(.(method_name), (n[neighbors] == .(n_neighbors_val)))))
    } else {
      x <- bquote(bold(.(method_name)~(n[neighbors] == .(n_neighbors_val))))
    }
  } else if (stringr::str_detect(x, "tSNE") && !stringr::str_detect(x, "Multi")) {
    x <- stringr::str_replace(x, "tSNE", "t-SNE")
    if (multiline) {
      x <- stringr::str_replace(x, "t-SNE ", "t-SNE\n")
    }
  } else if (stringr::str_detect(x, "Multi-SNE")) {
    if (stringr::str_detect(x, "30")) {
      x <- "Multi-SNE"
    }
  }
  return(x)
}

add_dr_plot_theme <- function(plt, method_name, title_size = 14, legend_text_size = 12, 
                              legend_point_size = 3, grid_width = 0.25) {
  plt <- plt +
    ggplot2::labs(
      title = method_name
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, 1, by = grid_width)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, 1, by = grid_width)
    ) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = legend_text_size),
      plot.title = ggplot2::element_text(hjust = 0.5, size = title_size, face = "bold")
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = legend_point_size)))
  return(plt)
}

# helper ggplot themes and objects
tag_plot_theme <- ggplot2::theme(
  plot.tag = ggplot2::element_text(
    face = "italic", size = 14, vjust = 1,
    margin = ggplot2::margin(r = 10, l = -10)
  ),
  plot.tag.position = c(0, 1)
)
divider <- ggplot2::ggplot() +
  ggplot2::theme_void() +
  ggplot2::geom_hline(yintercept = 0, color = "grey", linewidth = 0.5)
vdivider <- ggplot2::ggplot() +
  ggplot2::theme_void() +
  ggplot2::geom_vline(xintercept = 0, color = "grey", linewidth = 0.5)

# create figures directory if it doesn't exist
if (!dir.exists(figures_dir)) {
  dir.create(figures_dir, recursive = TRUE)
}

#### Load in data ####
data_ls <- list()
rm_datasets <- c("Gaussian Simulation", "Swiss Roll Simulation")
for (data_name in setdiff(names(datasets), rm_datasets)) {
  load(file.path(data_dir, paste0(datasets[[data_name]], ".RData")))
  if (data_name == "HIV") {
    info <- dplyr::case_when(
      info == "plasmablast" ~ "Plasmablast",
      info == "monocyte" ~ "Monocyte",
      TRUE ~ stringr::str_replace(info, "cell", "Cell")
    )
  } else if (data_name %in% c("4EQ", "8EQ")) {
    info <- dplyr::case_when(
      info == "b.cells" ~ "B Cells",
      info == "cd14.monocytes" ~ "CD14 Monocytes",
      info == "naive.cytotoxic" ~ "Naive Cytotoxic",
      info == "regulatory.t" ~ "Regulatory T",
      info == "cd4.t.helper" ~ "CD4 T Helper",
      info == "cd56.nk" ~ "CD56 NK",
      info == "memory.t" ~ "Memory T",
      info == "naive.t" ~ "Naive T"
    )
  } else if (data_name == "Star") {
    info <- dplyr::case_when(
      info == "GALAXY" ~ "Galaxy",
      info == "QSO" ~ "Quasar",
      info == "STAR" ~ "Star"
    )
  }
  data_ls[[data_name]] <- list(
    "X" = data,
    "info" = info
  )
}

data_ls[["Gaussian Simulation"]] <- readRDS(file.path(results_dir, "gaussian", "Gaussian.rds"))
data_ls[["Gaussian Simulation"]] <- setNames(data_ls[["Gaussian Simulation"]], c("X", "info"))

data_ls[["Swiss Roll Simulation"]] <- readRDS(file.path(results_dir, "swiss_roll", "Swiss_roll.rds"))
data_ls[["Swiss Roll Simulation"]] <- setNames(data_ls[["Swiss Roll Simulation"]], c("X", "info"))

#### Plot Original 3D data for simulation datasets ####
# Gaussian Simulation
info_vals <- data_ls[["Gaussian Simulation"]]$info
clusters <- sort(unique(info_vals), decreasing = FALSE)
point_colors <- cluster_colors[info_vals]
pdf(file = file.path(figures_dir, "Gaussian3d.pdf"), width = 8, height = 6)
par(mai = c(0.05, 0.05, 0.05, 0.05), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
scatterplot3d::scatterplot3d(
  data_ls[["Gaussian Simulation"]]$X,
  pch = 16,
  color = point_colors,
  grid = TRUE,
  box = FALSE,
  scale.y = 1,
  angle = 240,
  xlab = '', ylab = '', zlab = '',
  x.ticklabs = rep("", 10),
  y.ticklabs = rep("", 10),
  z.ticklabs = rep("", 10)
)
graphics::par(fig = c(0, 1, 0.88, 1), new = TRUE, mar = c(0,0,0,0))
plot.new()
legend(
  "center", 
  legend = clusters,
  col = cluster_colors[clusters],
  pch = 16,
  horiz = TRUE,
  bty = "n",
  cex = 2,
  pt.cex = 3 
)
dev.off()

# Swiss Roll Simulation
info_vals <- as.numeric(data_ls[["Swiss Roll Simulation"]]$info)
grad_fun <- scales::gradient_n_pal(ordered_cluster_colors)
color_vals <- grad_fun(scales::rescale(info_vals))
pdf(file = file.path(figures_dir, "SwissRoll3d_legend.pdf"), width = 8, height = 6)
par(mai = c(0.1, 0.05, 0.05, 0.05), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
s3d <- scatterplot3d::scatterplot3d(
  data_ls[["Swiss Roll Simulation"]]$X[, c(2, 3, 1)],
  pch = 16,
  color = color_vals,
  grid = TRUE,
  box = FALSE,
  scale.y = 1,
  angle = 105,
  xlab = '', ylab = '', zlab = '',
  x.ticklabs = rep("", 10),
  y.ticklabs = rep("", 10),
  z.ticklabs = rep("", 10)
)
graphics::par(fig = c(0, 1, 0.88, 1), new = TRUE, mar = c(0,0,0,0))
  fields::image.plot(
  legend.only = TRUE,
  zlim = range(info_vals),
  col = grad_fun(seq(0,1,length.out = 100)),
  horizontal = TRUE,
  legend.lab = "",
  axis.args = list(cex.axis = 2),
  smallplot = c(0.1, 0.9, 0.5, 0.99) #legend position
)
dev.off()


#### Get dimension reduction embedding results ####
dr_list_all <- list()
locomds_list_all <- list()
for (data_name in names(datasets)) {
  dr_out <- readRDS(
    file.path(results_dir, datasets[[data_name]], "dimension_reduction_methods.rds")
  )
  comds_out <- readRDS(
    file.path(results_dir, datasets[[data_name]], "comds.rds")
  )
  locomds_out <- readRDS(
    file.path(results_dir, datasets[[data_name]], "locomds.rds")
  )
  metaspec_out <- readRDS(
    file.path(results_dir, datasets[[data_name]], "metaspec.rds")
  )
  multitsne_out <- readRDS(
    file.path(results_dir, datasets[[data_name]], "multisne.rds")
  )
  locomds_list_all[[data_name]] <- locomds_out
  dr_list_all[[data_name]] <- c(
    dr_out,
    list(
      "CoMDS" = comds_out$embeddings,
      "LoCoMDS" = locomds_out$embeddings,
      "Meta-Spec (kPCA)" = metaspec_out$embeddings$kPCA,
      "Meta-Spec (UMAP)" = metaspec_out$embeddings$UMAP,
      "Multi-SNE (perplexity = 30)" = multitsne_out$embeddings$perplexity_30
    )
  )[method_names] |>
    purrr::map(
      function(.x) {
        # normalize to [0, 1] range for plotting
        for (j in 1:ncol(.x)) {
          .x[, j] <- (.x[, j] - min(.x[, j])) / (max(.x[, j]) - min(.x[, j]))
        }
        return(as.data.frame(.x))
      }
    )
}

#### Reorder samples for certain datasets for clearer plotting ####
reorder_idx <- purrr::map(
  c("B Cells", "Naive Cytotoxic", "Regulatory T", "CD14 Monocytes"),
  ~ which(data_ls[["4EQ"]]$info == .x)
) |>
  unlist()
data_ls[["4EQ"]]$X <- data_ls[["4EQ"]]$X[reorder_idx, ]
data_ls[["4EQ"]]$info <- data_ls[["4EQ"]]$info[reorder_idx]
dr_list_all[["4EQ"]] <- purrr::map(
  dr_list_all[["4EQ"]], ~ .x[reorder_idx, ]
)
locomds_list_all[["4EQ"]]$fit <- purrr::map(
  locomds_list_all[["4EQ"]]$fit,
  function(fit) {
    fit$gspace <- fit$gspace[reorder_idx, ]
    return(fit)
  }
)

#### Plot all dimension reduction methods ####
point_size <- 0.5
keep_methods <- method_names
for (data_name in names(datasets)) {
  if (data_name %in% c("Trajectory", "Swiss Roll Simulation")) {
    this_cluster_colors <- ordered_cluster_colors
  } else if (data_name == "8EQ") {
    this_cluster_colors <- eq_cluster_colors
  } else {
    this_cluster_colors <- cluster_colors
  }

  if (data_name != "Swiss Roll Simulation") {
    plt_ls <- purrr::imap(
      dr_list_all[[data_name]][keep_methods],
      function(plt_df, method_name) {
        plt <- plot_dr(
          plt_df,
          color = as.factor(data_ls[[data_name]]$info),
          point_size = point_size
        ) |>
          add_dr_plot_theme(
            method_name = clean_method_names(method_name),
          ) +
          ggplot2::scale_color_manual(
            values = this_cluster_colors
          )
        return(plt)
      }
    )
  } else {
    plt_ls <- purrr::imap(
      dr_list_all[[data_name]][keep_methods],
      function(plt_df, method_name) {
        color_vals <- as.numeric(data_ls[[data_name]]$info)
        plt <- plot_dr(
          plt_df,
          color = color_vals,   # continuous!
          point_size = point_size
        ) |>
          add_dr_plot_theme(
            method_name = clean_method_names(method_name)
          ) +
          ggplot2::scale_color_gradientn(
            colours = this_cluster_colors    # interpolates through your 5 colors
          ) +
          ggplot2::guides(color = ggplot2::guide_colorbar())

        return(plt)
      }  
    )
  }
  # add tag to first plot
  plt_ls[[1]] <- plt_ls[[1]] +
    ggplot2::labs(tag = "(A)") +
    tag_plot_theme
  plt_ls[[17]] <- plt_ls[[17]] +
    ggplot2::labs(tag = "(B)") +
    tag_plot_theme
  plt_ls <- c(
    plt_ls[1:16],
    divider, # add horizontal line between non-consensus and consensus methods
    plt_ls[17:length(plt_ls)]
  )
  design <- "##AAAABBBBCCCCDDDD##
             ##EEEEFFFFGGGGHHHH##
             ##IIIIJJJJKKKKLLLL##
             ##MMMMNNNNOOOOPPPP##
             QQQQQQQQQQQQQQQQQQQQ
             RRRRSSSSTTTTUUUUVVVV"
  plt <- patchwork::wrap_plots(
    plt_ls, byrow = TRUE, guides = "collect",
    heights = c(1, 1, 1, 1, 0.1, 1),
    design = design
  ) &
    ggplot2::theme(
      plot.margin = ggplot2::margin(l = 10, r = 0, t = 5, b = 5)
    )
  ggplot2::ggsave(
    plt,
    filename = file.path(
      figures_dir,
      sprintf("full_dimension_reduction_plots_%s.pdf", datasets[[data_name]])
    ),
    width = 13, height = 12.5
  )
}


#### Figure 1 ####
keep_methods <- c("PCA", "kPCA (sigma = 0.01)", "tSNE (perplexity = 30)", "UMAP (# neighbors = 30)")
point_size <- 0.75
title_size <- 15
for (data_name in c("HIV")) {
  if (data_name %in% c("Trajectory", "Swiss Roll Simulation")) {
    this_cluster_colors <- ordered_cluster_colors
  } else if (data_name == "8EQ") {
    this_cluster_colors <- eq_cluster_colors
  } else {
    this_cluster_colors <- cluster_colors
  }

  if (data_name != "Swiss Roll Simulation") {
    plt_ls <- purrr::imap(
      dr_list_all[[data_name]][keep_methods],
      function(plt_df, method_name) {
        plt <- plot_dr(
          plt_df,
          color = as.factor(data_ls[[data_name]]$info),
          point_size = point_size
        ) |>
          add_dr_plot_theme(
            method_name = clean_method_names(method_name, multiline = TRUE),
            title_size = title_size
          ) +
          ggplot2::scale_color_manual(
            values = this_cluster_colors
          )
        return(plt)
      }
    )
  } else {
    plt_ls <- purrr::imap(
      dr_list_all[[data_name]][keep_methods],
      function(plt_df, method_name) {
        color_vals <- as.numeric(data_ls[[data_name]]$info)
        plt <- plot_dr(
          plt_df,
          color = color_vals,   # continuous!
          point_size = point_size
        ) |>
          add_dr_plot_theme(
            method_name = clean_method_names(method_name, multiline = TRUE),
            title_size = title_size
          ) +
          ggplot2::scale_color_gradientn(
            colours = this_cluster_colors
          ) +
          ggplot2::guides(color = ggplot2::guide_colorbar())
        
        return(plt)
      }  
    )
  }
  plt <- patchwork::wrap_plots(
    plt_ls, nrow = 1, guides = "collect"
  )
  ggplot2::ggsave(
    plt,
    filename = file.path(
      figures_dir,
      paste0("figure1_", datasets[[data_name]], ".pdf")
    ),
    width = 10.5, height = 3
  )
}


#### Main Text DR Scatter Plot Figure ####
keep_methods <- c(
  "PCA", "kPCA (sigma = 0.01)", "HLLE", "tSNE (perplexity = 30)", "UMAP (# neighbors = 30)", 
  "CoMDS", "LoCoMDS", "Meta-Spec (UMAP)", "Meta-Spec (kPCA)", "Multi-SNE (perplexity = 30)"
)
title_size <- 15
# point_size <- 0.42
point_size <- 0.75
for (data_name in names(datasets)) {
  if (data_name %in% c("Trajectory", "Swiss Roll Simulation")) {
    this_cluster_colors <- ordered_cluster_colors
  } else if (data_name == "8EQ") {
    this_cluster_colors <- eq_cluster_colors
  } else {
    this_cluster_colors <- cluster_colors
  }

  if (data_name != "Swiss Roll Simulation") {
    plt_ls <- purrr::imap(
      dr_list_all[[data_name]][keep_methods],
      function(plt_df, method_name) {
        plt <- plot_dr(
          plt_df,
          color = as.factor(data_ls[[data_name]]$info),
          point_size = point_size
        ) |>
          add_dr_plot_theme(
            method_name = clean_method_names(method_name, multiline = TRUE),
            title_size = title_size
          ) +
          ggplot2::scale_color_manual(
            values = this_cluster_colors
          ) +
          ggplot2::theme(
            legend.position = "bottom",
            legend.text = ggplot2::element_text(margin = ggplot2::margin(l = 4, r = 10))
          )
        return(plt)
      }
    )
    plt <- patchwork::wrap_plots(
      plt_ls, ncol = 5, byrow = TRUE, guides = "collect"
    ) &
      ggplot2::theme(
        legend.position = "bottom"
      )
  } else {
    plt_ls <- purrr::imap(
      dr_list_all[[data_name]][keep_methods],
      function(plt_df, method_name) {
        color_vals <- as.numeric(data_ls[[data_name]]$info)
        plt <- plot_dr(
          plt_df,
          color = color_vals,   # continuous!
          point_size = point_size
        ) |>
          add_dr_plot_theme(
            method_name = clean_method_names(method_name, multiline = TRUE),
            title_size = title_size
          ) +
          ggplot2::scale_color_gradientn(
            colours = this_cluster_colors
          ) +
          ggplot2::guides(color = ggplot2::guide_colorbar()) +
          ggplot2::theme(
            legend.position = "none",
            legend.text = ggplot2::element_text(margin = ggplot2::margin(l = 4, r = 10))
          )
        
        return(plt)
      }  
    )
    plt <- patchwork::wrap_plots(
      plt_ls, ncol = 5, byrow = TRUE, guides = "collect"
    ) &
      ggplot2::theme(
        legend.position = "none"
      )
  }
  ggplot2::ggsave(
    plt,
    filename = file.path(
      figures_dir,
      sprintf("half_dimension_reduction_plots_%s.pdf", datasets[[data_name]])
    ),
    width = 12, height = 6
  )
}

#### Tuning appendix results ####
keep_percentiles <- c(0.2, 0.4, 0.6, 0.8)
keep_taus <- c(10, 1, 0.1, 0.01, 0.001)
point_size <- 0.5
title_size <- 16
rm_datasets <- c("Gaussian Simulation", "Swiss Roll Simulation")
simple_tuning_plt_ls <- list()
all_tuning_plt_ls <- list()
for (data_name in setdiff(names(datasets), rm_datasets)) {
  tuning_plt <- locomds_list_all[[data_name]]$tuning$plot +
    ggplot2::scale_shape_manual(
      values = rep(16, 81)
    ) +
    ggplot2::scale_color_manual(
      values = rep(scales::hue_pal()(9), each = 9)
    ) +
    ggplot2::labs(
      title = data_name
    ) +
    vthemes::theme_vmodern(size_preset = "medium") +
    ggplot2::theme(
      strip.text = ggplot2::element_blank()
    ) +
    ggplot2::guides(
      color = "none",
      shape = "none"
    )
  tuning_plt$data <- tuning_plt$data |>
    dplyr::filter(
      .mode == "All"
    ) |>
    dplyr::mutate(
      parameter_name = purrr::map_chr(
        parameter_name, 
        function(x) {
          parts <- stringr::str_split(x, "\\s*\\|\\s*", simplify = TRUE)
          return(stringr::str_c(parts[2], parts[1], sep = " | "))
        }
      )
    )
  all_tuning_plt_ls[[data_name]] <- tuning_plt

  simple_tuning_plt <- locomds_list_all[[data_name]]$tuning$plot +
    ggplot2::scale_shape_manual(
      values = rep(16, 9)
    ) +
    ggplot2::labs(
      title = data_name
    ) +
    vthemes::theme_vmodern(size_preset = "medium") +
    ggplot2::theme(
      strip.text = ggplot2::element_blank()
    )
  simple_tuning_plt$data <- simple_tuning_plt$data |>
    dplyr::filter(
      .mode == "Percentiles with best tau"
    ) |>
    dplyr::group_by(percentile) |>
    dplyr::mutate(
      is_smallest_tau = (tau == min(tau)),
      parameter_name = purrr::map_chr(
        parameter_name, 
        function(x) {
          parts <- stringr::str_split(x, "\\s*\\|\\s*", simplify = TRUE)
          return(parts[2])
        }
      )
    ) |>
    dplyr::filter(
      is_smallest_tau
    )
  simple_tuning_plt_ls[[data_name]] <- simple_tuning_plt

  if (data_name == "Trajectory") {
    this_cluster_colors <- ordered_cluster_colors
  } else if (data_name == "8EQ") {
    this_cluster_colors <- eq_cluster_colors
  } else {
    this_cluster_colors <- cluster_colors
  }
  plt_df <- purrr::map(
    locomds_list_all[[data_name]]$fit,
    function(locomds_fit) {
      if (!(locomds_fit$percentile %in% keep_percentiles &
            locomds_fit$tau %in% keep_taus)) {
        return(NULL)
      }
      plt_df <- as.data.frame(locomds_fit$gspace) |>
        setNames(c("Component 1", "Component 2")) |>
        dplyr::mutate(
          dplyr::across(
            c(`Component 1`, `Component 2`),
            ~ (.x - min(.x)) / (max(.x) - min(.x))
          ),
          Tau = locomds_fit$tau,
          Percentile = locomds_fit$percentile,
          .color = as.factor(data_ls[[data_name]]$info)
        )
      return(plt_df)
    }
  ) |>
    dplyr::bind_rows()
  plt <- plot_dr(
    plt_df,
    color = plt_df$.color,
    point_size = point_size
  ) +
    ggplot2::facet_grid(
      Percentile ~ Tau, 
      labeller = ggplot2::label_both
    )
  plt <- plt |>
    add_dr_plot_theme(
      method_name = data_name,
      title_size = title_size
    ) +
    ggplot2::scale_color_manual(
      values = this_cluster_colors
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(margin = ggplot2::margin(l = 4, r = 10)),
      strip.text = ggplot2::element_text(size = 14, face = "bold")
    )
  ggplot2::ggsave(
    plt,
    filename = file.path(
      figures_dir,
      sprintf("locomds_hyperparameter_%s.pdf", datasets[[data_name]])
    ),
    width = 12, height = 9
  )
}

all_tuning_plt <- patchwork::wrap_plots(
  all_tuning_plt_ls, ncol = 3, byrow = TRUE, guides = "collect"
) +
  patchwork::plot_layout(axis_titles = "collect") &
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "italic")
  )
ggplot2::ggsave(
  all_tuning_plt,
  filename = file.path(figures_dir, "locomds_hyperparameter_trace_all.pdf"),
  width = 10.25, height = 8
)

simple_tuning_plt <- patchwork::wrap_plots(
  simple_tuning_plt_ls, ncol = 3, byrow = TRUE, guides = "collect"
) +
  patchwork::plot_layout(axis_titles = "collect") &
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "italic")
  )
ggplot2::ggsave(
  simple_tuning_plt,
  filename = file.path(figures_dir, "locomds_hyperparameter_trace_simple.pdf"),
  width = 12, height = 8
)

#### Get evaluation results ####
keep_methods <- c(
  "CoMDS", "LoCoMDS", "Meta-Spec (UMAP)", "Meta-Spec (kPCA)", "Multi-SNE (perplexity = 30)"
)
rm_datasets <- c("Gaussian Simulation", "Swiss Roll Simulation")
eval_ls <- list()
for (data_name in setdiff(names(datasets), rm_datasets)) {
  eval_ls[[data_name]] <- readRDS(
    file.path(results_dir, datasets[[data_name]], "evaluation_results.rds")
  )
}
eval_tib <- dplyr::bind_rows(eval_ls, .id = ".dataset") |>
  dplyr::filter(
    !(.dataset %in% rm_datasets),
    .method_name %in% keep_methods
  ) |>
  dplyr::mutate(
    .method_name = ifelse(
      stringr::str_detect(.method_name, "Multi-SNE"), "Multi-SNE", .method_name
    )
  )

eval_summary_tib <- eval_tib |>
  dplyr::group_by(.dataset, .method_name) |>
  dplyr::summarise(
    dplyr::across(
      c(lcmc_local, triplet, spearman),
      mean
    ),
    .groups = "drop"
  )

#### Plot global/local preservation results ####
# flipped version
plt_ls <- eval_summary_tib |>
  tidyr::pivot_longer(
    cols = c(triplet, spearman, lcmc_local),
    names_to = "metric",
    values_to = "value"
  ) |>
  dplyr::mutate(
    .method_name = factor(.method_name, levels = rev(names(method_colors))),
    .dataset = factor(.dataset, levels = names(datasets)),
    metric = factor(metric, levels = names(metric_names))
  ) |>
  dplyr::group_by(metric) |>
  dplyr::group_map(
    function(.x, .y) {
      best_methods <- .x |>
        dplyr::group_by(.dataset) |>
        dplyr::mutate(
          is_best = ifelse(value == max(value), "**", ""),
          is_2nd_best = ifelse(value == sort(value, decreasing = TRUE)[2], "*", "")
        )
      plt <- ggplot2::ggplot(.x) +
        ggplot2::aes(
          x = .dataset, y = value, fill = .method_name
        ) +
        ggplot2::geom_bar(
          linewidth = 0.1,
          color = "grey98",
          stat = "identity", 
          position = ggplot2::position_dodge(width = 0.8),
          width = 0.8
        ) +
        ggplot2::geom_text(
          data = best_methods,
          ggplot2::aes(label = is_best, color = .method_name, group = .method_name),
          position = ggplot2::position_dodge(width = 0.82),
          hjust = -0.2,
          vjust = 0.8,
          size = 5
        ) +
        ggplot2::geom_text(
          data = best_methods,
          ggplot2::aes(label = is_2nd_best, color = .method_name, group = .method_name),
          position = ggplot2::position_dodge(width = 0.82),
          hjust = -0.2,
          vjust = 0.8,
          size = 5
        ) +
        ggplot2::scale_fill_manual(
          values = method_colors
        ) +
        ggplot2::scale_color_manual(
          values = method_colors
        ) +
        ggplot2::labs(
          x = "Dataset",
          y = metric_names[.y$metric[[1]]],
          fill = "Method"
        ) +
        ggplot2::scale_y_continuous(
          expand = ggplot2::expansion(mult = c(0, 0.1))
        ) +
        ggplot2::coord_flip() +
        vthemes::theme_vmodern(size_preset = "large") +
        ggplot2::theme(
          legend.position = "bottom",
          legend.text = ggplot2::element_text(margin = ggplot2::margin(l = 4, r = 10))
        ) +
        ggplot2::guides(
          fill = ggplot2::guide_legend(reverse = TRUE),
          color = "none"
        )
      if (.y$metric[[1]] != "triplet") {
        plt <- plt +
          ggplot2::theme(
            axis.title.y = ggplot2::element_blank()
          )
      }
      return(plt)
    }
  )
plt <- patchwork::wrap_plots(
  c(plt_ls[1:2], vdivider, plt_ls[3]),
  widths = c(1, 1, 0.1, 1),
  nrow = 1, guides = "collect"
) &
  ggplot2::theme(
    legend.position = "bottom"
  ) +
  patchwork::plot_layout(
    axis_titles = "collect"
  )
ggplot2::ggsave(
  plt,
  filename = file.path(figures_dir, "evaluation_structure_preservation_flipped.pdf"),
  width = 14, height = 6
)

# Calculate and save ranks
rank_df <- eval_summary_tib |>
  dplyr::group_by(.dataset) |>
  dplyr::mutate(
    dplyr::across(
      tidyselect::any_of(names(metric_names)),
      ~ rank(-.x)
    )
  )
for (metric_name in c("triplet", "spearman", "lcmc_local")) {
  rank_summary <- rank_df |>
    dplyr::rename(
      .metric := !!rlang::sym(metric_name)
    ) |>
    dplyr::group_by(.method_name) |>
    dplyr::summarise(
      `Rank 1` = sum(.metric == 1),
      `Rank 2` = sum(.metric == 2),
      `Rank 3` = sum(.metric == 3),
      `Rank 4` = sum(.metric == 4),
      `Rank 5` = sum(.metric == 5),
      `Average Rank` = mean(.metric),
      .groups = "drop"
    ) |>
    dplyr::rename(
      Method = .method_name
    )
  write.csv(
    rank_summary,
    file = file.path(figures_dir, sprintf("evaluation_structure_preservation_%s_ranks.csv", metric_name)),
    row.names = FALSE
  )
}

#### Plot RF results ####

# flipped version, with only one multi AUROC
plt_df_wide <- eval_tib |>
  dplyr::select(.dataset, .method_name, rf) |>
  dplyr::mutate(
    rf = purrr::map(rf, function(x) {
      if (is.data.frame(x)) {
        # convert 1×6 df -> numeric vector
        vec <- unlist(x[1, ])
        names(vec) <- colnames(x)
        vec
      } else if (is.numeric(x)) {
        # already a numeric vector
        x
      } else {
        # fallback
        unlist(x)
      }
    })
  ) |>
  tidyr::unnest_wider("rf") 

plt_ls <- plt_df_wide |>
  tidyr::pivot_longer(
    cols = c(multi_roc, diff_prob),
    names_to = "metric",
    values_to = "value"
  ) |>
  dplyr::mutate(
    .method_name = factor(.method_name, levels = rev(names(method_colors))),
    .dataset = factor(.dataset, levels = names(datasets)),
    metric = factor(metric, levels = names(metric_names))
  ) |>
  dplyr::group_by(metric) |>
  dplyr::group_map(
    function(.x, .y) {
      if (.y$metric[[1]] %in% c("diff_prob")) {
        best_methods <- .x |>
          dplyr::group_by(.dataset) |>
          dplyr::mutate(
            is_best = ifelse(value == min(value), "**", ""),
            is_2nd_best = ifelse(value == sort(value, decreasing = FALSE)[2], "*", "")
          )
      } else {
        best_methods <- .x |>
          dplyr::group_by(.dataset) |>
          dplyr::mutate(
            is_best = ifelse(value == max(value), "**", ""),
            is_2nd_best = ifelse(value == sort(value, decreasing = TRUE)[2], "*", "")
          )
      }
      plt <- ggplot2::ggplot(.x) +
        ggplot2::aes(
          x = .dataset, y = value, fill = .method_name
        ) +
        ggplot2::geom_bar(
          linewidth = 0.1,
          color = "grey98",
          stat = "identity", 
          position = ggplot2::position_dodge(width = 0.8),
          width = 0.8
        ) +
        ggplot2::geom_text(
          data = best_methods,
          ggplot2::aes(label = is_best, color = .method_name, group = .method_name),
          position = ggplot2::position_dodge(width = 0.82),
          hjust = -0.2,
          vjust = 0.8,
          size = 5
        ) +
        ggplot2::geom_text(
          data = best_methods,
          ggplot2::aes(label = is_2nd_best, color = .method_name, group = .method_name),
          position = ggplot2::position_dodge(width = 0.82),
          hjust = -0.2,
          vjust = 0.8,
          size = 5
        ) +
        ggplot2::scale_fill_manual(
          values = method_colors
        ) +
        ggplot2::scale_color_manual(
          values = method_colors
        ) +
        ggplot2::labs(
          x = "Dataset",
          y = metric_names[.y$metric[[1]]],
          fill = "Method"
        ) +
        ggplot2::scale_y_continuous(
          expand = ggplot2::expansion(mult = c(0, 0.1))
        ) +
        ggplot2::coord_flip() +
        vthemes::theme_vmodern(size_preset = "large") +
        ggplot2::theme(
          legend.position = "bottom",
          legend.text = ggplot2::element_text(margin = ggplot2::margin(l = 4, r = 10))
        ) +
        ggplot2::guides(
          fill = ggplot2::guide_legend(reverse = TRUE),
          color = "none"
        )
      if (.y$metric[[1]] != "multi_roc") {
        plt <- plt +
          ggplot2::theme(
            axis.title.y = ggplot2::element_blank()
          )
      }
      return(plt)
    }
  )
  
plt <- patchwork::wrap_plots(
  plt_ls,
  widths = c(1, 1),
  nrow = 1, guides = "collect"
) &
  ggplot2::theme(
    legend.position = "bottom"
  ) +
  patchwork::plot_layout(
    axis_titles = "collect"
  )

# Save
ggplot2::ggsave(
  plt,
  filename = file.path(figures_dir, "evaluation_rf_flipped.pdf"),
  width = 14, height = 6
)

# Calculate and save ranks
rank_df <- plt_df_wide |>
  dplyr::group_by(.dataset) |>
  dplyr::mutate(
    dplyr::across(
      tidyselect::any_of(setdiff(names(metric_names), c("diff_prob"))),
      ~ rank(-.x)
    ),
    diff_prob = rank(diff_prob)
  )
for (metric_name in c("multi_roc", "diff_prob")) {
  rank_summary <- rank_df |>
    dplyr::rename(
      .metric := !!rlang::sym(metric_name)
    ) |>
    dplyr::group_by(.method_name) |>
    dplyr::summarise(
      `Rank 1` = sum(.metric == 1),
      `Rank 2` = sum(.metric == 2),
      `Rank 3` = sum(.metric == 3),
      `Rank 4` = sum(.metric == 4),
      `Rank 5` = sum(.metric == 5),
      `Average Rank` = mean(.metric),
      .groups = "drop"
    ) |>
    dplyr::rename(
      Method = .method_name
    )
  write.csv(
    rank_summary,
    file = file.path(figures_dir, sprintf("evaluation_rf_%s_ranks.csv", metric_name)),
    row.names = FALSE
  )
}


#### Plot local LCMC curves ####
plt_ls <- eval_tib |>
  dplyr::select(.dataset, .method_name, lcmc_local = lcmc_local_raw) |>
  tidyr::unnest(lcmc_local) |>
  dplyr::mutate(
    .method_name = factor(.method_name, levels = names(method_colors)),
    .dataset = factor(.dataset, levels = names(datasets))
  ) |>
  dplyr::group_by(.dataset) |>
  dplyr::group_map(
    function(.x, .y) {
      ggplot2::ggplot(.x) +
        ggplot2::aes(
          x = k, y = lcmc, color = .method_name, group = .method_name
        ) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::geom_point(size = 1.5) +
        ggplot2::scale_color_manual(
          values = method_colors
        ) +
        ggplot2::labs(
          x = "Neighborhood Size (k)",
          y = "LCMC",
          color = "Method",
          title = .y$.dataset[[1]]
        ) +
        vthemes::theme_vmodern(size_preset = "medium") +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, face = "italic")
        )
    }
  )
plt <- patchwork::wrap_plots(
  plt_ls, ncol = 3, byrow = TRUE, guides = "collect"
) +
  patchwork::plot_layout(axis_titles = "collect")
ggplot2::ggsave(
  plt,
  filename = file.path(figures_dir, "evaluation_lcmc_curves.pdf"),
  width = 14, height = 8
)


##### plot method stability #####

######### Figure 7, method stability plots with 2 base methods, 8EQ data #########
#### Get embedding results ####
stability_dir = file.path(results_dir, "stability", "method") 
stability_method_names <- c(
  "PHATE (# neighbors = 50)", "PHATE (# neighbors = 30)", "HLLE",
  "CoMDS_base", "CoMDS_extra",
  "LoCoMDS_base", "LoCoMDS_extra",
  "Multi-SNE (perplexity = 30)_base", "Multi-SNE (perplexity = 30)_extra"
)

comds_out <- readRDS(
  file.path(stability_dir, datasets[["8EQ"]], "CoMDS_stability.rds")
)
locomds_out <- readRDS(
  file.path(stability_dir, datasets[["8EQ"]], "LoCoMDS_stability.rds")
)
multitsne_out <- readRDS(
  file.path(stability_dir, datasets[["8EQ"]], "Multi-SNE(perplexity=30)_stability.rds")
)

dr_list_all_stability <- c(
  list(
    "PHATE (# neighbors = 50)" = dr_list_all[["8EQ"]][["PHATE (# neighbors = 50)"]],
    "PHATE (# neighbors = 30)" = dr_list_all[["8EQ"]][["PHATE (# neighbors = 30)"]],
    "HLLE" = dr_list_all[["8EQ"]][["HLLE"]],
    "CoMDS_base" = comds_out$`base_number_2`,
    "CoMDS_extra" = comds_out$`extra_number_2`,
    "LoCoMDS_base" = locomds_out$`base_number_2`,
    "LoCoMDS_extra" = locomds_out$`extra_number_2`,
    "Multi-SNE (perplexity = 30)_base" = multitsne_out$`base_number_2`,
    "Multi-SNE (perplexity = 30)_extra" = multitsne_out$`extra_number_2`
  ) 
)[stability_method_names] |>
  purrr::map(
    function(.x) {
      # normalize to [0, 1] range for plotting
      for (j in 1:ncol(.x)) {
        .x[, j] <- (.x[, j] - min(.x[, j])) / (max(.x[, j]) - min(.x[, j]))
      }
      return(as.data.frame(.x))
    }
  )

point_size <- 0.5
plt_ls <- purrr::imap(
  dr_list_all_stability[stability_method_names],
  function(plt_df, method_name) {
    plt <- plot_dr(
      plt_df,
      color = as.factor(data_ls[["8EQ"]]$info),
      point_size = point_size
    ) |>
      add_dr_plot_theme(
        method_name = NULL,
      ) +
      ggplot2::scale_color_manual(
        values = eq_cluster_colors
      )
    return(plt)
  }
)
x1 <- clean_method_names("PHATE (# neighbors = 50)", multiline = FALSE)
x2 <- clean_method_names("PHATE (# neighbors = 30)", multiline = FALSE)
x3 <- clean_method_names("HLLE", multiline = FALSE)
input_label <- bquote(bold("Input Methods"))
row1_label <- bquote(atop(bold(.(x1) * ","), bold(.(x2))))
row2_label <- bquote(atop(bold(.(x3) * ", " * .(x1) * ","), bold(.(x2))))

all_row_grobs <- list(
  patchwork::wrap_elements(grid::textGrob(
    input_label, gp = grid::gpar(fontsize = 12, fontface = "bold"), rot = 90
    )),
  patchwork::wrap_elements(grid::textGrob(
    row1_label, gp = grid::gpar(fontsize = 12, fontface = "bold"), rot = 90
    )),
  patchwork::wrap_elements(grid::textGrob(
    row2_label, gp = grid::gpar(fontsize = 12, fontface = "bold"), rot = 90
    ))
  )
y1 <- clean_method_names("PHATE (# neighbors = 50)", multiline = TRUE)
y2 <- clean_method_names("PHATE (# neighbors = 30)", multiline = TRUE)
y3 <- clean_method_names("HLLE", multiline = TRUE)
y1_label <- bquote(bold(.(y1)))
y2_label <- bquote(bold(.(y2)))
y3_label <- bquote(bold(.(y3)))
input_col_grobs <- list(
  patchwork::wrap_elements(grid::textGrob(
      y1_label, gp = grid::gpar(fontsize = 12, fontface = "bold"), rot = 0
    )),
    patchwork::wrap_elements(grid::textGrob(
      y2_label, gp = grid::gpar(fontsize = 12, fontface = "bold"), rot = 0
    )),
    patchwork::wrap_elements(grid::textGrob(
      y3_label, gp = grid::gpar(fontsize = 12, fontface = "bold"), rot = 0
    ))
  )

consensus_column_titles <- list(
  expression(bold(atop("CoMDS", rho == 0.527))),
  expression(bold(atop("LoCoMDS", rho == 0.975))),
  expression(bold(atop("Multi-SNE", rho == 0.708)))
)

consensus_col_grobs <- lapply(consensus_column_titles, function(expr) {
  patchwork::wrap_elements(
    grid::textGrob(label = expr, gp = grid::gpar(fontsize = 12), rot = 0)
    )
  })

plt_ls_with_titles <- c(
  patchwork::wrap_elements(grid::nullGrob()),
  input_col_grobs[[1]], input_col_grobs[[2]], input_col_grobs[[3]],
  all_row_grobs[[1]],  plt_ls[1],  plt_ls[2],  plt_ls[3],
  patchwork::wrap_elements(grid::nullGrob()),
  consensus_col_grobs[[1]], consensus_col_grobs[[2]], consensus_col_grobs[[3]],
  all_row_grobs[[2]], plt_ls[4], plt_ls[6], plt_ls[8],
  all_row_grobs[[3]], plt_ls[5], plt_ls[7], plt_ls[9]
)
design <- "ABCD
           EFGH
           IJKL
           MNOP
           QRST"

plt <- patchwork::wrap_plots(
  plt_ls_with_titles, 
  byrow = TRUE, 
  guides = "collect",
  heights = c(0.3, 1, 0.3, 1, 1),
  widths = c(0.5, 1, 1, 1),
  design = design
  ) &
  ggplot2::theme(
    plot.margin = ggplot2::margin(l = 10, r = 0, t = 5, b = 5)
  )

ggplot2::ggsave(
  plt,
  filename = file.path(
    figures_dir,
    sprintf("method_stability_%s.pdf", datasets[["8EQ"]])
  ),
  width = 0.5 + 3 * 2.0 + 3,
  height = 0.5 + 3 * 2.0 
)



##### Appendix D.3, plot method stability for all base numbers #####
# Plot stability for all base numbers.
stability_dir = file.path(results_dir, "stability", "method") 
method_stability_datasets = c("Star", "8EQ")
dr_list_all_stability <- list()
for (data_name in method_stability_datasets) {
  comds_out <- readRDS(
    file.path(stability_dir, datasets[[data_name]], "CoMDS_stability.rds")
  )
  locomds_out <- readRDS(
    file.path(stability_dir, datasets[[data_name]], "LoCoMDS_stability.rds")
  )
  multitsne_out <- readRDS(
    file.path(stability_dir, datasets[[data_name]], "Multi-SNE(perplexity=30)_stability.rds")
  )
  if (data_name == "8EQ"){
    keep_methods <- c("PHATE (# neighbors = 50)", "PHATE (# neighbors = 30)",
    "kPCA (sigma = 0.001)", "HLLE")
  } else if (data_name == "Star"){
    keep_methods <- c("Isomap", "kPCA (sigma = 0.01)", "kPCA (sigma = 0.001)", "HLLE")
  }

  dr_list_all_stability[[data_name]] <- c(
    list(
      "Input Methods" = dr_list_all[[data_name]][keep_methods],
      "CoMDS" = comds_out[1:4],
      "LoCoMDS" = locomds_out[1:4],
      "Multi-SNE (perplexity = 30)" = multitsne_out[1:4]
      
    ) 
  )
  dr_list_all_stability <- purrr::imap(
    dr_list_all_stability,
    function(method_list, data_name) {
      purrr::map(
        method_list,
        function(.x) {
          # normalize to [0, 1] range for plotting
          purrr::map(
            .x,
            function(emb_df) {
              for (j in 1:ncol(emb_df)) {
                emb_df[, j] <- (emb_df[, j] - min(emb_df[, j])) / (max(emb_df[, j]) - min(emb_df[, j]))
              }
              return(as.data.frame(emb_df))
            }
          )
        }
      )
    }
  )
}




######### Method Stability Plots #########

point_size <- 0.5

for (data_name in method_stability_datasets) {
  
  this_cluster_colors <- if (data_name == "Trajectory") {
    ordered_cluster_colors
  } else if (data_name == "8EQ") {
    eq_cluster_colors
  } else {
    cluster_colors
  }
  if (data_name == "8EQ"){
    keep_methods <- c("PHATE (# neighbors = 50)", "PHATE (# neighbors = 30)",
    "kPCA (sigma = 0.001)", "HLLE")
  } else if (data_name == "Star"){
    keep_methods <- c("Isomap", "kPCA (sigma = 0.01)", "kPCA (sigma = 0.001)", "HLLE")
  }
  
  plt_ls <- purrr::imap(
      dr_list_all_stability[[data_name]],
      function(method_list, method_name) {
        purrr::map(
          method_list,
          function(plt_df) {
            plot_dr(
              plt_df,
              color = as.factor(data_ls[[data_name]]$info),
              point_size = point_size
            ) |>
          add_dr_plot_theme(method_name = NULL) +
          ggplot2::scale_color_manual(values = this_cluster_colors)
          }
        )
      }
  )
  n_rows <- length(plt_ls[[1]])
  n_cols <- length(plt_ls)
  plt_ls_colwise <- purrr::map(seq_len(n_rows), function(i) {
      purrr::map(plt_ls, ~ .x[[i]])
  }) |> purrr::flatten()
    
  base_counts <- 2:3
  if (data_name %in% c("Swiss Roll")) {
    row_labels <- as.vector(rbind(
    paste(base_counts, "base"),
    paste(base_counts, "base + PCA")
    ))
  } else {
    row_labels <- as.vector(rbind(
    paste(base_counts, "base"),
    paste(base_counts, "base + HLLE")
    ))
  }
  
  temp = sapply(keep_methods, clean_method_names, multiline = TRUE)
  left_row_grobs <- purrr::map(temp, ~ patchwork::wrap_elements(
    grid::textGrob(.x, gp = grid::gpar(fontsize = 12, fontface = "bold"), rot = 90)
  ))
  temp2 = sapply(row_labels, clean_method_names, multiline = TRUE)
  right_row_grobs <- purrr::map(temp2, ~ patchwork::wrap_elements(
    grid::textGrob(.x, gp = grid::gpar(fontsize = 12, fontface = "bold"), rot = 90)
  ))
  
  column_titles <- names(plt_ls)
  temp3 <- sapply(column_titles, clean_method_names, multiline = TRUE)
  col_grobs <- purrr::map(temp3, ~ patchwork::wrap_elements(
    grid::textGrob(.x, gp = grid::gpar(fontsize = 12, fontface = "bold"), rot = 0)
  ))

  n_rows <- length(plt_ls[[1]])
n_cols <- length(plt_ls)
plt_ls_colwise <- purrr::map(seq_len(n_rows), function(i) {
  purrr::map(plt_ls, ~ .x[[i]])
}) |> purrr::flatten()

# Combine left label + plots + right label
row_grobs <- purrr::map2(
  purrr::map2(left_row_grobs, split(plt_ls_colwise, rep(1:n_rows, each = n_cols)), ~ c(.x, .y)),
  right_row_grobs,
  ~ c(.x, .y)
)

# Move last column (right labels) to column 3
row_grobs_moved <- purrr::map(row_grobs, function(r) {
  if(length(r) >= 3) {
    right_label <- r[[length(r)]]
    r <- append(r[-length(r)], right_label, after = 2)  # insert after column 2
  }
  r
}) |> purrr::flatten()

plt_ls_with_titles <- c(
  patchwork::wrap_elements(grid::nullGrob()),  # top-left empty
  col_grobs[1],
  patchwork::wrap_elements(grid::nullGrob()),  # optional top-right empty
  col_grobs[2:4],
  row_grobs_moved
)

  plt <- patchwork::wrap_plots(
    plt_ls_with_titles,
    ncol = n_cols + 2,
    byrow = TRUE,
    heights = c(0.3, rep(1, n_rows)),
    widths = c(0.3, 1, 0.2, 1, 1, 1)
  ) &
    ggplot2::theme(plot.margin = ggplot2::margin(l = 10, r = 0, t = 5, b = 5))

  plt <- plt + patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(
      legend.position = "right",
      legend.box = "vertical"
    )

  ggplot2::ggsave(
    plt,
    filename = file.path(figures_dir, sprintf("method_stability_%s_all.pdf", datasets[[data_name]])),
    width = 0.3 + n_cols * 2.0 + 3,
    height = 0.3 + n_rows * 2.0
  )
}

######### Parameter Stability Plots #########

#### Get embedding results ####
stability_dir = file.path(results_dir, "stability", "param") 
method_stability_datasets = c("HIV", "8EQ")
stability_method_names <- c(
  "UMAP (# neighbors = 10)",
  "UMAP (# neighbors = 100)",
  "CoMDS",
  "LoCoMDS",
  "Meta-Spec (UMAP)",
  "Multi-SNE (perplexity = 30)"
)
dr_list_all_stability <- list()

for (data_name in method_stability_datasets) {
  all_out <- dr_out <- readRDS(
    file.path(stability_dir, datasets[[data_name]], "umap_plot.rds")
  )
  dr_list_all_stability[[data_name]] <- c(
    list(
      "UMAP (# neighbors = 10)" = all_out[["UMAP (# neighbors = 10)"]],
      "UMAP (# neighbors = 100)" = all_out[["UMAP (# neighbors = 100)"]],
      "CoMDS" = all_out[["CoMDS"]],
      "LoCoMDS" = all_out[["LoCoMDS"]],
      "Meta-Spec (UMAP)" = all_out[["Meta-Spec (UMAP)"]],
      "Multi-SNE (perplexity = 30)" = all_out[["Multi-SNE (perplexity = 30)"]]
    ) 
  )[stability_method_names] |>
    purrr::map(
      function(.x) {
        # normalize to [0, 1] range for plotting
        for (j in 1:(ncol(.x)-1)) {
          .x[, j] <- (.x[, j] - min(.x[, j])) / (max(.x[, j]) - min(.x[, j]))
        }
        return(as.data.frame(.x))
      }
    )
}


#### Plot stability results with respect to UMAP with different parameters ####
# Figure 8, and Figure A.11
point_size <- 0.8

keep_methods <- stability_method_names
for (data_name in method_stability_datasets) {
  if (data_name == "Trajectory") {
    this_cluster_colors <- ordered_cluster_colors
    plt_ls <- purrr::imap(
      dr_list_all_stability[[data_name]][keep_methods],
      function(plt_df, method_name) {
        plt <- plot_outlier(
          plt_df,
          color = as.factor(data_ls[[data_name]]$info),
          point_size = point_size,
          transparency = 0.2,
          label_size = point_size*5
        ) 
        plt$layers <- purrr::keep(
          plt$layers,
          ~ !inherits(.x$geom, "GeomTextRepel")
        )
        plt_df2 <- plt_df
        plt_df2$.color <- as.factor(data_ls[[data_name]]$info)
        plt_df2$.outlier <- as.numeric(plt_df2$outlier)
        plt <- plt +
        ggplot2::geom_point(
          data = dplyr::filter(plt_df2, .outlier == 0),
          ggplot2::aes(
            color = .color
          ),
          size = point_size,
          alpha = 0.4,
          shape = 16
        ) +
        ggplot2::geom_point(
          data = dplyr::filter(plt_df2, .outlier == 1),
          ggplot2::aes(
            color = .color
          ),
          size = point_size * 3,
          alpha = 1,
          shape = 8,
          show.legend = FALSE
        ) +
        ggplot2::geom_point(
          data = dplyr::filter(plt_df2, .outlier == 1),
          ggplot2::aes(color = .color),
          size = point_size * 10,
          shape = 1,
          stroke = 1
        )

        plt <- plt|>
          add_dr_plot_theme(
            method_name = clean_method_names(method_name),
          ) +
          ggplot2::scale_color_manual(
            values = this_cluster_colors
          ) +
          ggplot2::labs(
            color = NULL
          )
        return(plt)
      }
    )
  } else if (data_name == "8EQ") {
    this_cluster_colors <- eq_cluster_colors
    plt_ls <- purrr::imap(
      dr_list_all_stability[[data_name]][keep_methods],
      function(plt_df, method_name) {
        plt <- plot_outlier(
          plt_df,
          color = as.factor(data_ls[[data_name]]$info),
          point_size = point_size,
          transparency = 0.05,
          label_size = point_size*5
        ) 
        plt$layers <- purrr::keep(
          plt$layers,
          ~ !inherits(.x$geom, "GeomTextRepel")
        )
        plt_df2 <- plt_df
        plt_df2$.color <- as.factor(data_ls[[data_name]]$info)
        plt_df2$.outlier <- as.numeric(plt_df2$outlier)
        plt <- plt +
        ggplot2::geom_point(
          data = dplyr::filter(plt_df2, .outlier == 0),
          ggplot2::aes(
            color = .color
          ),
          size = point_size,
          alpha = 0.05,
          shape = 16
        ) +
        ggplot2::geom_point(
          data = dplyr::filter(plt_df2, .outlier == 1),
          ggplot2::aes(
            color = .color
          ),
          size = point_size * 3,
          alpha = 1,
          shape = 8,
          show.legend = FALSE
        ) +
        ggplot2::geom_point(
          data = dplyr::filter(plt_df2, .outlier == 1),
          ggplot2::aes(color = .color),
          size = point_size * 10,
          shape = 1,
          stroke = 1
        )

        plt <- plt|>
          add_dr_plot_theme(
            method_name = clean_method_names(method_name),
          ) +
          ggplot2::scale_color_manual(
            values = this_cluster_colors
          ) + 
          ggplot2::labs(
            color = NULL
          )
        return(plt)
      }
    )
  }  else if (data_name == "HIV") {
    this_cluster_colors <- cluster_colors
    plt_ls <- purrr::imap(
      dr_list_all_stability[[data_name]][keep_methods],
      function(plt_df, method_name) {
        plt <- plot_outlier(
          plt_df,
          color = as.factor(data_ls[[data_name]]$info),
          point_size = point_size,
          transparency = 0.05,
          label_size = point_size*5
        ) 
        plt$layers <- purrr::keep(
          plt$layers,
          ~ !inherits(.x$geom, "GeomTextRepel")
        )
        plt_df2 <- plt_df
        plt_df2$.color <- as.factor(data_ls[[data_name]]$info)
        plt_df2$.outlier <- as.numeric(plt_df2$outlier)
        plt <- plt +
        ggplot2::geom_point(
          data = dplyr::filter(plt_df2, .outlier == 0),
          ggplot2::aes(
            color = .color
          ),
          size = point_size,
          alpha = 0.05,
          shape = 16
        ) +
        ggplot2::geom_point(
          data = dplyr::filter(plt_df2, .outlier == 1),
          ggplot2::aes(
            color = .color
          ),
          size = point_size * 3,
          alpha = 1,
          shape = 8,
          show.legend = FALSE
        ) +
        ggplot2::geom_point(
          data = dplyr::filter(plt_df2, .outlier == 1),
          ggplot2::aes(color = .color),
          size = point_size * 10,
          shape = 1,
          stroke = 1
        )

        plt <- plt|>
          add_dr_plot_theme(
            method_name = clean_method_names(method_name),
          ) +
          ggplot2::scale_color_manual(
            values = this_cluster_colors
          ) + 
          ggplot2::labs(
            color = NULL
          )
        return(plt)
      }
    )
  } else {
    this_cluster_colors <- cluster_colors
  }
  
  plt_ls[[1]] <- plt_ls[[1]] +
    ggplot2::labs(tag = "(A)") +
    tag_plot_theme
  plt_ls[[3]] <- plt_ls[[3]] +
    ggplot2::labs(tag = "(B)") +
    tag_plot_theme
  plt_ls <- c(
    plt_ls[1:2],
    vdivider,vdivider,
    plt_ls[3:6]
  )
  design <- "ACEG
             BDFH"
  plt <- patchwork::wrap_plots(
    plt_ls, byrow = TRUE, guides = "collect",
    heights = c(1, 1),
    widths = c(1, 0.1, 1, 1),
    design = design
  ) &
    ggplot2::theme(
      plot.margin = ggplot2::margin(l = 10, r = 0, t = 5, b = 5)
    )
  ggplot2::ggsave(
    plt,
    filename = file.path(
      figures_dir,
      sprintf("umap_stability_%s_circle.pdf", datasets[[data_name]])
    ),
    width = 10, height = 6
  )
}
