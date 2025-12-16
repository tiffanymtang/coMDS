rm(list = ls())
set.seed(331)
library(future)
library(optparse)

for (fname in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) {
  source(fname)
}

# command line arguments
option_list <- list(
  make_option(
    "--data_fname", type = "character", default = "traj.RData",
    help = "Name of data file"
  ),
  make_option(
    "--data_dir", type = "character", default = "data",
    help = "Name of data directory"
  ),
  make_option(
    "--data_name", type = "character", default = NULL,
    help = "Name of data set (used for naming the results directory)"
  ),
  make_option(
    "--eps", type = "double", default = 1e-6,
    help = "Tolerance for CoMDS/LoCoMDS optimization convergence"
  ),
  make_option(
    "--itmax", type = "integer", default = 300,
    help = "Maximum number of iterations for CoMDS/LoCoMDS optimization"
  ),
  make_option(
    "--results_dir", type = "character", default = NULL,
    help = "Directory to save results"
  ),
  make_option(
    "--verbose", action = "store_true", default = FALSE,
    help = "Whether to print verbose output during CoMDS/LoCoMDS optimization"
  )
)
# parse the command line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# specify default data name
if (is.null(opt$data_name)) {
  opt$data_name <- tools::file_path_sans_ext(basename(opt$data_fname))
}

# number of cores
n_cores <- Sys.getenv("NSLOTS")
if (n_cores != "") {
  n_cores <- as.integer(n_cores)
  if (n_cores > 1) {
    plan(multicore, workers = n_cores)
  }
} else {
  n_cores <- 5
  plan(multisession, workers = n_cores)
}
cat(sprintf("Using %s cores\n", n_cores))
str(opt)

# results directory
if (is.null(opt$results_dir)) {
  RESULTS_DIR <- here::here("results", opt$data_name)
} else {
  RESULTS_DIR <- file.path(opt$results_dir, opt$data_name)
}
if (!dir.exists(RESULTS_DIR)) {
  dir.create(RESULTS_DIR, recursive = TRUE)
}

# dimension reduction methods list
dr_methods_list <- list(
  "PCA" = fit_pca,
  "MDS" = fit_mds,
  "Non-metric MDS" = fit_imds,
  "Sammon" = fit_sammon,
  "LLE" = fit_lle,
  "HLLE" = fit_hlle,
  "Isomap" = fit_isomap,
  "kPCA (sigma = 0.01)" = purrr::partial(fit_kpca, sigma = 0.01),
  "kPCA (sigma = 0.001)" = purrr::partial(fit_kpca, sigma = 0.001),
  "LEIM" = fit_leim,
  "UMAP (# neighbors = 30)" = purrr::partial(fit_umap, n_neighbors = 30),
  "UMAP (# neighbors = 50)" = purrr::partial(fit_umap, n_neighbors = 50),
  "tSNE (perplexity = 30)" = purrr::partial(fit_tsne, perplexity = 30),
  "tSNE (perplexity = 50)" = purrr::partial(fit_tsne, perplexity = 50),
  "PHATE (# neighbors = 30)" = purrr::partial(fit_phate, n_neighbors = 30),
  "PHATE (# neighbors = 50)" = purrr::partial(fit_phate, n_neighbors = 50)
)

# LoCoMDS hyperparameter grid
taus <- c(10, 5, 1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001)
percentiles <- seq(0.1, 0.9, 0.1)

# plotting options
point_size <- 1

################################ Run Application ###############################
# load in real data
dgp_out <- real_data_dgp(path = file.path(opt$data_dir, opt$data_fname))

# fit dimension reduction methods
dr_list <- fit_dimension_reduction_methods(
  data = dgp_out$X,
  methods_list = dr_methods_list
)
saveRDS(
  dr_list,
  file = file.path(RESULTS_DIR, "dimension_reduction_methods.rds")
)

# fit meta-spec
metaspec_out <- fit_metaspec(data_list = dr_list)
saveRDS(
  metaspec_out,
  file = file.path(RESULTS_DIR, "metaspec.rds")
)

# fit multi-SNE
multisne_out <- fit_multisne(
  data_list = dr_list,
  perplexity_vec = 30,
  max_iter = opt$itmax
)
saveRDS(
  multisne_out,
  file = file.path(RESULTS_DIR, "multisne.rds")
)

# fit CoMDS
dist_list <- purrr::map(dr_list, ~ dist(.x))
comds_out <- fit_comds(
  dist_list,
  eps = opt$eps, itmax = opt$itmax, verbose = opt$verbose
)
saveRDS(
  comds_out,
  file = file.path(RESULTS_DIR, "comds.rds")
)

# fit LoCoMDS
locomds_out <- fit_locomds(
  dist_list, data = dgp_out$X, taus = taus, percentiles = percentiles,
  eps = opt$eps, itmax = opt$itmax, verbose = opt$verbose
)
saveRDS(
  locomds_out,
  file = file.path(RESULTS_DIR, "locomds.rds")
)

############################### Evaluate Methods ###############################
eval_out <- evaluate_dr_methods(
  embed_list = c(
    dr_list,
    list(
      "CoMDS" = comds_out$embeddings,
      "LoCoMDS" = locomds_out$embeddings,
      "Meta-Spec (kPCA)" = metaspec_out$embeddings$kPCA,
      "Meta-Spec (UMAP)" = metaspec_out$embeddings$UMAP,
      "Multi-SNE (perplexity = 30)" = multisne_out$embeddings$perplexity_30
    )
  ),
  data = dgp_out$X,
  labels = dgp_out$labels
)
saveRDS(
  eval_out,
  file = file.path(RESULTS_DIR, "evaluation_results.rds")
)

################################## Plots #######################################
consensus_names <- c(
  "CoMDS", "LoCoMDS", 
  "Meta-Spec (kPCA)", "Meta-Spec (UMAP)", "Multi-SNE (perplexity = 30)"
)

# plot LoCoMDS hyperparameters trace
tuning_plt <- locomds_out$tuning$plot
ggplot2::ggsave(
  tuning_plt,
  filename = file.path(RESULTS_DIR, "locomds_hyperparameter_trace.pdf"),
  width = 12, height = 6
)

# plot LoCoMDS hyperparameters
locomds_plt <- purrr::map(
  locomds_out$fit,
  ~ as.data.frame(.x$gspace) |>
    setNames(paste("Component", 1:ncol(.x$gspace))) |>
    plot_dr(
      color = dgp_out$labels, point_size = point_size, blank_axes = TRUE
    ) +
    ggplot2::labs(
      title = sprintf("tau = %s, percentile = %s", .x$tau, .x$percentile)
    )
) |>
  patchwork::wrap_plots(
    ncol = length(percentiles), byrow = FALSE, guides = "collect"
  )
ggplot2::ggsave(
  locomds_plt,
  filename = file.path(RESULTS_DIR, "locomds_hyperparameter_plots.pdf"),
  width = 2.6 * length(percentiles) + 1, height = 2.4 * length(taus)
)

# plot all dimension reduction results
plt_data_ls <- c(
  dr_list,
  list(
    "CoMDS" = comds_out$embeddings,
    "LoCoMDS" = locomds_out$embeddings,
    "Meta-Spec (kPCA)" = metaspec_out$embeddings$kPCA,
    "Meta-Spec (UMAP)" = metaspec_out$embeddings$UMAP,
    "Multi-SNE (perplexity = 30)" = multisne_out$embeddings$perplexity_30
  )
)
plt <- purrr::imap(
  plt_data_ls,
  ~ plot_dr(
    .x, color = dgp_out$labels, point_size = point_size, blank_axes = TRUE
  ) +
    ggplot2::labs(title = .y)
) |>
  patchwork::wrap_plots(ncol = 4, byrow = TRUE, guides = "collect")
ggplot2::ggsave(
  plt,
  filename = file.path(RESULTS_DIR, "all_plots.pdf"),
  width = 12, height = 16
)

# plot full distribution of evaluation metrics
plt <- eval_out |>
  dplyr::filter(.method_name %in% !!consensus_names) |>
  plot_evaluation_distribution()
ggplot2::ggsave(
  plt,
  filename = file.path(RESULTS_DIR, "evaluation_distribution.pdf"),
  width = 14, height = 16
)

# plot mean of evaluation metrics
plt <- eval_out |>
  dplyr::filter(.method_name %in% !!consensus_names) |>
  plot_evaluation_means()
ggplot2::ggsave(
  plt,
  filename = file.path(RESULTS_DIR, "evaluation_means.pdf"),
  width = 14, height = 18
)

# plot lcmc curves
plt <- eval_out |>
  dplyr::filter(.method_name %in% !!consensus_names) |>
  plot_lcmc_curves()
ggplot2::ggsave(
  plt,
  filename = file.path(RESULTS_DIR, "evaluation_lcmc_curves.pdf"),
  width = 13, height = 6
)
