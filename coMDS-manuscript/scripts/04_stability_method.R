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

# DR methods results directory
if (is.null(opt$results_dir)) {
  RESULTS_DIR <- here::here("results", opt$data_name)
} else {
  RESULTS_DIR <- file.path(opt$results_dir, opt$data_name)
}

# stability results directory
if (is.null(opt$results_dir)) {
  STABILITY_DIR <- here::here("results", "stability", "method", opt$data_name)
} else {
  STABILITY_DIR <- file.path(opt$results_dir, "stability", opt$data_name)
}
if (!dir.exists(STABILITY_DIR)) {
  dir.create(STABILITY_DIR, recursive = TRUE)
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

# load in real data
if(opt$data_name == "swiss_roll"){
  dgp_out <- readRDS(file.path(RESULTS_DIR, "Swiss_roll.rds"))
} else{
  dgp_out <- real_data_dgp(path = file.path(opt$data_dir, opt$data_fname))
}

# load in dimension reduction results and metaspec results
dr_file <- file.path(RESULTS_DIR, "dimension_reduction_methods.rds")
if (file.exists(dr_file)) {
  dr_list <- readRDS(dr_file)
} else{
  dr_list <- fit_dimension_reduction_methods(
    data = dgp_out$X,
    methods_list = dr_methods_list
  )
}
metaspec_file <- file.path(RESULTS_DIR, "metaspec.rds")
if (file.exists(metaspec_file)) {
  metaspec_out <- readRDS(metaspec_file)
} else{
  metaspec_out <- fit_metaspec(data_list = dr_list)
}
method_names_ordered <- metaspec_out$fit$method_names

# Metaspec kPCA
run_stability(
  dgp_out = dgp_out,
  dr_methods_list = dr_methods_list,
  method_name = "Meta-Spec (kPCA)",
  fit_fun = fit_metaspec,
  method_ordered = method_names_ordered,
  embedding_accessor = function(res) res$embeddings$kPCA,
  stability_dir = STABILITY_DIR
)

# Metaspec UMAP
run_stability(
  dgp_out = dgp_out,
  dr_methods_list = dr_methods_list,
  method_name = "Meta-Spec (UMAP)",
  fit_fun = fit_metaspec,
  method_ordered = method_names_ordered,
  embedding_accessor = function(res) res$embeddings$UMAP,
  stability_dir = STABILITY_DIR
)

# multi-SNE (perplexity = 30)
run_stability(
  dgp_out = dgp_out,
  dr_methods_list = dr_methods_list,
  method_name = "Multi-SNE (perplexity = 30)",
  fit_fun = function(dr_list) fit_multisne(dr_list, perplexity_vec = c(30), max_iter = opt$itmax),
  method_ordered = method_names_ordered,
  embedding_accessor = function(res) res$embeddings$perplexity_30,
  stability_dir = STABILITY_DIR
)

# CoMDS
run_stability(
  dgp_out = dgp_out,
  dr_methods_list = dr_methods_list,
  method_name = "CoMDS",
  fit_fun = function(dist_list) fit_comds(dist_list, eps = opt$eps, itmax = opt$itmax, verbose = opt$verbose),
  method_ordered = method_names_ordered,
  embedding_accessor = function(res) res$embeddings,
  preprocess_fun = function(dr_list) purrr::map(dr_list, ~ dist(.x)),
  stability_dir = STABILITY_DIR
)

# LoCoMDS
run_stability(
  dgp_out = dgp_out,
  dr_methods_list = dr_methods_list,
  method_name = "LoCoMDS",
  fit_fun = function(dist_list) fit_locomds(
    dist_list, data = dgp_out$X, taus = taus, percentiles = percentiles,
    eps = opt$eps, itmax = opt$itmax, verbose = opt$verbose
  ),
  method_ordered = method_names_ordered,
  embedding_accessor = function(res) res$embeddings,
  preprocess_fun = function(dr_list) purrr::map(dr_list, ~ dist(.x)),
  stability_dir = STABILITY_DIR
)

plt_data_all_ls <- c(
  dr_list[c(method_names_ordered[1:5], method_names_ordered[length(method_names_ordered)])]
)

# plot all dimension reduction results
plt <- purrr::imap(
  plt_data_all_ls,
  ~ plot_dr(
    .x, color = dgp_out$labels, point_size = point_size, blank_axes = TRUE
  ) +
    ggplot2::labs(title = .y)
) |>
  patchwork::wrap_plots(ncol = 4, byrow = TRUE, guides = "collect")
ggplot2::ggsave(
  plt,
  filename = file.path(STABILITY_DIR, "stability_inputmethods.pdf"),
  width = 12, height = ceiling(length(plt_data_all_ls) / 4) * 3
)
