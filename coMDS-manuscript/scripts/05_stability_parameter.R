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
  ),
  make_option(
    "--method", type = "character", default = "UMAP",
    help = "Dimensionality reduction method to evaluate"
  ),
  make_option(
    "--stability_params", type = "character", default = "10,100",
    help = "Stability parameter value to evaluate"
  ),
  make_option(
    "--num_outliers", type = "double", default = 5,
    help = "Number of outliers to detect"
  )
)

# parse the command line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
opt$stability_params <- as.integer(strsplit(opt$stability_params, ",")[[1]])

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
  STABILITY_DIR <- here::here("results", "stability", "param", opt$data_name)
} else {
  STABILITY_DIR <- file.path(opt$results_dir, "stability", "param", opt$data_name)
}
if (!dir.exists(STABILITY_DIR)) {
  dir.create(STABILITY_DIR, recursive = TRUE)
}

# dimension reduction methods list
method_map <- list(
  kpca = list(fun = fit_kpca, param_name = "sigma", label = "kPCA (sigma = %s)"),
  umap = list(fun = fit_umap, param_name = "n_neighbors", label = "UMAP (# neighbors = %s)"),
  tsne = list(fun = fit_tsne, param_name = "perplexity", label = "tSNE (perplexity = %s)"),
  phate = list(fun = fit_phate, param_name = "n_neighbors", label = "PHATE (# neighbors = %s)")
)

method_key <- tolower(opt$method)
if (!(method_key %in% names(method_map))) stop("Unknown method: ", opt$method)
method_map <- method_map[[method_key]]

# dimension reduction methods list for input method and stability parameters
dr_methods_list <- purrr::set_names(
  lapply(opt$stability_params,
  function(p) {
    named_arg <- setNames(list(p), method_map$param_name)
    do.call(purrr::partial, c(list(method_map$fun), named_arg))
    }
  ),
  sprintf(method_map$label, opt$stability_params)
)

# LoCoMDS hyperparameter grid
taus <- c(10, 5, 1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001)
percentiles <- seq(0.1, 0.9, 0.1)

# load in original high-dim data
if(opt$data_name == "swiss_roll"){
  dgp_out <- readRDS(file.path(RESULTS_DIR, "Swiss_roll.rds")) 
} else{
  dgp_out <- real_data_dgp(path = file.path(opt$data_dir, opt$data_fname))
}

################################ Run Application ###############################
# fit dimension reduction methods
dr_list <- fit_dimension_reduction_methods(
  data = dgp_out$X,
  methods_list = dr_methods_list
)
saveRDS(
  dr_list,
  file = file.path(STABILITY_DIR, sprintf("dimension_reduction_%s_stability.rds", method_key))
)

# fit meta-spec
metaspec_out <- fit_metaspec(data_list = dr_list)
saveRDS(
  metaspec_out,
  file = file.path(STABILITY_DIR, sprintf("metaspec_%s_stability.rds", method_key))
)

# fit multi-SNE
multisne_out <- fit_multisne(
  data_list = dr_list,
  perplexity_vec = 30,
  max_iter = opt$itmax
)
saveRDS(
  multisne_out,
  file = file.path(STABILITY_DIR, sprintf("multisne_%s_stability.rds", method_key))
)

# fit CoMDS
dist_list <- purrr::map(dr_list, ~ dist(.x))
comds_out <- fit_comds(
  dist_list,
  eps = opt$eps, itmax = opt$itmax, verbose = opt$verbose
)
saveRDS(
  comds_out,
  file = file.path(STABILITY_DIR, sprintf("comds_%s_stability.rds", method_key))
)

# fit LoCoMDS
locomds_out <- fit_locomds(
  dist_list, data = dgp_out$X, taus = taus, percentiles = percentiles,
  eps = opt$eps, itmax = opt$itmax, verbose = opt$verbose
)
saveRDS(
  locomds_out,
  file = file.path(STABILITY_DIR, sprintf("locomds_%s_stability.rds", method_key))
)

dist_list_matrix <- purrr::map(dist_list, ~ as.matrix(.x))

get_outlier_index <- function(dist_list_matrix, num_outliers = opt$num_outliers){
  
  row_spearman <- function(X, Y, method = "spearman"){
    stopifnot(all(dim(X) == dim(Y))) 
    sapply(1:nrow(X), function(i) cor(X[i, ], Y[i, ], method = method))
  }

  pairs <- combn(length(dist_list_matrix), 2)
  corr_matrix <- apply(pairs, 2, function(idx) {
    row_spearman(dist_list_matrix[[idx[1]]], dist_list_matrix[[idx[2]]])
  })

  max_corr <- apply(corr_matrix, 1, max)
  c_val <- sort(max_corr)[num_outliers]
  result <- which(max_corr <= c_val)
  return(result)
}

outlier_index <- get_outlier_index(dist_list_matrix, num_outliers = opt$num_outliers)

plt_data_ls <- c(
  dr_list,
  list(
    "CoMDS" = comds_out$embeddings,
    "LoCoMDS" = locomds_out$embeddings,
    "Meta-Spec (UMAP)" = metaspec_out$embeddings$UMAP,
    "Multi-SNE (perplexity = 30)" = multisne_out$embeddings$perplexity_30
  )
)

# plotting options
point_size <- 2

plt_data_ls <- lapply(
  plt_data_ls,
  function(emb) {
    emb <- as.data.frame(emb)
    emb$outlier <- 0
    emb$outlier[outlier_index] <- 1
    return(emb)
  }
)

saveRDS(
  plt_data_ls,
  file = file.path(STABILITY_DIR, sprintf("%s_plot.rds", method_key))
)

plt <- purrr::imap(
  plt_data_ls,
  ~ plot_outlier(
    .x, point_size = point_size, blank_axes = TRUE, 
    label_size = point_size * 2, max_overlaps = Inf, color = dgp_out$labels
  ) +
    ggplot2::labs(title = .y)
) |>
  patchwork::wrap_plots(ncol = 2, byrow = TRUE, guides = "collect")

ggplot2::ggsave(
  plt,
  filename = file.path(STABILITY_DIR, sprintf("outlier_%s_stability.pdf", method_key)),
  width = 12, height = 16
)