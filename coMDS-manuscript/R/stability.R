check_correlation_Procrustes <- function(emb1, emb2){
  # check correlation between two embeddings using Procrustes analysis
  proc_fit <- vegan::procrustes(emb1, emb2)
  proc_test <- vegan::protest(emb1, emb2)
  return(list(statistic = proc_test$t0, p_value = proc_test$signif))
}

check_correlation_Mantel <- function(emb1, emb2){
  # check correlation between two embeddings using Mantel test
  dist1 <- stats::dist(emb1)
  dist2 <- stats::dist(emb2)
  mantel_test <- vegan::mantel(dist1, dist2, permutations = 999)
  return(list(statistic = mantel_test$statistic, p_value = mantel_test$signif))
}

check_correlation_subspace_corr <- function(emb1, emb2){
  svd_1 <- svd(emb1)
  svd_2 <- svd(emb2)
  rank_1 <- sum(svd_1$d > 1e-8)
  rank_2 <- sum(svd_2$d > 1e-8)
  if (rank_1 != rank_2) {
    stop("Embeddings have different ranks; cannot compute subspace correlation.")
  }
  basis_1 <- svd_1$u[, 1:rank_1, drop = FALSE]
  basis_2 <- svd_2$u[, 1:rank_2, drop = FALSE]
  M <- t(basis_1) %*% basis_2
  svd_M <- svd(M, nu = 0, nv = 0)$d
  subspace_corr <- mean(svd_M^2)
  return(list(statistic = subspace_corr, p_value = NA))
}

# stability assessment function
run_stability <- function(dgp_out, dr_methods_list, method_name, fit_fun, method_ordered, max_base = 3,
                          stability_dir = "./", point_size = 1, 
                          embedding_accessor = function(res) res$embeddings, 
                          preprocess_fun = function(dr_list) dr_list) {
  number_base <- 2
  plt_data_method_ls <- list()
  iter_rows <- list()

  while(number_base <= max_base){
    methods_base  <- method_ordered[1:number_base]
    methods_extra <- c(methods_base, method_ordered[length(method_ordered)])
    dr_list_base  <- fit_dimension_reduction_methods(data = dgp_out$X, methods_list = dr_methods_list[methods_base])
    dr_list_extra <- fit_dimension_reduction_methods(data = dgp_out$X, methods_list = dr_methods_list[methods_extra])
    base_input  <- preprocess_fun(dr_list_base)
    extra_input <- preprocess_fun(dr_list_extra)
    base_result  <- fit_fun(base_input)
    extra_result <- fit_fun(extra_input)
  
    base_embedding  <- embedding_accessor(base_result)
    extra_embedding <- embedding_accessor(extra_result)
    plt_data_method_ls[[paste0("base_number_", number_base)]]  <- base_embedding
    plt_data_method_ls[[paste0("extra_number_", number_base)]] <- extra_embedding
    procrustes_result <- check_correlation_Procrustes(base_embedding, extra_embedding)
    iter_rows[[number_base - 1]] <- data.frame(
      "Base Number" = number_base,
      "Procrustes Correlation" = procrustes_result$statistic,
      "Procrustes p_value"     = procrustes_result$p_value,
      check.names = FALSE
    )
    mantel_result <- check_correlation_Mantel(base_embedding, extra_embedding)
    iter_rows[[number_base - 1]] <- cbind(
      iter_rows[[number_base - 1]],
      data.frame(
        "Mantel Correlation" = mantel_result$statistic,
        "Mantel p_value"     = mantel_result$p_value,
        check.names = FALSE
      )
    )
    subspace_result <- check_correlation_subspace_corr(base_embedding, extra_embedding)
    iter_rows[[number_base - 1]] <- cbind(
      iter_rows[[number_base - 1]],
      data.frame(
        "Subspace Correlation" = subspace_result$statistic,
        check.names = FALSE
      )
    )
    number_base <- number_base + 1
  }
  saveRDS(plt_data_method_ls, file.path(stability_dir, sprintf("%s_stability.rds", (gsub(" ", "", method_name)))))
  corr_df <- do.call(rbind, iter_rows)

  write.csv(
    corr_df,
    file = file.path(stability_dir, paste0(tolower(method_name), sprintf("_%s_stability_correlation.csv", method_name))),
    row.names = FALSE
  )
  this_plt <- purrr::imap(
    plt_data_method_ls,
    ~ plot_dr(.x, color = dgp_out$labels, point_size = point_size, blank_axes = TRUE) +
      ggplot2::labs(title = .y)
  ) |>
    patchwork::wrap_plots(ncol = 2, byrow = TRUE, guides = "collect") +
    patchwork::wrap_table(corr_df[, c("Base Number", "Procrustes Correlation", "Mantel Correlation", "Subspace Correlation")]) +
    patchwork::plot_annotation(title = paste0(method_name, " Stability Assessment"))
  ggplot2::ggsave(
    this_plt,
    filename = file.path(stability_dir, paste0(tolower(method_name), "_plots.pdf")),
    width =  8, height = length(plt_data_method_ls) / 2 * 3 + 2
  )
  return()
}
