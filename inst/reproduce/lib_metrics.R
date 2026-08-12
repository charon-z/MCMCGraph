# lib_metrics.R -----------------------------------------------------------
# Clustering-quality metrics with label alignment (LSAP). Adapted from the
# project benchmark so the reproduce pipeline uses an identical metric
# definition to the paper's existing tables.
#
# Requires: mclust (adjustedRandIndex), clue (solve_LSAP).

suppressPackageStartupMessages({
  library(mclust)
  library(clue)
})

# Remap predicted labels onto truth labels to maximise agreement (LSAP).
mcg_remap_labels <- function(z_true, z_pred) {
  tbl <- table(z_true, z_pred)
  true_labs <- as.integer(rownames(tbl))
  pred_labs <- as.integer(colnames(tbl))
  K_true <- length(true_labs); K_pred <- length(pred_labs)
  map_to <- setNames(pred_labs, pred_labs)
  if (K_true <= K_pred) {
    assign <- clue::solve_LSAP(tbl, maximum = TRUE)
    for (i in seq_along(assign)) map_to[as.character(pred_labs[assign[i]])] <- true_labs[i]
  } else {
    assign <- clue::solve_LSAP(t(tbl), maximum = TRUE)
    for (j in seq_along(assign)) map_to[as.character(pred_labs[j])] <- true_labs[assign[j]]
  }
  as.integer(map_to[as.character(z_pred)])
}

# Full metric vector: ACC, ARI, NMI, Macro-F1, Balanced-ACC, SmallCluster-Recall.
mcg_metrics <- function(z_true, z_pred) {
  stopifnot(length(z_true) == length(z_pred))
  n <- length(z_true)
  z_true <- as.integer(as.factor(z_true))
  z_pred <- as.integer(as.factor(z_pred))
  z_pred_remap <- tryCatch(mcg_remap_labels(z_true, z_pred), error = function(e) z_pred)

  ari_val <- tryCatch(mclust::adjustedRandIndex(z_true, z_pred_remap), error = function(e) NA_real_)

  nmi_val <- tryCatch({
    tbl <- table(z_true, z_pred_remap)
    p_ij <- tbl / n; p_i <- rowSums(p_ij); p_j <- colSums(p_ij)
    H_T <- -sum(p_i * log2(p_i + 1e-16)); H_P <- -sum(p_j * log2(p_j + 1e-16))
    MI  <- sum(p_ij * log2(p_ij / (p_i %*% t(p_j)) + 1e-16), na.rm = TRUE)
    if (abs(H_T + H_P) < 1e-12) 0 else 2 * MI / (H_T + H_P)
  }, error = function(e) NA_real_)

  acc_val <- tryCatch({
    tbl <- table(z_true, z_pred)
    if (nrow(tbl) > ncol(tbl)) tbl <- t(tbl)
    assign <- clue::solve_LSAP(tbl, maximum = TRUE)
    sum(tbl[cbind(seq_along(assign), assign)]) / n
  }, error = function(e) NA_real_)

  macro_f1 <- tryCatch({
    levs <- sort(unique(z_true))
    f1s <- sapply(levs, function(cl) {
      tp <- sum(z_true == cl & z_pred_remap == cl)
      fp <- sum(z_true != cl & z_pred_remap == cl)
      fn <- sum(z_true == cl & z_pred_remap != cl)
      if (tp + fp + fn == 0) return(NA_real_)
      2 * tp / (2 * tp + fp + fn)
    })
    mean(f1s, na.rm = TRUE)
  }, error = function(e) NA_real_)

  bal_acc <- tryCatch({
    levs <- sort(unique(z_true))
    recalls <- sapply(levs, function(cl) sum(z_true == cl & z_pred_remap == cl) / max(1, sum(z_true == cl)))
    mean(recalls, na.rm = TRUE)
  }, error = function(e) NA_real_)

  small_recall <- tryCatch({
    props <- table(z_true) / n
    small <- as.integer(names(props)[props < 0.10])
    if (length(small) == 0) {
      sorted <- names(sort(props)); small <- as.integer(sorted[1:max(1, floor(length(props) / 4))])
    }
    recalls <- sapply(small, function(cl) sum(z_true == cl & z_pred_remap == cl) / max(1, sum(z_true == cl)))
    mean(recalls, na.rm = TRUE)
  }, error = function(e) NA_real_)

  c(ACC = unname(acc_val), ARI = unname(ari_val), NMI = unname(nmi_val),
    Macro_F1 = unname(macro_f1), Balanced_ACC = unname(bal_acc),
    SmallCluster_Recall = unname(small_recall))
}
