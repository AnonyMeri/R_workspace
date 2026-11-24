# ============================================================
# 相関分析の計算処理。選択・GUI・表示・p値の補正はhelperが担当。
# 独立した観測を行とする通常の二変量相関。全ペアをstats::cor.test()で検定する。
# Pearson / Spearman / Kendall。CIはcor.test()が返す場合だけ保持する。
# ============================================================

r_correlation <- function(data, variables, method = "pearson",
                          use = "pairwise.complete.obs", alternative = "two.sided",
                          ci = .95, exact = NULL) {
  if (!is.data.frame(data) || anyDuplicated(names(data)) ||
      !is.character(variables) || length(variables) < 2L || anyNA(variables) ||
      anyDuplicated(variables) || any(!variables %in% names(data))) {
    stop("データフレームと、異なる数値変数を2つ以上指定してください。", call. = FALSE)
  }
  if (any(!vapply(data[variables], function(x) is.numeric(x) && is.null(dim(x)), logical(1)))) stop("相関分析には数値変数を指定してください。", call. = FALSE)
  if (any(vapply(data[variables], function(x) any(is.infinite(x)), logical(1)))) {
    stop("分析対象の変数に無限大が含まれています。データを確認してください。", call. = FALSE)
  }
  method <- match.arg(method, c("pearson", "spearman", "kendall"))
  use <- match.arg(use, c("pairwise.complete.obs", "complete.obs"))
  alternative <- match.arg(alternative, c("two.sided", "greater", "less"))
  if (!is.numeric(ci) || length(ci) != 1L || !is.finite(ci) || ci <= 0 || ci >= 1) stop("ciは0より大きく1より小さい数値にしてください。", call. = FALSE)
  if (!is.null(exact) && (!is.logical(exact) || length(exact) != 1L || is.na(exact))) stop("exactにはNULL、TRUE、FALSEのいずれかを指定してください。", call. = FALSE)
  frame <- data[variables]
  common <- if (use == "complete.obs") stats::complete.cases(frame) else rep(TRUE, nrow(frame))
  pairs <- utils::combn(variables, 2L)
  tables <- tests <- diagnostics <- vector("list", ncol(pairs))
  correlation <- p_values <- counts <- matrix(NA_real_, length(variables), length(variables), dimnames = list(variables, variables))
  for (name in variables) {
    x <- frame[[name]][common & !is.na(frame[[name]])]
    counts[name, name] <- length(x)
    if (length(unique(x)) >= 2L) correlation[name, name] <- 1
  }
  for (i in seq_len(ncol(pairs))) {
    left <- pairs[1L, i]; right <- pairs[2L, i]
    rows <- which(common & stats::complete.cases(frame[c(left, right)]))
    x <- frame[[left]][rows]; y <- frame[[right]][rows]
    n <- length(rows)
    warnings <- character(); error <- NULL; test <- NULL
    if (n < if (method == "pearson") 3L else 2L) error <- "検定に必要な有効データ数が不足しています。"
    else if (length(unique(x)) < 2L || length(unique(y)) < 2L) error <- "一方または両方の変数が一定値のため、相関を計算できません。"
    else test <- tryCatch(withCallingHandlers(
      stats::cor.test(x, y, method = method, alternative = alternative, conf.level = ci, exact = exact),
      warning = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")}),
      error = function(e) {error <<- conditionMessage(e); NULL})
    row <- data.frame(variable1 = left, variable2 = right, N = n, estimate = NA_real_,
      statistic_name = NA_character_, statistic = NA_real_, df = NA_real_,
      lower.CL = NA_real_, upper.CL = NA_real_, p.value = NA_real_, stringsAsFactors = FALSE)
    if (!is.null(test)) {
      row$estimate <- unname(test$estimate[1L]); row$p.value <- test$p.value
      if (length(test$statistic)) {row$statistic <- unname(test$statistic[1L]); row$statistic_name <- names(test$statistic)[1L]}
      if (length(test$parameter)) row$df <- unname(test$parameter[1L])
      if (length(test$conf.int) == 2L) {row$lower.CL <- test$conf.int[1L]; row$upper.CL <- test$conf.int[2L]}
    }
    tables[[i]] <- row
    tests[i] <- list(test)
    diagnostics[[i]] <- list(variable1 = left, variable2 = right, rows_used = rows,
      error = error, warnings = unique(warnings))
    correlation[left, right] <- correlation[right, left] <- row$estimate
    p_values[left, right] <- p_values[right, left] <- row$p.value
    counts[left, right] <- counts[right, left] <- n
  }
  result <- list(analysis = "correlation", variables = variables, method = method,
    use = use, alternative = alternative, ci = ci, exact = exact, analysis_data = frame,
    correlation_table = do.call(rbind, tables), correlation_matrix = correlation,
    p.value_matrix = p_values, N_matrix = counts, tests = tests, diagnostics = diagnostics)
  class(result) <- c("analysis_correlation_result", "list")
  result
}
