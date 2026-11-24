# ============================================================
# 対応なし・1要因・2群のBrunner–Munzel検定（lawstatのt近似）。
# 表示・下位検定の選択・補正管理はtest_helper.Rが担当。
# このファイルの.t_relative_*はt_ats.Rでも使用する内部関数。
# ============================================================

.t_relative_design <- function(experiment, method) {
  valid <- function(x) is.character(x) && length(x) > 0L && !anyNA(x) && all(nzchar(x))
  label <- if (method == "brunner_munzel") "Brunner–Munzel" else "ATS＋MCTP"
  if (!is.list(experiment) || !is.data.frame(experiment$data) ||
      !valid(experiment$id) || length(experiment$id) != 1L ||
      !valid(experiment$iv) || length(experiment$iv) != 1L ||
      !valid(experiment$dv) || anyDuplicated(c(experiment$id, experiment$iv, experiment$dv))) {
    stop(label, "には前処理済みのexperimentと、重複しないID・1要因・従属変数を指定してください。", call. = FALSE)
  }
  if (length(experiment$within) || !identical(as.character(experiment$between), experiment$iv)) {
    stop(label, "は対応なし・1要因の比較に対応しています。", call. = FALSE)
  }
  data <- experiment$data
  columns <- c(experiment$id, experiment$iv, experiment$dv)
  if (anyDuplicated(names(data)) || !all(columns %in% names(data)) || !nrow(data)) {
    stop("必要な列名とデータの行数を確認してください。", call. = FALSE)
  }
  id <- data[[experiment$id]]
  group <- data[[experiment$iv]]
  if (anyNA(id) || anyNA(group)) stop("ID・要因列の欠損を確認してください。", call. = FALSE)
  group <- if (is.factor(group)) droplevels(group) else factor(group)
  k <- nlevels(group)
  if ((method == "brunner_munzel" && k != 2L) || (method == "ats" && k < 3L)) {
    stop(label, if (method == "brunner_munzel") "は2群を指定してください。" else "は3群以上を指定してください。", call. = FALSE)
  }
  if (any(vapply(split(as.character(group), id, drop = TRUE),
      function(x) length(unique(x)) > 1L, logical(1)))) {
    stop("対応なしの計画ですが、同じIDが複数群に含まれています。", call. = FALSE)
  }
  invisible(levels(group))
}

.t_relative_capture <- function(expression) {
  warnings <- messages <- character()
  error <- NULL
  value <- tryCatch(withCallingHandlers(force(expression), warning = function(w) {
    warnings <<- unique(c(warnings, conditionMessage(w)))
    invokeRestart("muffleWarning")
  }, message = function(m) {
    messages <<- unique(c(messages, conditionMessage(m)))
    invokeRestart("muffleMessage")
  }), error = function(e) { error <<- conditionMessage(e); NULL })
  list(value = value, warnings = warnings, messages = messages, error = error)
}

.t_relative_prepare <- function(experiment, alpha, followup_alpha, method) {
  level_names <- .t_relative_design(experiment, method)
  scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!scalar(alpha) || !scalar(followup_alpha) || alpha <= 0 ||
      alpha > followup_alpha || followup_alpha >= 1) {
    stop("0 < alpha ≤ followup_alpha < 1 にしてください。", call. = FALSE)
  }
  id <- experiment$id; iv <- experiment$iv; data <- experiment$data
  setNames(lapply(experiment$dv, function(target) {
    if (!is.numeric(data[[target]]) || any(is.infinite(data[[target]]))) {
      stop("従属変数は有限の数値またはNAにしてください: ", target, call. = FALSE)
    }
    used <- data[!is.na(data[[target]]), c(id, iv, target), drop = FALSE]
    if (!nrow(used)) stop("検定可能な観測がありません: ", target, call. = FALSE)
    aggregated <- stats::aggregate(used[[target]], by = used[c(id, iv)], FUN = mean)
    names(aggregated)[3L] <- target
    aggregated[[iv]] <- factor(aggregated[[iv]], levels = level_names)
    if (any(table(aggregated[[iv]]) < 2L)) {
      stop("すべての群に有効な参加者が2名以上必要です: ", target, call. = FALSE)
    }
    aggregated <- aggregated[order(aggregated[[iv]], as.character(aggregated[[id]])), , drop = FALSE]
    rownames(aggregated) <- NULL
    samples <- split(aggregated[[target]], aggregated[[iv]], drop = FALSE)
    sample_ids <- unique(as.character(aggregated[[id]]))
    excluded <- setdiff(unique(as.character(data[[id]])), sample_ids)
    repeats <- nrow(used) - nrow(unique(used[c(id, iv)]))
    notices <- c(if (length(experiment$covariate)) "この分析では共変量を使用していません。",
      if (repeats) "複数観測を参加者×水準の平均に集約しています。",
      if (anyNA(data[[target]])) paste0("欠損観測を除外: ", sum(is.na(data[[target]])), "行。"),
      if (length(excluded)) paste0("分析対象外の参加者: ", length(excluded), "名。"))
    list(analysis = method, dependent_variable = target, within = character(), between = iv,
      covariate = experiment$covariate, covariate_adjusted = FALSE, paired = FALSE, design = "between",
      levels = level_names, n_levels = length(level_names), n_participants = length(sample_ids),
      group_sizes = table(aggregated[[iv]]), excluded_ids = excluded,
      missing_rows = sum(is.na(data[[target]])), aggregated_rows = repeats,
      analysis_data = aggregated, samples = samples, contrasts = list(),
      alpha = alpha, followup_alpha = followup_alpha, confidence_level = 1 - alpha,
      messages = notices, warnings = character())
  }), experiment$dv)
}

.t_relative_collection <- function(results, experiment, method) {
  for (i in seq_along(results)) {
    results[[i]]$samples <- NULL
    class(results[[i]]) <- c(paste0("analysis_", method, "_result"), "analysis_test_result", "list")
  }
  attr(results, "stage") <- "test"
  attr(results, "type") <- method
  attr(results, "experiment") <- experiment
  class(results) <- c(paste0("analysis_", method), "analysis_test", "list")
  invisible(results)
}

t_brunner_munzel <- function(experiment, alpha = 0.05, followup_alpha = 0.10) {
  prepared <- .t_relative_prepare(experiment, alpha, followup_alpha, "brunner_munzel")
  if (!requireNamespace("lawstat", quietly = TRUE)) {
    stop('Brunner–Munzelにはlawstatが必要です。install.packages("lawstat")を実行してください。', call. = FALSE)
  }
  results <- lapply(prepared, function(result) {
    x <- result$samples[[1L]]; y <- result$samples[[2L]]
    computed <- .t_relative_capture(lawstat::brunner.munzel.test(x, y,
      alternative = "two.sided", alpha = alpha))
    test <- computed$value
    # 推定の向きは常に第1水準X、第2水準Y。退化例でも経験的推定値は保存する。
    pooled <- rank(c(x, y))
    probability <- (mean(pooled[length(x) + seq_along(y)]) - (length(y) + 1) / 2) / length(x)
    result$test <- test
    result$test_id <- "brunner_munzel"
    result$inference <- "t approximation (two-sided)"
    result$package_version <- as.character(utils::packageVersion("lawstat"))
    result$estimate_definition <- "P(X < Y) + 0.5 P(X = Y); X=level1, Y=level2"
    result$estimate_levels <- setNames(result$levels, c("X", "Y"))
    result$relative_probability <- probability
    result$test_table <- data.frame(Effect = experiment$iv, statistic_name = "t",
      statistic = if (is.null(test)) NA_real_ else unname(test$statistic),
      df = if (is.null(test)) NA_real_ else unname(test$parameter),
      relative_effect = probability,
      lower.CL = if (is.null(test)) NA_real_ else unname(test$conf.int[1L]),
      upper.CL = if (is.null(test)) NA_real_ else unname(test$conf.int[2L]),
      p.value = if (is.null(test)) NA_real_ else unname(test$p.value), stringsAsFactors = FALSE)
    result$followup_factors <- character()
    result$compute_followups <- NULL
    result$adjust <- "none"
    result$warnings <- computed$warnings
    result$messages <- c(result$messages, computed$messages)
    result$error_details <- if (!is.null(computed$error)) list(stage = "Brunner–Munzel", message = computed$error) else NULL
    if (!is.finite(result$test_table$p.value)) result$messages <- c(result$messages,
      "分散を推定できない等の理由で検定不能です。別の検定への自動切替は行っていません。",
      computed$error)
    result
  })
  .t_relative_collection(results, experiment, "brunner_munzel")
}
