# ============================================================
# 通常のt検定：対応なしStudent / 対応ありpaired（両側、2水準）。
# 計算だけを担当し、表示・方法選択はtest_helper.Rへ委譲する。
# .t_mean_*はt_welch.Rでも使用する内部関数。
# ============================================================

.t_mean_design <- function(experiment, method) {
  label <- if (method == "t") "t検定" else "Welch"
  valid <- function(x) is.character(x) && length(x) > 0L && !anyNA(x) && all(nzchar(x))
  if (!is.list(experiment) || !is.data.frame(experiment$data) ||
      !valid(experiment$id) || length(experiment$id) != 1L ||
      !valid(experiment$iv) || length(experiment$iv) != 1L ||
      !valid(experiment$dv) || anyDuplicated(c(experiment$id, experiment$iv, experiment$dv))) {
    stop(label, "には前処理済みのexperimentと、重複しないID・1要因・従属変数を指定してください。", call. = FALSE)
  }
  iv <- experiment$iv
  paired <- identical(as.character(experiment$within), iv) && !length(experiment$between)
  independent <- identical(as.character(experiment$between), iv) && !length(experiment$within)
  if ((!paired && !independent) || (method == "welch" && !independent)) {
    stop(if (method == "t") "t検定は1要因の対応あり・対応なしに対応しています。" else
      "Welchは対応なし・1要因の比較に対応しています。", call. = FALSE)
  }
  data <- experiment$data
  required <- c(experiment$id, iv, experiment$dv)
  if (anyDuplicated(names(data)) || !all(required %in% names(data)) || !nrow(data)) {
    stop("必要な列名とデータの行数を確認してください。", call. = FALSE)
  }
  if (anyNA(data[[experiment$id]]) || anyNA(data[[iv]])) stop("ID・要因列の欠損を確認してください。", call. = FALSE)
  group <- if (is.factor(data[[iv]])) droplevels(data[[iv]]) else factor(data[[iv]])
  if ((method == "t" && nlevels(group) != 2L) || (method == "welch" && nlevels(group) < 2L)) {
    stop(if (method == "t") "t検定は2水準のみです。3水準以上はanovaを指定してください。" else
      "Welchには2群以上が必要です。", call. = FALSE)
  }
  if (independent && any(vapply(split(as.character(group), data[[experiment$id]], drop = TRUE),
      function(x) length(unique(x)) > 1L, logical(1)))) {
    stop("対応なしの計画ですが、同じIDが複数群に含まれています。", call. = FALSE)
  }
  list(paired = paired, levels = levels(group))
}

.t_mean_capture <- function(expression) {
  warnings <- messages <- character()
  error <- NULL
  value <- tryCatch(withCallingHandlers(force(expression), warning = function(w) {
    warnings <<- unique(c(warnings, conditionMessage(w))); invokeRestart("muffleWarning")
  }, message = function(m) {
    messages <<- unique(c(messages, conditionMessage(m))); invokeRestart("muffleMessage")
  }), error = function(e) { error <<- conditionMessage(e); NULL })
  list(value = value, warnings = warnings, messages = messages, error = error)
}

.t_mean_prepare <- function(experiment, alpha, followup_alpha, method) {
  design <- .t_mean_design(experiment, method)
  scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!scalar(alpha) || !scalar(followup_alpha) || alpha <= 0 ||
      alpha > followup_alpha || followup_alpha >= 1) {
    stop("0 < alpha ≤ followup_alpha < 1 にしてください。", call. = FALSE)
  }
  data <- experiment$data; id <- experiment$id; iv <- experiment$iv
  setNames(lapply(experiment$dv, function(target) {
    if (!is.numeric(data[[target]]) || any(is.infinite(data[[target]]))) {
      stop("従属変数は有限の数値またはNAにしてください: ", target, call. = FALSE)
    }
    used <- data[!is.na(data[[target]]), c(id, iv, target), drop = FALSE]
    if (!nrow(used)) stop("検定可能な観測がありません: ", target, call. = FALSE)
    aggregated <- stats::aggregate(used[[target]], by = used[c(id, iv)], FUN = mean)
    names(aggregated)[3L] <- target
    aggregated[[iv]] <- factor(aggregated[[iv]], levels = design$levels)
    if (any(table(aggregated[[iv]]) == 0L)) {
      stop("従属変数に有効な観測のない水準があります: ", target, call. = FALSE)
    }
    if (design$paired) {
      counts <- table(as.character(aggregated[[id]]))
      complete_ids <- names(counts)[counts == length(design$levels)]
      aggregated <- aggregated[as.character(aggregated[[id]]) %in% complete_ids, , drop = FALSE]
    }
    if (any(table(aggregated[[iv]]) < 2L)) {
      stop(if (design$paired) "完全な対応を持つ参加者が2名以上必要です: " else
        "すべての群に有効な参加者が2名以上必要です: ", target, call. = FALSE)
    }
    aggregated <- aggregated[order(aggregated[[iv]], as.character(aggregated[[id]])), , drop = FALSE]
    rownames(aggregated) <- NULL
    sample_ids <- unique(as.character(aggregated[[id]]))
    samples <- setNames(lapply(design$levels, function(level) {
      cell <- aggregated[aggregated[[iv]] == level, , drop = FALSE]
      if (design$paired) cell <- cell[match(sample_ids, as.character(cell[[id]])), , drop = FALSE]
      cell[[target]]
    }), design$levels)
    excluded <- setdiff(unique(as.character(data[[id]])), sample_ids)
    repeats <- nrow(used) - nrow(unique(used[c(id, iv)]))
    descriptive <- data.frame(level = design$levels, n = lengths(samples),
      mean = vapply(samples, mean, numeric(1)), SD = vapply(samples, stats::sd, numeric(1)), row.names = NULL)
    notes <- c(if (length(experiment$covariate)) "この分析では共変量を使用していません。",
      if (repeats) "複数観測を参加者×水準の平均に集約しています。",
      if (anyNA(data[[target]])) paste0("欠損観測を除外: ", sum(is.na(data[[target]])), "行。"),
      if (length(excluded)) paste0("分析対象外の参加者: ", length(excluded), "名。"),
      if (design$paired && length(excluded)) "両水準の有効値が揃った参加者のみで分析しています。")
    list(analysis = method, dependent_variable = target, paired = design$paired,
      design = if (design$paired) "within" else "between", within = experiment$within, between = experiment$between,
      levels = design$levels, n_levels = length(design$levels),
      n_participants = length(sample_ids), group_sizes = table(aggregated[[iv]]),
      excluded_ids = excluded, missing_rows = sum(is.na(data[[target]])), aggregated_rows = repeats,
      covariate = experiment$covariate, covariate_adjusted = FALSE,
      analysis_data = aggregated, samples = samples, descriptive = descriptive,
      contrasts = list(), followup_factors = character(), compute_followups = NULL,
      alpha = alpha, followup_alpha = followup_alpha, confidence_level = 1 - alpha,
      alternative = "two.sided", null_difference = 0, adjust = "none",
      package_version = as.character(utils::packageVersion("stats")),
      messages = notes, warnings = character())
  }), experiment$dv)
}

.t_mean_two <- function(result, iv, welch = FALSE) {
  x <- result$samples[[1L]]; y <- result$samples[[2L]]
  computed <- .t_mean_capture(stats::t.test(x, y, paired = result$paired,
    var.equal = !welch, alternative = "two.sided", mu = 0, conf.level = result$confidence_level))
  test <- computed$value
  result$test <- test
  result$test_id <- if (result$paired) "paired_t" else if (welch) "welch_t" else "student_t"
  result$variance_assumption <- if (result$paired) "normality of paired differences" else if (welch)
    "unequal group variances allowed" else "equal group variances"
  result$estimate_definition <- "mean(level1) - mean(level2)"
  result$estimate_levels <- result$levels
  result$mean_difference <- if (result$paired) mean(x - y) else mean(x) - mean(y)
  result$test_table <- data.frame(Effect = iv, statistic_name = "t",
    statistic = if (is.null(test)) NA_real_ else unname(test$statistic),
    df = if (is.null(test)) NA_real_ else unname(test$parameter), estimate = result$mean_difference,
    SE = if (is.null(test)) NA_real_ else unname(test$stderr),
    lower.CL = if (is.null(test)) NA_real_ else unname(test$conf.int[1L]),
    upper.CL = if (is.null(test)) NA_real_ else unname(test$conf.int[2L]),
    p.value = if (is.null(test)) NA_real_ else unname(test$p.value), stringsAsFactors = FALSE)
  result$messages <- c(result$messages, computed$messages)
  result$warnings <- computed$warnings
  result$error_details <- if (!is.null(computed$error)) list(stage = result$test_id, message = computed$error) else NULL
  if (!is.finite(result$test_table$p.value)) result$messages <- c(result$messages,
    "分散を推定できない等の理由で検定不能です。別の検定への自動切替は行っていません。", computed$error)
  result
}

.t_mean_collection <- function(results, experiment, method) {
  for (i in seq_along(results)) {
    results[[i]]$samples <- NULL
    class(results[[i]]) <- c(paste0("analysis_", method, "_result"), "analysis_test_result", "list")
  }
  attr(results, "stage") <- "test"; attr(results, "type") <- method; attr(results, "experiment") <- experiment
  class(results) <- c(paste0("analysis_", method), "analysis_test", "list")
  invisible(results)
}

t_t <- function(experiment, alpha = 0.05, followup_alpha = 0.10) {
  results <- .t_mean_prepare(experiment, alpha, followup_alpha, "t")
  results <- lapply(results, .t_mean_two, iv = experiment$iv, welch = FALSE)
  .t_mean_collection(results, experiment, "t")
}
