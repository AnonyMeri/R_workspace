# ============================================================
# Welch系：対応なし2群はWelch t、3群以上はWelch一元配置ANOVA。
# helperからt_t.Rに続けて読み込む。共変量・球面性補正は使用しない。
# 下位検定はGames–Howell（全ペア）。p値はPMCMRplus公式実装。
# CIは同じstudentized-range分布の分位点を用いたGames–Howell区間。
# https://search.r-project.org/CRAN/refmans/PMCMRplus/html/gamesHowellTest.html
# 計算時は丸めない。表示・実行基準・補正の管理はhelperが担当。
# ============================================================

.t_welch_followups <- function(samples, levels, iv, confidence_level) {
  if (!requireNamespace("PMCMRplus", quietly = TRUE)) {
    stop('Games–HowellにはPMCMRplusが必要です。install.packages("PMCMRplus")を実行してください。', call. = FALSE)
  }
  # リスト入力でパッケージ内の群名解釈を避け、元の水準順は別途保持。
  computed <- .t_mean_capture(PMCMRplus::gamesHowellTest(unname(samples)))
  if (!is.null(computed$error)) {
    return(setNames(list(list(kind = "games_howell", source_effect = iv, error = computed$error,
      error_details = list(stage = "Games–Howell", message = computed$error),
      warnings = computed$warnings, messages = computed$messages)), iv))
  }
  raw <- computed$value
  pairs <- utils::combn(seq_along(samples), 2L)
  sizes <- lengths(samples); means <- vapply(samples, mean, numeric(1))
  variances <- vapply(samples, stats::var, numeric(1)); k <- length(samples)
  rows <- lapply(seq_len(ncol(pairs)), function(index) {
    left <- pairs[1L, index]; right <- pairs[2L, index]
    components <- variances[c(left, right)] / sizes[c(left, right)]
    SE <- sqrt(sum(components))
    df <- sum(components)^2 / sum(components^2 / (sizes[c(left, right)] - 1))
    difference <- means[left] - means[right]
    # packageの下三角の統計量はright-left。表示・保存はleft-rightで統一。
    q <- -unname(raw$statistic[as.character(right), as.character(left)])
    p <- unname(raw$p.value[as.character(right), as.character(left)])
    half_width <- if (is.finite(df) && df > 0 && is.finite(SE))
      stats::qtukey(confidence_level, nmeans = k, df = df) * SE / sqrt(2) else NA_real_
    data.frame(level1 = levels[left], level2 = levels[right],
      pair_key = paste(left, right, sep = ":"), comparison = paste(levels[left], levels[right], sep = " - "),
      estimate = difference, SE = SE, df = df, q.ratio = q,
      lower.CL = difference - half_width, upper.CL = difference + half_width,
      p.value = p, stringsAsFactors = FALSE)
  })
  table <- do.call(rbind, rows)
  rownames(table) <- NULL
  if (any(!is.finite(table$p.value))) computed$warnings <- c(computed$warnings,
    "一部のGames–Howell比較は分散を推定できない等の理由で検定不能です。")
  setNames(list(list(kind = "games_howell", raw = raw, formatted = table,
    intrinsic_adjustment = "Games–Howell", confidence_level = confidence_level,
    package_version = as.character(utils::packageVersion("PMCMRplus")),
    compared_factor = iv, source_effect = iv, test_id = "games_howell",
    interval_method = "Games-Howell studentized-range interval (approximate simultaneous coverage)",
    warnings = computed$warnings, messages = computed$messages, error = NULL)), iv)
}

t_welch <- function(experiment, alpha = 0.05, followup_alpha = 0.10) {
  if (!exists(".t_mean_prepare", mode = "function")) {
    stop("t_welch.Rはtest_helper.R経由で読み込んでください。", call. = FALSE)
  }
  results <- .t_mean_prepare(experiment, alpha, followup_alpha, "welch")
  results <- lapply(results, function(result) {
    if (result$n_levels == 2L) return(.t_mean_two(result, experiment$iv, welch = TRUE))
    samples <- result$samples
    model_data <- data.frame(response = result$analysis_data[[result$dependent_variable]],
      group = result$analysis_data[[experiment$iv]])
    computed <- .t_mean_capture(stats::oneway.test(response ~ group, data = model_data, var.equal = FALSE))
    test <- computed$value
    result$test <- test; result$test_id <- "welch_anova"
    result$variance_assumption <- "unequal group variances allowed"
    result$followup_test_id <- "games_howell"
    result$followup_factors <- experiment$iv
    result$adjust <- "Games–Howell"
    result$followup_settings <- list(test = "Games-Howell", contrast = "all pairs",
      confidence_level = 1 - alpha, additional_adjustment = "none")
    result$test_table <- data.frame(Effect = experiment$iv,
      df1 = if (is.null(test)) NA_real_ else unname(test$parameter[1L]),
      df2 = if (is.null(test)) NA_real_ else unname(test$parameter[2L]),
      F = if (is.null(test)) NA_real_ else unname(test$statistic),
      p.value = if (is.null(test)) NA_real_ else unname(test$p.value), stringsAsFactors = FALSE)
    result$anova_table <- result$test_table
    result$warnings <- computed$warnings; result$messages <- c(result$messages, computed$messages)
    result$error_details <- if (!is.null(computed$error)) list(stage = "Welch ANOVA", message = computed$error) else NULL
    if (!is.finite(result$test_table$p.value)) result$messages <- c(result$messages,
      "群内分散が0である等の理由でWelch ANOVAを算出できませんでした。", computed$error)
    iv <- experiment$iv; levels <- result$levels; confidence_level <- 1 - alpha
    # 全体検定を表示した後にhelperが呼び出す。
    result$compute_followups <- function(effects, forced = character(), main_effects = character()) {
      if (!iv %in% effects) return(list())
      .t_welch_followups(samples, levels, iv, confidence_level)
    }
    result
  })
  .t_mean_collection(results, experiment, "welch")
}
