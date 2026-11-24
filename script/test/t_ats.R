# ============================================================
# 対応なし・1要因・3群以上：rankFDのATS（F近似）＋MCTP。
# t_brunner_munzel.Rの内部データ検証関数を使用。helperが先に読み込む。
# H0p・unweighted pseudo-ranks。MCTPは全ペア・Fisher近似。
# p値と同時信頼区間はrankFD内で補正済み。追加補正しない。
# Konietschke et al. (2023), https://doi.org/10.32614/RJ-2023-029
# ============================================================

# MCTPの多変量t積分は乱数を用いる。明示seed時は呼出元の乱数状態を復元。
.t_relative_rng <- function(seed, expression) {
  if (is.null(seed)) return(force(expression))
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) previous <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) assign(".Random.seed", previous, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(list = ".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  force(expression)
}

t_ats <- function(experiment, alpha = 0.05, followup_alpha = 0.10, seed = NULL) {
  if (!exists(".t_relative_prepare", mode = "function")) {
    stop("t_ats.Rはtest_helper.R経由で読み込んでください。", call. = FALSE)
  }
  prepared <- .t_relative_prepare(experiment, alpha, followup_alpha, "ats")
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != floor(seed))) {
    stop("seedには0以上の整数を指定してください。", call. = FALSE)
  }
  if (!requireNamespace("rankFD", quietly = TRUE)) {
    stop('ATS＋MCTPにはrankFDが必要です。install.packages("rankFD")を実行してください。', call. = FALSE)
  }
  results <- lapply(prepared, function(result) {
    iv <- experiment$iv; target <- result$dependent_variable
    # 日本語・空白・記号を含む列名も式の解釈に影響しないよう内部名で計算。
    model_data <- data.frame(response = result$analysis_data[[target]],
      group = result$analysis_data[[iv]])
    fit <- function(contrast = NULL) rankFD::rankFD(response ~ group, data = model_data,
      alpha = alpha, CI.method = "logit", effect = "unweighted", hypothesis = "H0p",
      contrast = contrast, sci.method = "fisher", info = FALSE, covariance = TRUE, rounds = Inf)
    computed <- .t_relative_capture(fit())
    test <- computed$value
    result$test <- test
    result$test_id <- "ats"
    result$followup_test_id <- "mctp"
    result$package_version <- as.character(utils::packageVersion("rankFD"))
    result$inference <- "ATS F approximation"
    result$hypothesis <- "H0p"
    result$effect <- "unweighted"
    result$estimate_definition <- "p_i = mean_j[P(X_j < X_i) + 0.5 P(X_j = X_i)]"
    result$seed <- seed
    result$adjust <- "MCTP (Fisher, Tukey)"
    result$followup_settings <- list(contrast = "Tukey all pairs", sci.method = "fisher",
      additional_adjustment = "none", confidence_level = 1 - alpha, seed = seed)
    ats <- if (is.null(test)) rep(NA_real_, 4L) else unname(test$ANOVA.Type.Statistic[1L, ])
    result$test_table <- data.frame(Effect = iv, statistic_name = "ATS", statistic = ats[1L],
      df1 = ats[2L], df2 = ats[3L], p.value = ats[4L], stringsAsFactors = FALSE)
    result$descriptive <- if (is.null(test)) NULL else test$Descriptive
    result$warnings <- computed$warnings
    result$messages <- c(result$messages, computed$messages)
    result$error_details <- if (!is.null(computed$error)) list(stage = "ATS", message = computed$error) else NULL
    if (!is.finite(ats[4L])) result$messages <- c(result$messages,
      "分散を推定できない等の理由でATSを算出できませんでした。", computed$error)
    result$followup_factors <- iv
    level_names <- result$levels
    # 全体検定の表示が済み、helperがこの要因を選んだ時点でMCTPを計算する。
    result$compute_followups <- function(effects, forced = character(), main_effects = character()) {
      if (!iv %in% effects) return(list())
      pairs <- utils::combn(seq_along(level_names), 2L)
      matrix <- matrix(0, nrow = ncol(pairs), ncol = length(level_names))
      for (i in seq_len(nrow(matrix))) {
        matrix[i, pairs[1L, i]] <- 1
        matrix[i, pairs[2L, i]] <- -1
      }
      computed <- .t_relative_capture(.t_relative_rng(seed, fit(list("group", matrix))))
      if (!is.null(computed$error)) {
        return(setNames(list(list(kind = "mctp", error = computed$error,
          error_details = list(stage = "MCTP", message = computed$error),
          warnings = computed$warnings, messages = computed$messages,
          source_effect = iv)), iv))
      }
      raw <- computed$value
      local <- raw$MCTP$Local.Results
      table <- data.frame(level1 = level_names[pairs[1L, ]], level2 = level_names[pairs[2L, ]],
        pair_key = paste(pairs[1L, ], pairs[2L, ], sep = ":"),
        comparison = paste(level_names[pairs[1L, ]], level_names[pairs[2L, ]], sep = " - "),
        estimate = local$Effect, SE = local$Std.Error, df = raw$MCTP$DF,
        t.ratio = local$T, lower.CL = local$Lower, upper.CL = local$Upper,
        p.value = local$p.value, stringsAsFactors = FALSE)
      if (any(!is.finite(table$p.value))) computed$warnings <- c(computed$warnings,
        "一部のMCTP比較は分散を推定できない等の理由で検定不能です。")
      setNames(list(list(kind = "mctp", raw = raw, formatted = table,
        intrinsic_adjustment = "MCTP (Fisher, Tukey)", compared_factor = iv,
        source_effect = iv, test_id = "mctp", confidence_level = 1 - alpha,
        warnings = computed$warnings, messages = computed$messages, error = NULL)), iv)
    }
    result
  })
  .t_relative_collection(results, experiment, "ats")
}
