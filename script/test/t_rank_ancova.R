# ============================================================
# 順位ANCOVA：npANCOVA::Quade（対応なし・1要因・共変量1列以上）。
# 順位化した応答を順位化した共変量で調整し、残差の群間差を検定。
# 下位検定はnpANCOVA::pairwise_npANCOVA(method="Quade")と同じ
# 各2群での再順位化・再調整。未補正の全ペアをhelperへ返す。
# 群を含まない補助回帰の係数・p値は、群効果の表とは別に保存する。
# 表示・実行基準・補正管理はtest_helper.Rが担当。
# Quade (1967), doi:10.1080/01621459.1967.10500925
# https://cran.r-project.org/web/packages/npANCOVA/refman/npANCOVA.html
# ============================================================

.t_rank_ancova_design <- function(experiment) {
  valid <- function(x) is.character(x) && length(x) > 0L && !anyNA(x) && all(nzchar(x))
  if (!is.list(experiment) || !is.data.frame(experiment$data) ||
      !valid(experiment$id) || length(experiment$id) != 1L ||
      !valid(experiment$iv) || length(experiment$iv) != 1L || !valid(experiment$dv)) {
    stop("順位ANCOVAには前処理済みのexperimentと、ID・1要因・従属変数を指定してください。", call. = FALSE)
  }
  if (length(experiment$within) || !identical(as.character(experiment$between), experiment$iv)) {
    stop("順位ANCOVA（Quade）は対応なし・1要因に対応しています。", call. = FALSE)
  }
  covariate <- experiment$covariate
  if (!valid(covariate)) stop("順位ANCOVAには共変量を1列以上指定してください。", call. = FALSE)
  columns <- c(experiment$id, experiment$iv, experiment$dv, covariate)
  data <- experiment$data
  if (anyDuplicated(columns) || anyDuplicated(names(data)) || !all(columns %in% names(data)) || !nrow(data)) {
    stop("順位ANCOVAの列名・行数を確認してください。ID・要因・DV・共変量は重複できません。", call. = FALSE)
  }
  if (anyNA(data[[experiment$id]]) || anyNA(data[[experiment$iv]])) {
    stop("ID・要因列の欠損を確認してください。", call. = FALSE)
  }
  group <- if (is.factor(data[[experiment$iv]])) droplevels(data[[experiment$iv]]) else factor(data[[experiment$iv]])
  if (nlevels(group) < 2L) stop("順位ANCOVAには2群以上が必要です。", call. = FALSE)
  if (any(vapply(split(as.character(group), data[[experiment$id]], drop = TRUE),
      function(x) length(unique(x)) > 1L, logical(1)))) {
    stop("対応なしの計画ですが、同じIDが複数群に含まれています。", call. = FALSE)
  }
  for (name in covariate) {
    values <- data[[name]]
    if (!is.numeric(values) || any(is.infinite(values))) {
      stop("共変量は有限の数値またはNAにしてください: ", name, call. = FALSE)
    }
    counts <- vapply(split(values, data[[experiment$id]], drop = TRUE),
      function(x) length(unique(x[!is.na(x)])), integer(1))
    if (any(counts > 1L)) stop("共変量は参加者内で一定の値を指定してください: ", name, call. = FALSE)
  }
  invisible(levels(group))
}

.t_rank_ancova_capture <- function(expression) {
  warnings <- messages <- character(); error <- NULL
  value <- tryCatch(withCallingHandlers(force(expression), warning = function(w) {
    warnings <<- unique(c(warnings, conditionMessage(w))); invokeRestart("muffleWarning")
  }, message = function(m) {
    messages <<- unique(c(messages, conditionMessage(m))); invokeRestart("muffleMessage")
  }), error = function(e) { error <<- conditionMessage(e); NULL })
  list(value = value, warnings = warnings, messages = messages, error = error)
}

.t_rank_ancova_fit <- function(formula, data) {
  # 定数・線形従属の共変量を黙って除外しない。ペア比較にも同じ検証を行う。
  covariates <- setdiff(names(data), c("response", "group"))
  ranked_cov <- as.data.frame(lapply(data[covariates], rank, ties.method = "average"))
  if (any(vapply(ranked_cov, function(x) length(unique(x)) < 2L, logical(1)))) {
    stop("分析対象内で一定値の共変量があります。", call. = FALSE)
  }
  cov_design <- cbind(1, as.matrix(ranked_cov))
  if (nrow(data) <= ncol(cov_design) || qr(cov_design)$rank < ncol(cov_design)) {
    stop("共変量の順位が線形従属しているか、補助回帰の残差自由度が不足しています。", call. = FALSE)
  }
  joint_design <- cbind(cov_design, stats::model.matrix(~ group, data)[, -1L, drop = FALSE])
  if (qr(joint_design)$rank < ncol(joint_design)) {
    stop("要因と順位化した共変量が完全に交絡しています。", call. = FALSE)
  }
  npANCOVA::Quade(formula = formula, data = data)
}

.t_rank_ancova_table <- function(test, effect) {
  if (is.null(test)) return(data.frame(Effect = effect, df1 = NA_real_, df2 = NA_real_, F = NA_real_, p.value = NA_real_))
  table <- test$anova[[1L]]
  result <- data.frame(Effect = effect, df1 = unname(table[1L, "Df"]),
    df2 = unname(table[nrow(table), "Df"]), F = unname(table[1L, "F value"]),
    p.value = unname(table[1L, "Pr(>F)"]))
  # 完全に調整し尽くされた応答に残る丸め誤差を、群の効果として検定しない。
  y <- test$data$ranked_response; residuals <- test$data$residuals
  scale <- max(1, sum((y - mean(y))^2))
  if (sum(residuals^2) <= .Machine$double.eps * scale) {
    result$F <- result$p.value <- NA_real_
  }
  result
}

t_rank_ancova <- function(experiment, alpha = 0.05, followup_alpha = 0.10,
                          covariate_interactions = FALSE) {
  level_names <- .t_rank_ancova_design(experiment)
  scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!scalar(alpha) || !scalar(followup_alpha) || alpha <= 0 || alpha > followup_alpha || followup_alpha >= 1) {
    stop("0 < alpha ≤ followup_alpha < 1 にしてください。", call. = FALSE)
  }
  if (!identical(covariate_interactions, FALSE)) {
    stop("順位ANCOVA（Quade）は共変量との交互作用を扱いません。covariate_interactions = FALSEを指定してください。", call. = FALSE)
  }
  if (!requireNamespace("npANCOVA", quietly = TRUE) || utils::packageVersion("npANCOVA") < "0.2.0") {
    stop('順位ANCOVAにはnpANCOVA 0.2.0以上が必要です。install.packages("npANCOVA")を実行してください。', call. = FALSE)
  }
  data <- experiment$data; id <- experiment$id; iv <- experiment$iv; covariate <- experiment$covariate
  results <- setNames(lapply(experiment$dv, function(target) {
    if (!is.numeric(data[[target]]) || any(is.infinite(data[[target]]))) {
      stop("従属変数は有限の数値またはNAにしてください: ", target, call. = FALSE)
    }
    columns <- c(id, iv, target, covariate)
    complete <- stats::complete.cases(data[columns])
    used <- data[complete, columns, drop = FALSE]
    if (!nrow(used)) stop("順位ANCOVAに使用できる完全ケースがありません: ", target, call. = FALSE)
    aggregated <- stats::aggregate(used[[target]], by = used[c(id, iv)], FUN = mean)
    names(aggregated)[3L] <- target
    aggregated[[iv]] <- factor(aggregated[[iv]], levels = level_names)
    if (any(table(aggregated[[iv]]) < 2L)) stop("各群に有効な参加者が2名以上必要です: ", target, call. = FALSE)
    cov_rows <- match(as.character(aggregated[[id]]), as.character(used[[id]]))
    for (name in covariate) aggregated[[name]] <- used[[name]][cov_rows]
    aggregated <- aggregated[order(aggregated[[iv]], as.character(aggregated[[id]])), , drop = FALSE]
    rownames(aggregated) <- NULL
    model_data <- data.frame(response = aggregated[[target]])
    cov_names <- paste0("cv", seq_along(covariate))
    for (i in seq_along(covariate)) model_data[[cov_names[i]]] <- aggregated[[covariate[i]]]
    model_data$group <- aggregated[[iv]]
    formula <- stats::reformulate(c(cov_names, "group"), response = "response")
    # 構造上の問題は停止し、統計量の退化は結果に記録する。
    computed <- .t_rank_ancova_capture(.t_rank_ancova_fit(formula, model_data))
    if (!is.null(computed$error)) stop(computed$error, " [", target, "]", call. = FALSE)
    test <- computed$value
    table <- .t_rank_ancova_table(test, iv)
    coefficients <- test$regression_equation$coefficients
    rows <- match(paste0("ranked_", cov_names), rownames(coefficients))
    cov_table <- data.frame(Effect = covariate, estimate = coefficients[rows, "Estimate"],
      SE = coefficients[rows, "Std. Error"], df = test$regression_equation$df[2L],
      t.ratio = coefficients[rows, "t value"], p.value = coefficients[rows, "Pr(>|t|)"], row.names = NULL)
    sample_ids <- unique(as.character(aggregated[[id]]))
    excluded_ids <- setdiff(unique(as.character(data[[id]])), sample_ids)
    repeated <- nrow(used) - nrow(unique(used[c(id, iv)]))
    notes <- c(if (repeated) "複数観測を参加者×水準の平均に集約しています。",
      if (sum(!complete)) paste0("DV・共変量の欠損により除外: ", sum(!complete), "行。"),
      if (length(excluded_ids)) paste0("分析対象外の参加者: ", length(excluded_ids), "名。"),
      if (!is.finite(table$p.value)) "応答が共変量で完全に説明される等の理由で、群効果を検定できませんでした。")
    compute_followups <- if (length(level_names) > 2L) function(effects, forced = character(), main_effects = character()) {
      if (!iv %in% effects) return(list())
      pairs <- utils::combn(seq_along(level_names), 2L)
      raw_tests <- vector("list", ncol(pairs)); warnings <- character(); errors <- character()
      rows <- lapply(seq_len(ncol(pairs)), function(index) {
        left <- pairs[1L, index]; right <- pairs[2L, index]
        subset <- model_data[model_data$group %in% level_names[c(left, right)], , drop = FALSE]
        subset$group <- droplevels(subset$group)
        computed <- .t_rank_ancova_capture(.t_rank_ancova_fit(formula, subset))
        raw_tests[index] <<- list(computed$value)
        warnings <<- unique(c(warnings, computed$warnings, computed$messages))
        if (!is.null(computed$error)) {
          errors[index] <<- paste(level_names[left], level_names[right], computed$error, sep = " / ")
          warnings <<- unique(c(warnings, errors[index]))
        }
        pair <- .t_rank_ancova_table(computed$value, iv)
        if (is.null(computed$error) && !is.finite(pair$p.value)) {
          warnings <<- unique(c(warnings, paste0(level_names[left], " - ", level_names[right], "：群効果を検定できませんでした。")))
        }
        data.frame(level1 = level_names[left], level2 = level_names[right],
          pair_key = paste(left, right, sep = ":"), comparison = paste(level_names[left], level_names[right], sep = " - "),
          statistic_name = "F", statistic = pair$F, df1 = pair$df1, df2 = pair$df2,
          p.value = pair$p.value, stringsAsFactors = FALSE)
      })
      names(raw_tests) <- vapply(rows, `[[`, character(1), "pair_key")
      table <- do.call(rbind, rows)
      setNames(list(list(kind = "rank_ancova", raw = raw_tests, unadjusted = table, formatted = table,
        compared_factor = iv, source_effect = iv, test_id = "quade_pairwise",
        strategy = "Refit Quade in every two-group subset; rerank within each subset",
        reference = "npANCOVA::pairwise_npANCOVA(method = 'Quade')",
        warnings = warnings, pair_errors = errors[nzchar(errors) & !is.na(errors)], error = NULL)), iv)
    } else NULL
    result <- list(analysis = "rank_ancova", test_id = "quade", dependent_variable = target,
      design = "between", paired = FALSE, within = character(), between = iv,
      covariate = covariate, covariate_adjusted = TRUE, covariate_interactions = FALSE,
      covariate_table = cov_table, covariate_regression = test$regression_equation,
      covariate_inference = "Auxiliary regression of ranked response on ranked covariates, without group",
      levels = level_names, n_levels = length(level_names), n_participants = length(sample_ids),
      group_sizes = table(aggregated[[iv]]), excluded_ids = excluded_ids,
      missing_rows = sum(!complete), missing_dv_rows = sum(is.na(data[[target]])),
      missing_covariate_rows = sum(!stats::complete.cases(data[covariate])), aggregated_rows = repeated,
      analysis_data = aggregated, model_data = model_data, ranked_data = test$data,
      column_map = data.frame(role = c("dv", rep("covariate", length(covariate)), "iv"),
        original = c(target, covariate, iv), internal = c("response", cov_names, "group")),
      test = test, test_table = table, anova_table = table,
      followup_test_id = if (length(level_names) > 2L) "quade_pairwise" else NULL,
      followup_factors = iv, compute_followups = compute_followups, contrasts = list(),
      alpha = alpha, followup_alpha = followup_alpha, adjust = "none",
      package_version = as.character(utils::packageVersion("npANCOVA")),
      assumptions = c("Independent groups", "Common marginal distribution of covariates across populations"),
      assumptions_automatically_verified = FALSE,
      inference_status = if (is.finite(table$p.value)) "ok" else "not_estimable",
      error_details = if (is.finite(table$p.value)) NULL else list(stage = "Quade",
        message = "応答の調整後変動がない等の理由で、群効果の検定は定義できません。"),
      messages = notes, warnings = computed$warnings)
    class(result) <- c("analysis_rank_ancova_result", "analysis_test_result", "list")
    result
  }), experiment$dv)
  attr(results, "stage") <- "test"; attr(results, "type") <- "rank_ancova"; attr(results, "experiment") <- experiment
  class(results) <- c("analysis_rank_ancova", "analysis_test", "list")
  invisible(results)
}
