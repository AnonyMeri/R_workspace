# ============================================================
# 順位検定の計算。方法の受付・下位検定の選択・補正・表示はhelperが担当。
#
# 1要因のみ：
#   between・2水準   stats::wilcox.test（順位和）
#   within・2水準    stats::wilcox.test（符号付順位）
#   between・3水準〜 stats::kruskal.test
#   within・3水準〜  stats::friedman.test
#
# 下位検定は未補正の全ペアを返す：
#   Kruskal–Wallis後 PMCMRplus::kwAllPairsDunnTest(p.adjust.method = "none")
#   Friedman後      IDで対応させたstats::wilcox.test(paired = TRUE)
# 共変量は使用しない。複数観測は参加者×水準の平均に集約する。
# 対応ありは全水準の揃った参加者を用い、下位検定にも同じ参加者を使う。
# Wilcoxonのexact = NULL / correct = TRUEはR標準の既定に従う。
# 信頼区間・元尺度の平均差・偏η²は計算しない。
# ============================================================

t_rank <- function(experiment, alpha = 0.05, followup_alpha = 0.10) {
  if (!is.list(experiment) || !is.data.frame(experiment$data)) {
    stop("experimentには前処理済みのデータを指定してください。", call. = FALSE)
  }
  scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!scalar(alpha) || !scalar(followup_alpha) || alpha <= 0 ||
      alpha > followup_alpha || followup_alpha >= 1) {
    stop("0 < alpha ≤ followup_alpha < 1 にしてください。", call. = FALSE)
  }
  id <- experiment$id
  iv <- experiment$iv
  dv <- experiment$dv
  valid_names <- function(x) is.character(x) && length(x) > 0L && !anyNA(x) && all(nzchar(x))
  if (!valid_names(id) || length(id) != 1L || !valid_names(iv) || length(iv) != 1L) {
    stop("ランクは1要因の比較に対応しています。id・ivはそれぞれ1列を指定してください。", call. = FALSE)
  }
  if (!valid_names(dv) || anyDuplicated(c(id, iv, dv))) {
    stop("id・iv・dvには重複しない列を指定してください。", call. = FALSE)
  }
  data <- experiment$data
  required <- c(id, iv, dv)
  missing_columns <- setdiff(required, names(data))
  if (length(missing_columns)) stop("必要な列がありません: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  if (anyDuplicated(names(data)) || !nrow(data) || anyNA(data[[id]]) || anyNA(data[[iv]])) {
    stop("データの列名・行数、ID・要因列の欠損を確認してください。", call. = FALSE)
  }
  paired <- identical(as.character(experiment$within), iv) && !length(experiment$between)
  independent <- identical(as.character(experiment$between), iv) && !length(experiment$within)
  if (!paired && !independent) stop("ランクには1要因のwithinまたはbetween計画が必要です。", call. = FALSE)
  data[[iv]] <- if (is.factor(data[[iv]])) droplevels(data[[iv]]) else factor(data[[iv]])
  level_names <- levels(data[[iv]])
  k <- length(level_names)
  if (k < 2L) stop("ランクには2水準以上が必要です。", call. = FALSE)
  if (independent && any(vapply(split(as.character(data[[iv]]), data[[id]], drop = TRUE),
      function(x) length(unique(x)) > 1L, logical(1)))) {
    stop("対応なしの計画ですが、同じIDが複数群に含まれています。", call. = FALSE)
  }
  test_id <- if (k == 2L) {
    if (paired) "wilcoxon_signed_rank" else "wilcoxon_rank_sum"
  } else if (paired) "friedman" else "kruskal_wallis"
  followup_id <- if (k == 2L) NULL else if (paired) "wilcoxon_signed_rank" else "dunn"

  # Rの標準警告（同順位・ゼロ差など）も結果に保持し、表示はhelperへ委譲する。
  calculate <- function(expression) {
    notices <- character()
    value <- withCallingHandlers(force(expression), warning = function(w) {
      notices <<- unique(c(notices, conditionMessage(w)))
      invokeRestart("muffleWarning")
    })
    list(value = value, warnings = notices)
  }
  wilcoxon <- function(x, y, paired) {
    calculate(stats::wilcox.test(x, y, paired = paired, alternative = "two.sided",
      exact = NULL, correct = TRUE, conf.int = FALSE))
  }
  results <- setNames(lapply(dv, function(target) {
    if (!is.numeric(data[[target]]) || any(is.infinite(data[[target]]))) {
      stop("従属変数は有限の数値またはNAにしてください: ", target, call. = FALSE)
    }
    source_ids <- unique(as.character(data[[id]]))
    used <- data[!is.na(data[[target]]), c(id, iv, target), drop = FALSE]
    if (!nrow(used)) stop("検定可能な観測がありません: ", target, call. = FALSE)
    aggregated <- stats::aggregate(used[[target]], by = used[c(id, iv)], FUN = mean)
    names(aggregated)[3L] <- target
    aggregated[[iv]] <- factor(aggregated[[iv]], levels = level_names)
    # 値を持たない水準を落として別の検定へ切り替えることはしない。
    if (any(table(aggregated[[iv]]) == 0L)) {
      stop("従属変数に有効な観測のない水準があります: ", target, call. = FALSE)
    }
    if (paired) {
      counts <- table(as.character(aggregated[[id]]))
      complete_ids <- names(counts)[counts == k]
      aggregated <- aggregated[as.character(aggregated[[id]]) %in% complete_ids, , drop = FALSE]
      if (length(complete_ids) < 2L) stop("全水準の揃った参加者が2名以上必要です: ", target, call. = FALSE)
    }
    aggregated <- aggregated[order(aggregated[[iv]], as.character(aggregated[[id]])), , drop = FALSE]
    rownames(aggregated) <- NULL
    sample_ids <- unique(as.character(aggregated[[id]]))
    excluded_ids <- setdiff(source_ids, sample_ids)
    samples <- lapply(level_names, function(level) {
      cell <- aggregated[aggregated[[iv]] == level, , drop = FALSE]
      if (paired) cell <- cell[match(sample_ids, as.character(cell[[id]])), , drop = FALSE]
      cell[[target]]
    })
    names(samples) <- level_names
    calculated <- if (k == 2L) wilcoxon(samples[[1L]], samples[[2L]], paired) else if (paired) {
      calculate(stats::friedman.test(do.call(cbind, samples)))
    } else calculate(stats::kruskal.test(aggregated[[target]], aggregated[[iv]]))
    test <- calculated$value
    statistic_name <- if (k == 2L) { if (paired) "V" else "W" } else if (paired) "χ²" else "H"
    table <- data.frame(Effect = iv, statistic_name = statistic_name,
      statistic = unname(test$statistic),
      df = if (length(test$parameter)) unname(test$parameter) else NA_real_,
      p.value = unname(test$p.value), stringsAsFactors = FALSE)
    notes <- character()
    if (length(experiment$covariate)) notes <- c(notes, "ランクでは共変量を使用していません。")
    if (nrow(used) > nrow(unique(used[c(id, iv)]))) {
      notes <- c(notes, "複数観測を参加者×水準の平均に集約しています。")
    }
    if (anyNA(data[[target]])) notes <- c(notes, paste0("欠損観測を除外: ", sum(is.na(data[[target]])), "行。"))
    if (length(excluded_ids)) notes <- c(notes, paste0("分析対象外の参加者: ", length(excluded_ids), "名。"))
    if (paired && length(excluded_ids)) notes <- c(notes, "全水準の揃った参加者のみで、全体検定と下位検定を行います。")
    if (!is.finite(test$p.value)) notes <- c(notes, "順位の変動がない等の理由で、p値を算出できませんでした。")

    # 各DVのデータを閉包に保持し、helperが要求した段階で全ペアを計算する。
    compute_followups <- if (k > 2L) function(effects, forced = character(), main_effects = character()) {
      if (!iv %in% effects) return(list())
      pairs <- utils::combn(seq_len(k), 2L)
      raw_tests <- list()
      warnings <- character()
      dunn <- NULL
      if (!paired) {
        if (!requireNamespace("PMCMRplus", quietly = TRUE)) {
          stop('Dunn検定にはPMCMRplusが必要です。install.packages("PMCMRplus")を実行してください。', call. = FALSE)
        }
        calculated <- calculate(PMCMRplus::kwAllPairsDunnTest(aggregated[[target]],
          aggregated[[iv]], p.adjust.method = "none"))
        dunn <- calculated$value
        warnings <- calculated$warnings
      }
      rows <- lapply(seq_len(ncol(pairs)), function(index) {
        left <- pairs[1L, index]
        right <- pairs[2L, index]
        a <- level_names[left]
        b <- level_names[right]
        if (paired) {
          calculated <- wilcoxon(samples[[left]], samples[[right]], TRUE)
          raw_tests[[index]] <<- calculated$value
          warnings <<- unique(c(warnings, calculated$warnings))
          statistic <- unname(calculated$value$statistic)
          p <- calculated$value$p.value
          inference <- calculated$value$method
        } else {
          statistic <- unname(dunn$statistic[b, a])
          p <- unname(dunn$p.value[b, a])
          inference <- dunn$method
        }
        data.frame(level1 = a, level2 = b, pair_key = paste(left, right, sep = ":"),
          comparison = paste(a, b, sep = " - "), statistic_name = if (paired) "V" else "z",
          statistic = statistic, p.value = unname(p), inference = inference,
          stringsAsFactors = FALSE)
      })
      unadjusted <- do.call(rbind, rows)
      setNames(list(list(kind = "rank", raw = if (paired) raw_tests else dunn,
        unadjusted = unadjusted, formatted = unadjusted, compared_factor = iv,
        source_effect = iv, test_id = followup_id, warnings = warnings, error = NULL)), iv)
    } else NULL
    result <- list(analysis = "rank", dependent_variable = target, test_id = test_id,
      followup_test_id = followup_id, design = if (paired) "within" else "between",
      within = experiment$within, between = experiment$between,
      covariate = experiment$covariate, covariate_adjusted = FALSE,
      paired = paired, levels = level_names, n_levels = k, n_participants = length(sample_ids),
      group_sizes = table(aggregated[[iv]]), excluded_ids = excluded_ids,
      missing_rows = sum(is.na(data[[target]])),
      aggregated_rows = nrow(used) - nrow(unique(used[c(id, iv)])),
      analysis_data = aggregated, test = test, test_table = table,
      followup_factors = if (k > 2L) iv else character(),
      compute_followups = compute_followups, contrasts = list(),
      alpha = alpha, followup_alpha = followup_alpha, adjust = "none",
      messages = notes, warnings = calculated$warnings)
    class(result) <- c("analysis_rank_result", "analysis_test_result", "list")
    result
  }), dv)
  attr(results, "stage") <- "test"
  attr(results, "type") <- "rank"
  attr(results, "experiment") <- experiment
  class(results) <- c("analysis_rank", "analysis_test", "list")
  invisible(results)
}
