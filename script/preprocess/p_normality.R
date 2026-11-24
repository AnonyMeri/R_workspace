# ============================================================
# 前処理：正規性検定（preprocessとは独立して実行）
#
# 公開関数：
#   p_normality()（内部）
#   normality()（preprocess_helper.Rから公開）
#
# p_prepare()が返すexperimentを受け取り、全独立変数を
# 組み合わせた各条件セルで正規性を確認します。
#
# 既定の検定：
#   Shapiro-Wilk検定
#
# 任意で追加できる検定：
#   Anderson-Darling検定、Lilliefors検定
#   追加検定にはnortestパッケージが必要です。
#
# 注意：この関数は分析前のスクリーニングです。ANOVAの正規性は
# 最終的にはモデル残差のQ-Qプロット等も含めて判断してください。
# 参考候補は主要検定（methodsの先頭）で判定します。
# 検定可能セルの半数以上が棄却ならART / 置換ANCOVA、それ未満なら
# 通常ANOVA / ANCOVAを参考候補とします。検定不可セルは分母から除きます。
# 棄却が一部にある場合や検定不可セルがある場合は注意を付記します。
# 結果には生の検定値、検定別集計、参考候補と判定根拠を保存します。
# ============================================================

p_normality <- function(
  experiment,
  methods = "shapiro",
  alpha = 0.05,
  digits = 3,
  p_digits = 3,
  show = TRUE
) {

  significance <- function(p) {
    ifelse(
      is.na(p), "",
      ifelse(
        p < 0.001, "***",
        ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ifelse(p < 0.10, ".", "")))
      )
    )
  }

  # 添付のtest_helperと同じ表示規則。数値の丸めは表示時だけ行う。
  format_p <- function(p) {
    threshold <- 10^(-p_digits)
    output <- rep("—", length(p))
    valid <- is.finite(p)
    output[valid] <- formatC(p[valid], format = "f", digits = p_digits)
    output[valid & p < threshold] <- paste0("<", formatC(threshold, format = "f", digits = p_digits))
    sub("^(<?)0\\.", "\\1.", output)
  }

  format_table <- function(table) {
    table <- as.data.frame(table, stringsAsFactors = FALSE)
    right <- intersect(c("p.value", "significance"), names(table))
    table <- table[, c(setdiff(names(table), right), right), drop = FALSE]
    for (name in names(table)) {
      values <- table[[name]]
      if (name == "p.value") {
        table[[name]] <- format_p(values)
      } else if (is.numeric(values)) {
        rendered <- rep("—", length(values))
        valid <- is.finite(values)
        rendered[valid] <- formatC(values[valid], format = "f",
          digits = if (name %in% c("n", "total", "testable", "rejected", "not_rejected", "untestable")) 0 else digits)
        table[[name]] <- rendered
      } else {
        table[[name]] <- as.character(values)
        table[[name]][is.na(values)] <- "—"
      }
    }
    labels <- c(method = "方法", n = "N", statistic = "統計量", judgement = "判定",
      p.value = "p", significance = "", total = "条件数", testable = "検定可能",
      rejected = "棄却", not_rejected = "棄却せず", untestable = "検定不可", rate = "棄却率")
    names(table) <- vapply(names(table), function(name) {
      if (name %in% names(labels)) labels[[name]] else name
    }, character(1))
    table
  }

  print_table <- function(table) {
    if (!nrow(table)) return(invisible(NULL))
    table <- format_table(table)
    widths <- vapply(seq_along(table), function(j) {
      max(nchar(c(names(table)[j], table[[j]]), type = "width"))
    }, numeric(1))
    old <- options(width = min(10000L, max(getOption("width"), sum(widths) + 3L * ncol(table) + 4L)))
    on.exit(options(old), add = TRUE)
    print(table, row.names = FALSE, right = TRUE, quote = FALSE)
    invisible(NULL)
  }

  heading <- function(section, labels = character()) {
    detail <- if (length(labels)) paste0("（", paste(labels, collapse = " | "), "）") else ""
    title <- paste0("──【", section, detail, "】")
    if (nchar(title, type = "width") > 38L) {
      title <- paste0("──【", section, "】")
      detail <- paste0("  ", detail)
    } else detail <- NULL
    cat("\n", title, strrep("─", max(2L, 40L - nchar(title, type = "width"))), "\n", sep = "")
    if (!is.null(detail)) cat(detail, "\n", sep = "")
  }

  normalize_method <- function(method) {
    aliases <- c(
      shapiro = "shapiro",
      anderson = "anderson",
      ad = "anderson",
      lilliefors = "lilliefors",
      lillie = "lilliefors"
    )

    key <- tolower(method)
    if (!key %in% names(aliases)) {
      stop(
        "methodsにはshapiro, anderson, lillieforsのいずれかを指定してください。",
        call. = FALSE
      )
    }

    unname(aliases[key])
  }

  run_test <- function(values, method) {
    x <- values[!is.na(values)]
    n <- length(x)
    unique_n <- length(unique(x))

    if (unique_n < 3) {
      return(list(
        statistic = NA_real_,
        p.value = NA_real_,
        judgement = "検定不可",
        note = "異なる値が3種類未満です。"
      ))
    }

    if (method == "shapiro") {
      if (n < 3 || n > 5000) {
        return(list(
          statistic = NA_real_,
          p.value = NA_real_,
          judgement = "検定不可",
          note = "Shapiro-Wilk検定はN=3～5000で実行します。"
        ))
      }

      test <- stats::shapiro.test(x)
    } else if (method == "anderson") {
      if (n < 8) {
        return(list(
          statistic = NA_real_,
          p.value = NA_real_,
          judgement = "検定不可",
          note = "Anderson-Darling検定には8件以上を必要とします。"
        ))
      }

      test <- nortest::ad.test(x)
    } else {
      if (n < 5) {
        return(list(
          statistic = NA_real_,
          p.value = NA_real_,
          judgement = "検定不可",
          note = "Lilliefors検定には5件以上を必要とします。"
        ))
      }

      test <- nortest::lillie.test(x)
    }

    p <- test$p.value

    list(
      statistic = unname(test$statistic),
      p.value = p,
      judgement = if (p < alpha) "正規性を棄却" else "正規性を棄却せず",
      note = ""
    )
  }

  participant_values <- function(data, id, dv, index) {
    ids <- data[[id]][index]
    values <- data[[dv]][index]
    positions <- split(seq_along(index), ids, drop = TRUE)

    vapply(
      positions,
      function(position) {
        x <- values[position]
        if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
      },
      numeric(1)
    )
  }

  required <- c("data", "id", "iv", "dv", "within", "between", "design")
  missing_elements <- setdiff(required, names(experiment))

  if (length(missing_elements) > 0) {
    stop(
      "experimentに次の要素がありません: ",
      paste(missing_elements, collapse = ", "),
      call. = FALSE
    )
  }

  data <- experiment$data
  id <- experiment$id
  iv <- experiment$iv
  dv <- experiment$dv
  within <- experiment$within
  between <- experiment$between
  covariate <- if (!is.null(experiment$covariate)) experiment$covariate else character(0)

  if (!is.data.frame(data)) {
    stop("experiment$dataはデータフレームである必要があります。", call. = FALSE)
  }
  if (!nrow(data)) stop("正規性検定に使用するデータがありません。", call. = FALSE)
  if (length(id) != 1) {
    stop("experiment$idには参加者ID列を1つ指定してください。", call. = FALSE)
  }
  if (length(iv) == 0) {
    stop("正規性検定には少なくとも1つの独立変数が必要です。", call. = FALSE)
  }
  if (length(dv) == 0) {
    stop("従属変数が指定されていません。", call. = FALSE)
  }
  if (!all(within %in% iv) || !all(between %in% iv)) {
    stop(
      "withinとbetweenにはivに含まれる要因を指定してください。",
      call. = FALSE
    )
  }
  if (length(intersect(within, between)) > 0) {
    stop(
      "同じ要因をwithinとbetweenの両方に指定できません。",
      call. = FALSE
    )
  }

  required_columns <- unique(c(id, iv, dv))
  missing_columns <- setdiff(required_columns, names(data))

  if (length(missing_columns) > 0) {
    stop(
      "dataに次の列がありません: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  non_numeric <- dv[!vapply(data[dv], is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop(
      "次の従属変数が数値型ではありません: ",
      paste(non_numeric, collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.character(methods) || length(methods) == 0 || anyNA(methods)) {
    stop("methodsには正規性検定名を文字列で指定してください。", call. = FALSE)
  }

  if (any(tolower(methods) == "all")) {
    methods <- c("shapiro", "anderson", "lilliefors")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alphaには0より大きく1より小さい数値を指定してください。", call. = FALSE)
  }
  if (!is.numeric(digits) || length(digits) != 1 || !is.finite(digits) ||
      digits < 0 || digits != floor(digits)) {
    stop("digitsには0以上の整数を指定してください。", call. = FALSE)
  }
  if (!is.numeric(p_digits) || length(p_digits) != 1 || !is.finite(p_digits) ||
      p_digits < 1 || p_digits != floor(p_digits)) {
    stop("p_digitsには1以上の整数を指定してください。", call. = FALSE)
  }
  if (!is.logical(show) || length(show) != 1 || is.na(show)) {
    stop("showにはTRUEまたはFALSEを指定してください。", call. = FALSE)
  }

  methods <- unique(vapply(methods, normalize_method, character(1)))

  if (any(methods != "shapiro") && !requireNamespace("nortest", quietly = TRUE)) {
    stop(
      "Anderson-Darling検定またはLilliefors検定にはnortestパッケージが必要です。\n",
      'install.packages("nortest")を実行してください。',
      call. = FALSE
    )
  }

  method_labels <- c(
    shapiro = "Shapiro-Wilk",
    anderson = "Anderson-Darling",
    lilliefors = "Lilliefors"
  )

  groups <- do.call(
    interaction,
    c(
      data[iv],
      list(
        drop = TRUE,
        lex.order = TRUE,
        sep = "\r"
      )
    )
  )
  rows <- split(seq_len(nrow(data)), groups, drop = TRUE)
  results <- list()

  for (dv_name in dv) {
    test_rows <- list()
    row_index <- 1L

    for (index in rows) {
      levels <- data[index[1], iv, drop = FALSE]
      values <- participant_values(data, id, dv_name, index)
      n <- sum(!is.na(values))

      for (method in methods) {
        test <- run_test(values, method)

        test_rows[[row_index]] <- cbind(
          levels,
          data.frame(
            method = method_labels[method],
            n = as.integer(n),
            statistic = test$statistic,
            p.value = test$p.value,
            significance = significance(test$p.value),
            judgement = test$judgement,
            note = test$note,
            stringsAsFactors = FALSE
          )
        )

        row_index <- row_index + 1L
      }
    }

    tests <- do.call(rbind, test_rows)
    rownames(tests) <- NULL

    primary_label <- unname(method_labels[methods[1]])
    summary_rows <- lapply(methods, function(method) {
      selected <- tests[tests$method == method_labels[method], , drop = FALSE]
      testable <- is.finite(selected$p.value)
      tested_n <- sum(testable)
      rejected_n <- sum(selected$p.value[testable] < alpha)
      data.frame(method = unname(method_labels[method]), total = nrow(selected),
        testable = tested_n, rejected = rejected_n, not_rejected = tested_n - rejected_n,
        untestable = sum(!testable),
        rejection_rate = if (tested_n) rejected_n / tested_n else NA_real_,
        stringsAsFactors = FALSE)
    })
    normality_summary <- do.call(rbind, summary_rows)
    rownames(normality_summary) <- NULL
    primary <- normality_summary[1L, , drop = FALSE]
    tested_n <- primary$testable
    rejected_n <- primary$rejected
    untestable_n <- primary$untestable
    # ユーザー指定の基準：主要検定の検定可能セルの半数以上。
    rejection_threshold <- 0.5
    covariate_present <- length(covariate) > 0L
    has_many_rejections <- tested_n > 0L && primary$rejection_rate >= rejection_threshold
    status <- if (!tested_n) "untestable" else if (has_many_rejections) "many_rejections" else
      if (rejected_n > 0L) "some_rejections" else "no_rejection"
    recommended_method <- if (!tested_n) NA_character_ else
      if (has_many_rejections) {
        if (covariate_present) "permuco_ancova" else "art"
      } else if (covariate_present) "ancova" else "anova"
    label <- if (is.na(recommended_method)) "判定保留" else switch(recommended_method,
      anova = "通常ANOVA", ancova = "通常ANCOVA",
      art = "ART", permuco_ancova = "置換ANCOVA（permuco）")
    if (!tested_n) {
      reason <- "主要検定で検定可能な条件がないため、分析方法の参考候補を提示できません。"
    } else {
      relation <- if (has_many_rejections) "半数以上" else "半数未満"
      reason <- paste0("主要検定で", rejected_n, " / ", tested_n,
        "条件の正規性が棄却され、検定可能条件の", relation, "です。共変量は",
        if (covariate_present) "指定されています。" else "指定されていません。")
    }
    cautions <- character()
    if (rejected_n > 0L && !has_many_rejections) {
      cautions <- c(cautions, "半数未満ですが、正規性が棄却された条件があります。該当条件の分布も確認してください。")
    }
    if (untestable_n > 0L) {
      cautions <- c(cautions, paste0("検定不可の", untestable_n,
        "条件は棄却率の分母から除外しています。未検定の条件を正規性ありとは扱いません。"))
    }
    cautions <- c(cautions,
      "正規性が棄却されないことは正規分布の証明ではありません。半数という基準は今回の参考表示用の運用基準です。",
      "正規性検定だけで分析方法を確定せず、分布・外れ値・標本数を確認し、ANOVA・ANCOVAではモデル推定後の残差も確認してください。")
    recommendation <- list(
      method = label, analysis_method = recommended_method,
      candidate = if (is.na(recommended_method)) "分布と標本数を確認" else
        paste0('analyze(experiment, method = "', recommended_method, '")'),
      reason = reason, caution = paste(cautions, collapse = "\n"),
      status = status, primary_method = primary_label,
      testable_cells = tested_n, rejected_cells = rejected_n, untestable_cells = untestable_n,
      rejection_rate = primary$rejection_rate, rejection_threshold = rejection_threshold)

    normality_data <- data[
      ,
      unique(c(id, iv, covariate, dv_name)),
      drop = FALSE
    ]

    results[[dv_name]] <- list(
      analysis = "normality",
      dependent_variable = dv_name,
      design = experiment$design,
      iv = iv,
      within = within,
      between = between,
      covariate = covariate,
      analysis_data = normality_data,
      tests = tests,
      normality_summary = normality_summary,
      primary_method = primary_label,
      recommendation = recommendation,
      alpha = alpha
    )

    class(results[[dv_name]]) <- c(
      "analysis_normality_result",
      "analysis_preprocess_result",
      "analysis_test_result",
      "list"
    )

    if (show) {
      cat("\n========================================\n", dv_name,
        "\n========================================\n", sep = "")
      heading("正規性検定")
      cat("集計単位：条件セル内の参加者平均\n")
      cat("棄却基準：p < ", sub("^0\\.", ".", format(alpha, scientific = FALSE, trim = TRUE)), "\n", sep = "")
      for (index in seq_along(methods)) {
        method <- methods[index]
        selected <- tests[tests$method == method_labels[method], , drop = FALSE]
        cat("\n", unname(method_labels[method]), if (index == 1L) "（主要検定）" else "", "\n", sep = "")
        displayed <- selected[, c(iv, "n", "statistic", "judgement", "p.value", "significance"), drop = FALSE]
        statistic_label <- switch(method, shapiro = "W", anderson = "A", lilliefors = "D")
        names(displayed)[names(displayed) == "statistic"] <- statistic_label
        print_table(displayed)
        notes <- which(nzchar(selected$note))
        if (length(notes)) {
          text <- vapply(notes, function(j) paste0(
            paste(paste0(iv, "=", vapply(iv, function(name) as.character(selected[[name]][j]), character(1))), collapse = ", "), "：", selected$note[j]), character(1))
          cat("\n注記：", paste(text, collapse = "\n      "), "\n", sep = "")
        }
      }
      heading("分析方法の参考")
      displayed <- normality_summary[, c("method", "testable", "rejected", "untestable"), drop = FALSE]
      displayed$rate <- ifelse(is.finite(normality_summary$rejection_rate),
        paste0(formatC(100 * normality_summary$rejection_rate, format = "f", digits = 1), "%"), "—")
      print_table(displayed)
      cat("\n主要検定：", primary_label,
        "\n切替基準：検定可能な条件の半数以上で棄却（50%以上）",
        "\n参考候補：", recommendation$method,
        if (is.na(recommended_method)) "\n確認事項：" else "\n実行例：", recommendation$candidate,
        "\n理由：", recommendation$reason, "\n", sep = "")
      cat("\n注記：", paste(cautions, collapse = "\n      "), "\n", sep = "")
    }
  }

  attr(results, "stage") <- "preprocess"
  attr(results, "type") <- "normality"
  attr(results, "methods") <- methods
  attr(results, "experiment") <- experiment
  class(results) <- c("analysis_normality", "analysis_preprocess", "analysis_test", "list")
  invisible(results)
}
