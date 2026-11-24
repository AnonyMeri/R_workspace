# ============================================================
# 関係性分析の統合インターフェース（相関・回帰・SEM）
# 公開関数：relate()。計算処理はr_correlation.R / r_regression.R / r_sem.Rへ委譲する。
# 4ファイルを同じフォルダに置く。相関はR標準機能、SEMはlavaanを使う。
# 傾き分析にはemmeans、置換検定にはpermucoを使う。
#
# source("script/relation/relation_helper.R", encoding = "UTF-8")
# result <- relate(data)                       # GUIで変数・設定を選択
# result <- relate(data, analysis = "correlation")  # 相関分析のGUI
# result <- relate(data, variables = c("X", "Y", "Z"), method = "spearman")
# result <- relate(data, analysis = "sem")  # SEMのモデル式・設定をGUIで指定
# result <- relate(data, sem_model = "Y ~ X + Z", estimator = "MLR")
# result <- relate(data, formula = Y ~ X + Z)  # 引数で直接実行
# result <- relate(data, formula = Y ~ A * Z, trend = "Z", by = "A")
# result <- relate(model = fitted_model, trend = "Z", by = c("A", "B"))
# followup <- relate(result)  # analyze()の結果から有意な共変量との交互作用を抽出
# followup <- relate(result, select = "gui")  # 抽出した対象と設定をGUIで確認
# 既存のafex_aovモデルも同じ入口で利用する。再推定せず元の反復測定構造を使う。
# 結果入力では単一結果・従属変数別・方法別の入れ子を受け付ける。
# ANCOVAと置換ANCOVAはtest側で排他選択し、選択済みの方法で下位分析する。
# 過去の結果に両方が含まれる場合は、relate(result$ancova)などで一方を渡す。
# 通常ANCOVA・置換ANCOVAとも、有意な交互作用の条件別傾き・傾き差を検定する。
# 順位ANCOVA（Quade）の結果入力は拒否する。データからの分析にはこの制限を適用しない。
# 置換回数とseedは親結果から引き継ぐ。対応なしでは親の置換法を使う。
# 反復測定では参加者ごとの条件平均・条件差を応答とするFreedman–Lane検定を使う。
# 置換分析のSE・信頼区間は線形モデルに基づき、p値は置換検定に基づく。
# 計算の入口はhelperで選択。通常回帰と置換回帰へ確定済みの引数を渡す。
# 相関は選択した全ペアを1つの補正対象とする。反復測定相関・偏相関は今後追加する。
# SEMは単一集団の個票データを対象に、パス・潜在変数・共分散・定義した効果を扱う。
# SEMのp値補正は既定でnone。固定係数を除く検定可能な全パラメータを1つの補正対象とする。
# ============================================================

.relation_source <- NULL
for (.relation_frame in rev(sys.frames())) {
  if (!is.null(.relation_frame$ofile) && nzchar(.relation_frame$ofile)) {
    .relation_source <- .relation_frame$ofile
    break
  }
}
if (is.null(.relation_source)) stop("relation_helper.Rをsource()で読み込んでください。", call. = FALSE)
.relation_directory <- dirname(normalizePath(.relation_source, winslash = "/", mustWork = TRUE))
.relation_child <- file.path(.relation_directory, c("r_regression.R", "r_correlation.R", "r_sem.R"))
if (any(!file.exists(.relation_child))) stop("relation_helper.R・r_regression.R・r_correlation.R・r_sem.Rを同じフォルダに配置してください。", call. = FALSE)
.relation_private_environment <- new.env(parent = environment())
for (.relation_script in .relation_child) source(.relation_script, local = .relation_private_environment, encoding = "UTF-8")

evalq({
  validate_options <- function(alpha, followup_alpha, adjust, ci, digits, p_digits, show, compare_slopes) {
    scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
    if (!scalar(alpha) || !scalar(followup_alpha) || alpha <= 0 || alpha > followup_alpha || followup_alpha >= 1) {
      stop("0 < alpha ≤ followup_alpha < 1 にしてください。", call. = FALSE)
    }
    if (!is.character(adjust) || length(adjust) != 1L || is.na(adjust) || !adjust %in% stats::p.adjust.methods) {
      stop("adjustにはp.adjust()で利用できる補正法を指定してください。", call. = FALSE)
    }
    if (!scalar(ci) || ci <= 0 || ci >= 1) stop("ciは0より大きく1より小さくしてください。", call. = FALSE)
    if (!scalar(digits) || digits < 0 || digits != floor(digits) ||
        !scalar(p_digits) || p_digits < 1 || p_digits != floor(p_digits)) {
      stop("digitsは0以上、p_digitsは1以上の整数にしてください。", call. = FALSE)
    }
    for (flag in list(show, compare_slopes)) if (!is.logical(flag) || length(flag) != 1L || is.na(flag)) {
      stop("showとcompare_slopesにはTRUEまたはFALSEを指定してください。", call. = FALSE)
    }
    invisible(NULL)
  }
  make_formula <- function(dv, predictors, interactions = FALSE) {
    if (!is.character(dv) || length(dv) != 1L || is.na(dv) ||
        !is.character(predictors) || !length(predictors) || anyNA(predictors) || anyDuplicated(predictors) || dv %in% predictors) {
      stop("目的変数を1つ、異なる説明変数を1つ以上指定してください。", call. = FALSE)
    }
    if (any(grepl("`", c(dv, predictors), fixed = TRUE))) stop("列名にバッククォートは使用できません。", call. = FALSE)
    quote_name <- function(x) paste0("`", x, "`")
    stats::as.formula(paste(quote_name(dv), "~", paste(quote_name(predictors), collapse = if (interactions) " * " else " + ")))
  }
  significance <- function(p) {
    ifelse(is.na(p), "", ifelse(p < .001, "***", ifelse(p < .01, "**",
      ifelse(p < .05, "*", ifelse(p < .1, ".", "")))))
  }
  adjust_table <- function(table, adjust = "none") {
    if (is.null(table)) return(NULL)
    table$p.uncorrected <- table$p.value
    table$p.value <- stats::p.adjust(table$p.value, method = adjust)
    table$significance <- significance(table$p.value)
    table
  }
  finalize <- function(result, adjust, alpha, followup_alpha) {
    result$coefficient_table <- adjust_table(result$coefficient_table)
    result$model_table <- adjust_table(result$model_table)
    result$slope_table <- adjust_table(result$slope_table, adjust)
    sizes <- vapply(result$slope_comparisons, function(x) nrow(x$table), integer(1))
    if (sum(sizes)) {
      original <- unlist(lapply(result$slope_comparisons, function(x) x$table$p.value), use.names = FALSE)
      corrected <- stats::p.adjust(original, method = adjust)
      offset <- 0L
      for (name in names(result$slope_comparisons)) {
        table <- adjust_table(result$slope_comparisons[[name]]$table)
        table$p.value <- corrected[offset + seq_len(nrow(table))]
        table$significance <- significance(table$p.value)
        result$slope_comparisons[[name]]$table <- table
        offset <- offset + nrow(table)
      }
    }
    result$adjust <- adjust
    result$alpha <- alpha
    result$followup_alpha <- followup_alpha
    result
  }
  print_table <- function(table, digits, p_digits, estimate_label = "推定値") {
    if (is.null(table) || !nrow(table)) return(invisible(NULL))
    table <- table[, setdiff(names(table), "p.uncorrected"), drop = FALSE]
    last <- intersect(c("p.value", "significance"), names(table))
    table <- table[, c(setdiff(names(table), last), last), drop = FALSE]
    for (name in names(table)) if (is.numeric(table[[name]])) {
      x <- table[[name]]
      p_column <- name == "p.value"
      whole <- name %in% c("N", "df", "df1", "df2") && all(x[is.finite(x)] == floor(x[is.finite(x)]))
      rendered <- rep("—", length(x))
      valid <- is.finite(x)
      rendered[valid] <- formatC(x[valid], format = "f", digits = if (p_column) p_digits else if (whole) 0 else digits)
      if (!p_column) rendered[is.infinite(x)] <- ifelse(x[is.infinite(x)] > 0, "Inf", "-Inf")
      if (p_column) {
        threshold <- 10^(-p_digits)
        rendered[valid & x < threshold] <- paste0("<", formatC(threshold, format = "f", digits = p_digits))
        rendered <- sub("^(<?)0\\.", "\\1.", rendered)
      }
      table[[name]] <- rendered
    } else {
      x <- as.character(table[[name]])
      x[is.na(x)] <- "—"
      table[[name]] <- x
    }
    labels <- c(term = "項", estimate = estimate_label, p.value = "p", significance = "", contrast = "比較",
      lower.CL = "CI下限", upper.CL = "CI上限", adjusted.R2 = "調整済R²", R2 = "R²", t.ratio = "t", z.ratio = "z",
      variable1 = "変数1", variable2 = "変数2", statistic_name = "統計量", statistic = "値")
    names(table) <- vapply(names(table), function(x) if (x %in% names(labels)) labels[[x]] else x, character(1))
    width <- sum(vapply(seq_along(table), function(i) max(nchar(c(names(table)[i], table[[i]]), type = "width")), numeric(1))) + 3L * ncol(table)
    old <- options(width = min(10000L, max(getOption("width"), width)))
    on.exit(options(old), add = TRUE)
    print(table, row.names = FALSE, right = TRUE, quote = FALSE)
    invisible(NULL)
  }
  summaries <- function(table, labels, alpha, threshold) {
    text <- function(index) if (any(index)) paste(labels[index], collapse = ", ") else "なし"
    p <- table$p.value
    cat("有意：", text(is.finite(p) & p < alpha), "\n",
      "有意傾向：", text(is.finite(p) & p >= alpha & p < threshold), "\n\n", sep = "")
  }
  display <- function(result, digits, p_digits, context = NULL) {
    target <- result$dependent_variable
    if (!length(target) || !nzchar(target)) target <- "回帰分析"
    if (is.null(context)) cat("\n========================================\n", target, "\n========================================\n", sep = "")
    heading <- function(text) {
      title <- paste0("──【", text, "】")
      if (nchar(title, type = "width") > 38L) {
        cat("\n", strrep("─", 40), "\n", text, "\n", sep = "")
      } else cat("\n", title, strrep("─", 40L - nchar(title, type = "width")), "\n", sep = "")
    }
    if (is.null(context)) heading(if (is.null(result$coefficient_table)) "既存モデルの傾き分析" else "線形回帰")
    else cat("\n【", context, "】\n", sep = "")
    threshold <- function(x) sub("^0\\.", ".", format(x, scientific = FALSE, trim = TRUE))
    cat("信頼水準 ", 100 * result$ci, "%（信頼区間は未補正）\n",
      "有意 p < ", threshold(result$alpha), " / 有意傾向 ", threshold(result$alpha),
      " ≤ p < ", threshold(result$followup_alpha), "\n", sep = "")
    if (identical(result$inference, "permutation")) {
      info <- result$permutation
      cat("検定：置換検定（両側） / 方法 ", info$method, " / ", info$np, "回\n", sep = "")
      if (!is.null(info$seed)) cat("乱数シード：", info$seed, "\n", sep = "")
    }
    if (!is.null(result$coefficient_table) && is.null(context)) {
      cat("モデル：", paste(deparse(stats::formula(result$model)), collapse = " "), "\n", sep = "")
      cat("\n【モデルの概要】\n")
      print_table(result$model_table, digits, p_digits)
      cat("\n【回帰係数】\n")
      print_table(result$coefficient_table, digits, p_digits)
    }
    if (!is.null(result$slope_table)) {
      heading("条件別の傾き")
      cat("対象変数：", result$trend, "\n",
        "p値の補正：", result$adjust, "（傾きの検定と傾き差の比較を別々に補正）\n",
        sep = "")
      cat("\n【傾きの推定・検定】\n")
      print_table(result$slope_table, digits, p_digits, estimate_label = "傾き b")
      for (name in names(result$slope_comparisons)) {
        comparison <- result$slope_comparisons[[name]]
        table <- comparison$table
        label <- as.character(table$contrast)
        fixed <- comparison$fixed_factors
        if (length(fixed)) {
          conditions <- apply(table[, fixed, drop = FALSE], 1L, function(row) paste(paste0(fixed, "=", row), collapse = ", "))
          label <- paste(conditions, label, sep = "：")
        }
        cat("\n【傾きの比較：", name, "】\n", sep = "")
        summaries(table, label, result$alpha, result$followup_alpha)
        shown <- table[is.finite(table$p.value) & table$p.value < result$followup_alpha, , drop = FALSE]
        if (nrow(shown)) print_table(shown, digits, p_digits, estimate_label = "傾き差 Δb") else cat("表示基準を満たす比較はありません。\n")
      }
    }
    if (identical(result$inference, "permutation")) {
      cat("\n注記：p値は置換検定、SE・信頼区間は線形モデルに基づきます。\n")
      if (result$permutation$repeated_measures) cat("注記：参加者ごとの条件平均・条件差を用い、参加者単位で残差を置換します。\n",
        "注記：親検定の置換法 ", result$permutation$parent_method,
        " に対し、下位分析ではFreedman–Lane法を用います。\n", sep = "")
    }
    translate <- function(x) {
      x <- sub("^NOTE: *", "", x)
      x <- gsub("Results may be misleading due to involvement in interactions",
        "比較対象の要因が交互作用に含まれるため、解釈には注意が必要です。", x, fixed = TRUE)
      x <- gsub("The number of permutations is below 2000, p-values might be unreliable.",
        "置換回数が2000回未満のため、p値の精度に注意が必要です。", x, fixed = TRUE)
      x <- sub("^the distribution of .+ may be discrete\\.$",
        "一部の置換分布が離散的です。p値の精度に注意が必要です。", x)
      x
    }
    messages <- unique(translate(sub("[[:space:]]+$", "", result$messages)))
    warnings <- unique(translate(sub("[[:space:]]+$", "", result$warnings)))
    if (length(messages) || length(warnings)) cat("\n")
    if (length(messages)) cat(paste0("通知：", messages, collapse = "\n"), "\n", sep = "")
    if (length(warnings)) cat(paste0("警告：", warnings, collapse = "\n"), "\n", sep = "")
    invisible(NULL)
  }

  # 親結果の構造・設定の継承はhelperで解決し、計算側に親結果を渡さない。
  permutation_arguments <- function(source) {
    error_call <- function(x) {
      if (!is.call(x)) return(NULL)
      if (identical(x[[1L]], as.name("Error"))) return(x)
      for (part in as.list(x)[-1L]) {out <- error_call(part); if (!is.null(out)) return(out)}
      NULL
    }
    id <- source$id
    error <- if (inherits(source$formula, "formula")) error_call(source$formula[[3L]]) else NULL
    if (is.null(id) && !is.null(error)) id <- all.vars(error[[2L]])[1L]
    if (is.null(id)) {
      needed <- unique(c(source$dependent_variable, source$within, source$between, source$covariate))
      candidates <- setdiff(names(source$analysis_data), needed)
      if (length(candidates) == 1L) id <- candidates
    }
    np <- source$requested_np
    if (is.null(np)) np <- source$np
    if (is.null(np)) np <- source$model$np
    np <- np[1L]
    if (!is.numeric(np) || length(np) != 1L || !is.finite(np) ||
        np < 2 || np > .Machine$integer.max || np != floor(np)) {
      stop("親結果の置換回数を取得できません。", call. = FALSE)
    }
    parent_method <- source$permutation_method
    if (is.null(parent_method)) parent_method <- source$model$method
    method <- if (length(source$within)) "freedman_lane" else parent_method
    if (is.null(method)) method <- "freedman_lane"
    type <- attr(source$model$P, "type")
    if (is.null(type)) type <- "permutation"
    list(data = source$analysis_data, formula = source$formula, model = source$model,
      design = list(dv = source$dependent_variable, id = id, within = source$within,
        between = source$between, covariate = source$covariate),
      permutation = list(np = as.integer(np), seed = source$seed, method = method,
        parent_method = parent_method, type = type, P = source$model$P))
  }

  run_regression <- function(data, formula, model, settings) {
    validate_options(settings$alpha, settings$followup_alpha, settings$adjust, settings$ci,
      settings$digits, settings$p_digits, settings$show, settings$compare_slopes)
    messages <- warnings <- character()
    result <- withCallingHandlers({
      calculator <- r_regression
      arguments <- list(data = data, formula = formula, model = model)
      if (!is.null(settings$permutation)) {
        calculator <- r_permutation_regression
        arguments <- permutation_arguments(settings$permutation)
      }
      do.call(calculator, c(arguments, list(trend = settings$trend, by = settings$by,
        compare_slopes = settings$compare_slopes, ci = settings$ci)))
    },
      message = function(m) {messages <<- c(messages, conditionMessage(m)); invokeRestart("muffleMessage")},
      warning = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")})
    result <- finalize(result, settings$adjust, settings$alpha, settings$followup_alpha)
    result$messages <- unique(messages)
    result$warnings <- unique(warnings)
    class(result) <- unique(c("analysis_relation_result", class(result)))
    result
  }
  analysis_label <- function(method) {
    switch(method, ancova = "ANCOVA", permuco_ancova = "permuco ANCOVA",
      anova = "ANOVA", art = "ART", permuco = "permuco", method)
  }

  result_leaves <- function(input) {
    leaves <- list()
    walk <- function(x, path = character(), method = NULL, depth = 0L) {
      if (!is.list(x) || is.data.frame(x) || depth > 4L ||
          inherits(x, c("lm", "afex_aov", "lmperm", "aovperm", "analysis_relation_result"))) return(invisible(NULL))
      if (length(x$dependent_variable) == 1L &&
          (!is.null(x$test_table) || !is.null(x$anova_table))) {
        kind <- x[["analysis_method"]]
        if (is.null(kind)) kind <- x[["analysis"]]
        if (is.null(kind)) kind <- method
        if (is.null(kind)) kind <- if (inherits(x, "analysis_rank_ancova_result")) "rank_ancova" else
          if (inherits(x, "analysis_ancova_result")) "ancova" else "不明"
        leaves[[length(leaves) + 1L]] <<- list(result = x, method = as.character(kind)[1L],
          target = as.character(x$dependent_variable), path = path)
        return(invisible(NULL))
      }
      module <- attr(x, "method", exact = TRUE)
      if (length(module) == 1L) method <- module
      for (i in seq_along(x)) {
        name <- names(x)[i]
        if (is.null(name) || is.na(name) || !nzchar(name)) name <- as.character(i)
        child_method <- if (name %in% c("ancova", "permuco_ancova", "rank_ancova", "anova", "art", "permuco")) name else method
        walk(x[[i]], c(path, name), child_method, depth + 1L)
      }
      invisible(NULL)
    }
    walk(input)
    leaves
  }
  is_result_input <- function(input) {
    !is.data.frame(input) && (inherits(input, c("analysis_results", "analysis_multi_results",
      "analysis_test", "analysis_test_result")) || length(result_leaves(input)) > 0L)
  }
  inherited_settings <- function(settings, result, inherit) {
    for (name in names(inherit)[inherit]) {
      if (!is.null(result[[name]])) settings[name] <- list(result[[name]])
    }
    settings
  }
  prepare_followups <- function(input, settings, inherit) {
    leaves <- result_leaves(input)
    if (any(vapply(leaves, function(leaf) identical(leaf$method, "rank_ancova") ||
        inherits(leaf$result, "analysis_rank_ancova_result"), logical(1)))) {
      stop("順位ANCOVA（Quade）の結果はrelate()の結果入力に対応していません。共変量との交互作用を検定しないため、条件別回帰は実行できません。",
        call. = FALSE)
    }
    source_analysis <- unique(vapply(leaves, `[[`, character(1), "method"))
    source_analysis <- intersect(source_analysis, c("ancova", "permuco_ancova"))
    if (length(source_analysis) > 1L) {
      stop('ANCOVAと置換ANCOVA（permuco）の両方が含まれています。relate(result$ancova) または relate(result$permuco_ancova) のように、一方の結果を指定してください。',
        call. = FALSE)
    }
    checks <- list(); jobs <- list()
    add_check <- function(method, target, effect = "—", trend = "—", p = NA_real_,
                          alpha = settings$alpha, status = "対象外", reason = "") {
      row <- length(checks) + 1L
      checks[[row]] <<- data.frame(id = row, method = method, dependent_variable = target,
        effect = effect, trend = trend, alpha = alpha, p.value = p, status = status, reason = reason,
        stringsAsFactors = FALSE)
      row
    }
    for (leaf in leaves) {
      source <- leaf$result
      options <- inherited_settings(settings, source, inherit)
      covariates <- source$covariate
      if (!length(covariates) || isFALSE(source$covariate_adjusted) ||
          !leaf$method %in% c("ancova", "permuco_ancova")) {
        add_check(leaf$method, leaf$target, reason = "共変量を調整したANCOVA結果ではありません。")
        next
      }
      table <- source$test_table
      if (is.null(table)) table <- source$anova_table
      if (!is.data.frame(table) || !all(c("Effect", "p.value") %in% names(table)) || !is.numeric(table$p.value)) {
        add_check(leaf$method, leaf$target, status = "確認失敗", reason = "効果名とp値を含む検定表がありません。")
        next
      }
      factors <- unique(c(source$within, source$between, source$followup_factors))
      local_jobs <- list(); found <- FALSE
      for (i in seq_len(nrow(table))) {
        effect <- as.character(table$Effect[i])
        if (is.na(effect)) next
        parts <- tryCatch(all.vars(parse(text = effect)), error = function(e)
          trimws(strsplit(gsub("`", "", effect, fixed = TRUE), ":", fixed = TRUE)[[1L]]))
        cov <- intersect(parts, covariates)
        if (!length(cov)) next
        found <- TRUE
        by <- intersect(parts, factors)
        p <- table$p.value[i]
        status <- "対象"
        reason <- ""
        if (!length(by)) {
          # 共変量の主効果は通常の対象外処理とし、定型注記は生成しない。
          status <- "対象外"
        } else if (!is.finite(p) || p < 0 || p > 1) {
          status <- "判定不可"; reason <- "親検定のp値を取得できません。"
        } else if (p >= options$alpha) status <- "非有意"
        else if (length(cov) != 1L || length(setdiff(parts, c(by, cov)))) {
          status <- "未対応"; reason <- "共変量同士の交互作用など、現在の条件別分析で扱えない効果です。"
        }
        row <- add_check(leaf$method, leaf$target, effect, paste(cov, collapse = ", "), p,
          options$alpha, status, reason)
        if (status != "対象") next
        local_jobs[[length(local_jobs) + 1L]] <- list(row = row, source = source,
          method = leaf$method, target = leaf$target, effect = effect, trend = cov, by = by, settings = options)
      }
      if (!found) add_check(leaf$method, leaf$target, reason = "検定表に共変量を含む効果がありません。")
      # 同じ共変量で高次の交互作用が有意なら、その条件別分析を優先する。
      for (job in local_jobs) {
        covering <- vapply(local_jobs, function(other) identical(other$trend, job$trend) &&
          length(other$by) > length(job$by) && all(job$by %in% other$by), logical(1))
        if (any(covering)) {
          checks[[job$row]]$status <- "高次効果を優先"
          checks[[job$row]]$reason <- paste0("対象の高次効果：", paste(vapply(local_jobs[covering], `[[`, character(1), "effect"), collapse = ", "))
          next
        }
        supported <- if (job$method == "permuco_ancova") "lmperm" else c("lm", "afex_aov")
        if (is.null(job$source$model) || !inherits(job$source$model, supported)) {
          checks[[job$row]]$status <- "未対応"
          checks[[job$row]]$reason <- "傾き分析に対応した推定済みモデルが保存されていません。"
          next
        }
        jobs[[length(jobs) + 1L]] <- job
      }
    }
    empty <- data.frame(id = integer(), method = character(), dependent_variable = character(),
      effect = character(), trend = character(), alpha = numeric(), p.value = numeric(), status = character(), reason = character())
    list(checks = if (length(checks)) do.call(rbind, checks) else empty, jobs = jobs,
      source_analysis = if (length(source_analysis)) source_analysis else NULL)
  }

  choose_followups <- function(plan, initial) {
    if (!length(plan$jobs)) return(list(rows = integer(), settings = initial))
    if (!isTRUE(capabilities("tcltk")) || !suppressWarnings(requireNamespace("tcltk", quietly = TRUE))) {
      stop("GUIにはtcltkが必要です。select = \"arguments\"で自動実行できます。", call. = FALSE)
    }
    window <- tryCatch(tcltk::tktoplevel(), error = function(e) stop("GUIを開けませんでした。\n", conditionMessage(e), call. = FALSE))
    on.exit({try(tcltk::tkgrab.release(window), silent = TRUE); try(tcltk::tkdestroy(window), silent = TRUE)}, add = TRUE)
    tcltk::tkwm.title(window, "ANCOVA結果からの回帰分析")
    tcltk::tkwm.minsize(window, 720, 520)
    tcltk::tkgrid.columnconfigure(window, 1, weight = 1)
    tcltk::tkgrid.rowconfigure(window, 1, weight = 1)
    tcltk::tkgrid(tcltk::tklabel(window, text = "有意な共変量との交互作用から抽出した分析対象です。Ctrl / Shiftで複数選択できます。",
      anchor = "w", wraplength = 680), row = 0, column = 0, columnspan = 2, sticky = "w", padx = 12, pady = 8)
    box <- tcltk::tklistbox(window, selectmode = "extended", exportselection = FALSE, height = 7, width = 75)
    for (job in plan$jobs) tcltk::tkinsert(box, "end", paste(job$target, analysis_label(job$method), job$effect, sep = " / "))
    tcltk::tkselection.set(box, 0, "end")
    tcltk::tkgrid(box, row = 1, column = 0, columnspan = 2, sticky = "nsew", padx = 12, pady = 8)
    vars <- list()
    for (i in seq_along(c("ci", "digits", "p_digits"))) {
      name <- c("ci", "digits", "p_digits")[i]
      vars[[name]] <- tcltk::tclVar(as.character(initial[[name]]))
      label <- c(ci = "信頼水準", digits = "数値の桁数", p_digits = "p値の桁数")[[name]]
      tcltk::tkgrid(tcltk::tklabel(window, text = label), row = i + 1L, column = 0, sticky = "w", padx = 12, pady = 6)
      tcltk::tkgrid(tcltk::tkentry(window, textvariable = vars[[name]]), row = i + 1L, column = 1, sticky = "ew", padx = 12, pady = 6)
    }
    adjust_var <- tcltk::tclVar(initial$adjust)
    tcltk::tkgrid(tcltk::tklabel(window, text = "p値の補正法"), row = 5, column = 0, sticky = "w", padx = 12)
    tcltk::tkgrid(tcltk::ttkcombobox(window, textvariable = adjust_var, values = stats::p.adjust.methods, state = "readonly"),
      row = 5, column = 1, sticky = "ew", padx = 12)
    compare_var <- tcltk::tclVar(if (initial$compare_slopes) "1" else "0")
    show_var <- tcltk::tclVar(if (initial$show) "1" else "0")
    tcltk::tkgrid(tcltk::tkcheckbutton(window, text = "条件間の傾きも比較する", variable = compare_var), row = 6, column = 1, sticky = "w")
    tcltk::tkgrid(tcltk::tkcheckbutton(window, text = "コンソールに結果を表示", variable = show_var), row = 7, column = 1, sticky = "w")
    error_var <- tcltk::tclVar("")
    tcltk::tkgrid(tcltk::tklabel(window, textvariable = error_var, foreground = "#a32020", wraplength = 680),
      row = 8, column = 0, columnspan = 2, sticky = "w", padx = 12)
    chosen <- NULL
    confirm <- function() {
      settings <- tryCatch({
        options <- initial
        for (name in names(vars)) options[[name]] <- suppressWarnings(as.numeric(tcltk::tclvalue(vars[[name]])))
        options$adjust <- tcltk::tclvalue(adjust_var)
        options$compare_slopes <- tcltk::tclvalue(compare_var) == "1"
        options$show <- tcltk::tclvalue(show_var) == "1"
        validate_options(options$alpha, options$followup_alpha, options$adjust, options$ci,
          options$digits, options$p_digits, options$show, options$compare_slopes)
        options
      }, error = function(e) {tcltk::tclvalue(error_var) <- conditionMessage(e); NULL})
      if (is.null(settings)) return(invisible(NULL))
      index <- as.integer(tcltk::tkcurselection(box)) + 1L
      chosen <<- list(rows = vapply(plan$jobs[index], `[[`, integer(1), "row"), settings = settings)
      tcltk::tkdestroy(window)
    }
    buttons <- tcltk::tkframe(window)
    tcltk::tkgrid(buttons, row = 9, column = 1, sticky = "e", padx = 12, pady = 10)
    cancel <- function() tcltk::tkdestroy(window)
    tcltk::tkpack(tcltk::tkbutton(buttons, text = "キャンセル", command = cancel), side = "left", padx = 4)
    tcltk::tkpack(tcltk::tkbutton(buttons, text = "実行", command = confirm), side = "left", padx = 4)
    tcltk::tkwm.protocol(window, "WM_DELETE_WINDOW", cancel)
    tcltk::tkbind(window, "<Escape>", cancel)
    if (!as.logical(as.integer(tcltk::tkwinfo("viewable", window)))) tcltk::tkwait.visibility(window)
    tcltk::tkgrab.set(window)
    tcltk::tkwait.window(window)
    if (is.null(chosen)) stop("設定をキャンセルしました。回帰分析は実行していません。", call. = FALSE)
    chosen
  }

  followup_results <- function(input, settings, inherit, select, result_rows = NULL) {
    plan <- prepare_followups(input, settings, inherit)
    if (select == "gui" && length(plan$jobs)) {
      chosen <- choose_followups(plan, settings)
      for (i in seq_along(plan$jobs)) {
        job <- plan$jobs[[i]]
        if (!job$row %in% chosen$rows) plan$checks$status[job$row] <- "選択対象外"
        # 対象判定の閾値は親結果から保持し、GUIでは表示・補正を変更する。
        for (name in c("ci", "digits", "p_digits", "adjust", "compare_slopes", "show")) {
          plan$jobs[[i]]$settings[name] <- chosen$settings[name]
        }
      }
      for (name in c("ci", "digits", "p_digits", "adjust", "compare_slopes", "show")) settings[name] <- chosen$settings[name]
    }
    if (!is.null(result_rows)) {
      if (!is.numeric(result_rows) || anyNA(result_rows) || any(!is.finite(result_rows)) ||
          any(result_rows != floor(result_rows)) || any(!result_rows %in% plan$checks$id)) {
        stop("result_rowsには確認表に含まれる行番号を指定してください。", call. = FALSE)
      }
      omitted <- plan$checks$status == "対象" & !plan$checks$id %in% result_rows
      plan$checks$status[omitted] <- "選択対象外"
    }
    results <- list()
    for (job in plan$jobs) {
      if (plan$checks$status[job$row] != "対象") next
      options <- job$settings
      options$trend <- job$trend; options$by <- job$by
      permutation <- identical(job$method, "permuco_ancova")
      options$permutation <- if (permutation) job$source else NULL
      model <- if (permutation) NULL else job$source$model
      computed <- tryCatch(run_regression(NULL, NULL, model, options), error = function(e) {
        plan$checks$status[job$row] <<- "実行失敗"
        plan$checks$reason[job$row] <<- conditionMessage(e)
        NULL
      })
      if (is.null(computed)) next
      computed$dependent_variable <- job$target
      computed$source_analysis <- job$method
      computed$parent_effect <- job$effect
      computed$parent_p.value <- plan$checks$p.value[job$row]
      computed$parent_alpha <- options$alpha
      computed$source_row <- job$row
      computed$display_options <- options[c("digits", "p_digits")]
      key <- paste(job$method, job$target, job$effect, sep = "/")
      key <- tail(make.unique(c(names(results), key)), 1L)
      results[[key]] <- computed
      plan$checks$status[job$row] <- "実行済み"
    }
    output <- list(analysis = "ancova_followup", checks = plan$checks, results = results,
      source_analysis = plan$source_analysis,
      input_type = "analysis_results", settings = settings,
      selected_rows = plan$checks$id[plan$checks$status %in% c("実行済み", "実行失敗")], show = settings$show)
    class(output) <- c("analysis_relation_followups", "list")
    output
  }
  display_followups <- function(output, digits, p_digits) {
    # 共変量の確認ではANCOVA系だけを表示し、他の方法の対象外行は内部記録に残す。
    checks <- output$checks[output$checks$method %in% c("ancova", "permuco_ancova"), , drop = FALSE]
    if (!nrow(checks)) {
      cat("\nANCOVAの分析結果がありません。回帰分析は実行していません。\n")
      return(invisible(NULL))
    }
    for (target in unique(checks$dependent_variable)) {
      cat("\n========================================\n", target, "\n========================================\n", sep = "")
      text <- "──【共変量の確認】"
      cat("\n", text, strrep("─", 40L - nchar(text, type = "width")), "\n", sep = "")
      rows <- checks[checks$dependent_variable == target, , drop = FALSE]
      groups <- unique(rows[, c("method", "alpha"), drop = FALSE])
      for (i in seq_len(nrow(groups))) {
        group <- groups[i, , drop = FALSE]
        selected <- rows$method == group$method & rows$alpha == group$alpha
        threshold <- sub("^0\\.", ".", format(group$alpha, scientific = FALSE, trim = TRUE))
        cat("\n方法：", analysis_label(group$method), "\n実行基準：p < ", threshold, "\n\n", sep = "")
        shown <- rows[selected, c("effect", "trend", "p.value"), drop = FALSE]
        shown$significance <- significance(shown$p.value)
        names(shown)[1:2] <- c("効果", "共変量")
        print_table(shown, digits, p_digits)
        notes <- rows[selected & nzchar(rows$reason), , drop = FALSE]
        if (nrow(notes)) cat("\n")
        for (j in seq_len(nrow(notes))) cat("注記 [", analysis_label(notes$method[j]), " / ", notes$effect[j], "]：", notes$reason[j], "\n", sep = "")
      }
      for (result in output$results) if (identical(result$dependent_variable, target)) {
        display(result, result$display_options$digits, result$display_options$p_digits,
          context = paste0(result$parent_effect, "（", analysis_label(result$source_analysis), "）"))
      }
      if (!any(rows$status == "実行済み")) cat("\n実行した条件別回帰分析はありません。\n")
    }
    invisible(NULL)
  }

  correlation_label <- function(method) switch(method, pearson = "Pearson", spearman = "Spearman", kendall = "Kendall")
  correlation_symbol <- function(method) switch(method, pearson = "r", spearman = "ρ", kendall = "τ")
  correlation_choices <- function() list(
    method = c("Pearson（積率相関）" = "pearson", "Spearman（順位相関）" = "spearman", "Kendall（順位相関）" = "kendall"),
    use = c("ペアごとに欠測を除外" = "pairwise.complete.obs", "選択した全変数で欠測を除外" = "complete.obs"),
    alternative = c("両側" = "two.sided", "片側：正の相関" = "greater", "片側：負の相関" = "less"),
    exact = list("自動" = NULL, "小標本向けの計算法" = TRUE, "漸近近似" = FALSE))

  validate_correlation_variables <- function(data, variables) {
    if (!is.data.frame(data) || anyDuplicated(names(data))) stop("列名が重複しないデータフレームを指定してください。", call. = FALSE)
    if (!is.character(variables) || length(variables) < 2L || anyNA(variables) ||
        anyDuplicated(variables) || any(!variables %in% names(data)) ||
        any(!vapply(data[variables], function(x) is.numeric(x) && is.null(dim(x)), logical(1)))) {
      stop("異なる数値変数を2つ以上選択してください。", call. = FALSE)
    }
    invisible(NULL)
  }

  display_correlation <- function(result, digits, p_digits) {
    cat("\n========================================\n相関分析\n========================================\n")
    title <- "──【相関係数の検定】"
    cat("\n", title, strrep("─", 40L - nchar(title, type = "width")), "\n", sep = "")
    choices <- correlation_choices()
    threshold <- function(x) sub("^0\\.", ".", format(x, scientific = FALSE, trim = TRUE))
    cat("方法：", correlation_label(result$method), " / 検定：", names(choices$alternative)[match(result$alternative, choices$alternative)],
      "\n欠測の処理：", names(choices$use)[match(result$use, choices$use)],
      "\np値の補正：", result$adjust, "（選択した全ペアをまとめて補正）\n",
      "有意 p < ", threshold(result$alpha), " / 有意傾向 ", threshold(result$alpha), " ≤ p < ", threshold(result$followup_alpha), "\n", sep = "")
    if (result$method == "pearson") cat("信頼水準 ", 100 * result$ci, "%（信頼区間は未補正）\n", sep = "")
    else {
      index <- which(vapply(choices$exact, identical, logical(1), result$exact))
      cat("順位相関のp値計算：", names(choices$exact)[index], "\n", sep = "")
    }
    table <- result$correlation_table
    labels <- paste(table$variable1, table$variable2, sep = " - ")
    cat("\n")
    summaries(table, labels, result$alpha, result$followup_alpha)
    shown <- table
    for (name in c("df", "lower.CL", "upper.CL")) if (all(is.na(shown[[name]]))) shown[[name]] <- NULL
    print_table(shown, digits, p_digits, estimate_label = correlation_symbol(result$method))
    notes <- if (result$method != "pearson") "順位相関の信頼区間は、この版では算出しません。" else if (any(table$N < 4L))
      "Pearson相関の信頼区間は、有効件数Nが4以上のペアで算出します。" else character()
    if (result$use == "pairwise.complete.obs" && length(unique(table$N)) > 1L) notes <- c(notes,
      "ペアごとに欠測を除外しているため、有効件数Nが異なります。")
    translate <- function(text) {
      text <- gsub("Cannot compute exact p-value with ties", "同順位があるため、小標本向けのp値を計算できず、近似p値を使用しました。", text, fixed = TRUE)
      text
    }
    warnings <- unlist(lapply(result$diagnostics, function(x) {
      label <- paste(x$variable1, x$variable2, sep = " - ")
      c(if (!is.null(x$error)) paste0(label, "：", x$error),
        if (length(x$warnings)) paste0(label, "：", translate(x$warnings)))
    }), use.names = FALSE)
    if (length(notes) || length(warnings)) cat("\n")
    if (length(notes)) cat(paste0("注記：", notes, collapse = "\n"), "\n", sep = "")
    if (length(warnings)) cat(paste0("警告：", unique(warnings), collapse = "\n"), "\n", sep = "")
    invisible(NULL)
  }

  choose_correlation <- function(data, variables, initial) {
    if (!is.data.frame(data) || anyDuplicated(names(data))) stop("列名が重複しないデータフレームを指定してください。", call. = FALSE)
    numeric <- names(data)[vapply(data, function(x) is.numeric(x) && is.null(dim(x)), logical(1))]
    if (length(numeric) < 2L) stop("相関分析には数値列が2つ以上必要です。", call. = FALSE)
    if (!isTRUE(capabilities("tcltk")) || !suppressWarnings(requireNamespace("tcltk", quietly = TRUE))) stop('GUIにはtcltkが必要です。variablesを指定しselect = "arguments"で実行できます。', call. = FALSE)
    window <- tryCatch(tcltk::tktoplevel(), error = function(e) stop("GUIを開けませんでした。\n", conditionMessage(e), call. = FALSE))
    on.exit({try(tcltk::tkgrab.release(window), silent = TRUE); try(tcltk::tkdestroy(window), silent = TRUE)}, add = TRUE)
    tcltk::tkwm.title(window, "関係性分析：相関分析の設定")
    tcltk::tkwm.minsize(window, 720, 620)
    tcltk::tkgrid.columnconfigure(window, 1, weight = 1)
    tcltk::tkgrid.rowconfigure(window, 1, weight = 1)
    tcltk::tkgrid(tcltk::tklabel(window, text = "相関を調べる数値変数を2つ以上選択します（Ctrl / Shiftで複数選択）。\n1行を1観測として扱い、選択した全ペアを検定します。",
      anchor = "w", justify = "left", wraplength = 680), row = 0, column = 0, columnspan = 2, sticky = "w", padx = 12, pady = 8)
    box <- tcltk::tklistbox(window, selectmode = "extended", exportselection = FALSE, height = 6, width = 45)
    for (name in numeric) tcltk::tkinsert(box, "end", name)
    for (i in match(variables, numeric, nomatch = 0L)) if (i > 0L) tcltk::tkselection.set(box, i - 1L)
    tcltk::tkgrid(box, row = 1, column = 0, columnspan = 2, sticky = "nsew", padx = 12, pady = 4)
    maps <- correlation_choices()
    maps$adjust <- setNames(as.list(stats::p.adjust.methods), stats::p.adjust.methods)
    fields <- vars <- list()
    labels <- c(method = "相関係数の種類", use = "欠測の処理", alternative = "検定の方向", exact = "順位相関のp値計算", adjust = "p値の補正法",
      alpha = "有意水準", followup_alpha = "有意傾向の基準", ci = "信頼水準（Pearson）", digits = "数値の桁数", p_digits = "p値の桁数")
    for (i in seq_along(labels)) {
      name <- names(labels)[i]
      value <- initial[[name]]
      if (name %in% names(maps)) {
        values <- as.list(maps[[name]])
        index <- which(vapply(values, function(x) isTRUE(all.equal(x, value)), logical(1)))[1L]
        value <- names(maps[[name]])[index]
      } else value <- as.character(value)
      vars[[name]] <- tcltk::tclVar(value)
      tcltk::tkgrid(tcltk::tklabel(window, text = labels[[name]], anchor = "w"), row = i + 1L, column = 0, sticky = "w", padx = 12, pady = 4)
      fields[[name]] <- if (name %in% names(maps)) tcltk::ttkcombobox(window, textvariable = vars[[name]], values = names(maps[[name]]), state = "readonly") else
        tcltk::tkentry(window, textvariable = vars[[name]])
      tcltk::tkgrid(fields[[name]], row = i + 1L, column = 1, sticky = "ew", padx = 12, pady = 4)
    }
    refresh <- function() {
      pearson <- tcltk::tclvalue(vars$method) == names(maps$method)[1L]
      tcltk::tkconfigure(fields$ci, state = if (pearson) "normal" else "disabled")
      tcltk::tkconfigure(fields$exact, state = if (pearson) "disabled" else "readonly")
    }
    tcltk::tkbind(fields$method, "<<ComboboxSelected>>", refresh)
    refresh()
    show_var <- tcltk::tclVar(if (initial$show) "1" else "0")
    tcltk::tkgrid(tcltk::tkcheckbutton(window, text = "コンソールに結果を表示", variable = show_var), row = 12, column = 1, sticky = "w", padx = 12)
    error_var <- tcltk::tclVar("")
    tcltk::tkgrid(tcltk::tklabel(window, textvariable = error_var, foreground = "#a32020", wraplength = 680), row = 13, column = 0, columnspan = 2, sticky = "w", padx = 12)
    result <- NULL
    confirm <- function() {
      result_candidate <- tryCatch({
        selected <- numeric[as.integer(tcltk::tkcurselection(box)) + 1L]
        validate_correlation_variables(data, selected)
        settings <- initial
        for (name in names(vars)) {
          value <- tcltk::tclvalue(vars[[name]])
          settings[name] <- list(if (name %in% names(maps)) as.list(maps[[name]])[[value]] else suppressWarnings(as.numeric(value)))
        }
        settings$show <- tcltk::tclvalue(show_var) == "1"
        validate_options(settings$alpha, settings$followup_alpha, settings$adjust, settings$ci, settings$digits, settings$p_digits, settings$show, TRUE)
        list(variables = selected, settings = settings)
      }, error = function(e) {tcltk::tclvalue(error_var) <- conditionMessage(e); NULL})
      if (is.null(result_candidate)) return(invisible(NULL))
      result <<- result_candidate
      tcltk::tkdestroy(window)
    }
    buttons <- tcltk::tkframe(window)
    tcltk::tkgrid(buttons, row = 14, column = 1, sticky = "e", padx = 12, pady = 8)
    cancel <- function() tcltk::tkdestroy(window)
    tcltk::tkpack(tcltk::tkbutton(buttons, text = "キャンセル", command = cancel), side = "left", padx = 4)
    tcltk::tkpack(tcltk::tkbutton(buttons, text = "実行", command = confirm), side = "left", padx = 4)
    tcltk::tkwm.protocol(window, "WM_DELETE_WINDOW", cancel)
    tcltk::tkbind(window, "<Escape>", cancel)
    if (!as.logical(as.integer(tcltk::tkwinfo("viewable", window)))) tcltk::tkwait.visibility(window)
    tcltk::tkgrab.set(window)
    tcltk::tkwait.window(window)
    if (is.null(result)) stop("設定をキャンセルしました。相関分析は実行していません。", call. = FALSE)
    result
  }

  correlate <- function(data, variables, settings, select, data_expression) {
    if (select == "gui") {
      selected <- choose_correlation(data, variables, settings)
      variables <- selected$variables; settings <- selected$settings
    }
    validate_correlation_variables(data, variables)
    result <- r_correlation(data, variables, settings$method, settings$use, settings$alternative, settings$ci, settings$exact)
    result$p.uncorrected_matrix <- result$p.value_matrix
    result$correlation_table <- adjust_table(result$correlation_table, settings$adjust)
    for (i in seq_len(nrow(result$correlation_table))) {
      row <- result$correlation_table[i, ]
      result$p.value_matrix[row$variable1, row$variable2] <- result$p.value_matrix[row$variable2, row$variable1] <- row$p.value
    }
    result$adjust <- settings$adjust; result$alpha <- settings$alpha; result$followup_alpha <- settings$followup_alpha
    args <- c(list(data = data_expression, analysis = "correlation", variables = variables), settings, list(select = "arguments"))
    call <- as.call(c(list(as.name("relate")), args))
    code <- deparse(call, width.cutoff = 80L); code[1L] <- paste0("relation_result <- ", code[1L])
    result$selection <- list(method = select, call = call, code = paste(code, collapse = "\n"))
    class(result) <- c("analysis_relation_result", class(result))
    if (settings$show) {
      display_correlation(result, settings$digits, settings$p_digits)
      if (select == "gui") cat("\n【再実行用コード】\n", result$selection$code, "\n", sep = "")
    }
    invisible(result)
  }

  sem_choices <- function() list(
    estimator = c("自動（連続値：ML／順序尺度：WLSMV）" = "auto", "ML" = "ML", "MLR（頑健推定）" = "MLR", "WLSMV（順序尺度）" = "WLSMV"),
    missing_data = c("欠測のある行を除外" = "listwise", "完全情報最尤法（FIML）" = "fiml"),
    se = c("推定法に対応する標準誤差" = "default", "ブートストラップ" = "bootstrap"))

  validate_sem_options <- function(settings) {
    validate_options(settings$alpha, settings$followup_alpha, settings$adjust, settings$ci, settings$digits, settings$p_digits, settings$show, TRUE)
    for (name in c("std_lv", "fixed_x", "meanstructure")) if (!is.logical(settings[[name]]) || length(settings[[name]]) != 1L || is.na(settings[[name]])) stop(name, "にはTRUEまたはFALSEを指定してください。", call. = FALSE)
    for (name in c("estimator", "missing_data", "se")) {
      if (!is.character(settings[[name]]) || length(settings[[name]]) != 1L || is.na(settings[[name]]) || !settings[[name]] %in% unname(sem_choices()[[name]])) stop(name, "の選択値を確認してください。", call. = FALSE)
    }
    scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
    if (!scalar(settings$bootstrap) || settings$bootstrap < 2 || settings$bootstrap > .Machine$integer.max || settings$bootstrap != floor(settings$bootstrap)) stop("bootstrapは2以上の整数にしてください。", call. = FALSE)
    if (!is.null(settings$seed) && (!scalar(settings$seed) || settings$seed < 0 || settings$seed > .Machine$integer.max || settings$seed != floor(settings$seed))) stop("seedは空欄、または0以上の整数にしてください。", call. = FALSE)
    invisible(NULL)
  }

  prepare_sem <- function(data, model, settings) {
    if (!requireNamespace("lavaan", quietly = TRUE)) stop('SEMにはlavaanが必要です。install.packages("lavaan")で導入してください。', call. = FALSE)
    if (!is.data.frame(data) || !nrow(data) || anyDuplicated(names(data))) stop("列名が重複しない個票データフレームを指定してください。", call. = FALSE)
    if (!is.character(model) || !length(model) || anyNA(model) || !nzchar(trimws(paste(model, collapse = "\n")))) stop("lavaan形式のモデル式を入力してください。", call. = FALSE)
    validate_sem_options(settings)
    parsed <- tryCatch(lavaan::lavaanify(paste(model, collapse = "\n"), model.type = "sem", auto = TRUE),
      error = function(e) stop("モデル式を解釈できません。\n", conditionMessage(e), call. = FALSE))
    observed <- lavaan::lavNames(parsed, "ov")
    if (!length(observed)) stop("モデル式に観測変数を指定してください。", call. = FALSE)
    absent <- setdiff(observed, names(data))
    if (length(absent)) stop("モデル式の変数がデータにありません：", paste(absent, collapse = ", "), call. = FALSE)
    ordered <- settings$ordered
    if (!is.null(ordered) && (!is.character(ordered) || anyNA(ordered) || anyDuplicated(ordered) || any(!ordered %in% observed))) stop("orderedには、モデルで使用する順序尺度の列名を指定してください。", call. = FALSE)
    ordered <- unique(c(ordered, observed[vapply(data[observed], is.ordered, logical(1))]))
    for (name in observed) {
      x <- data[[name]]
      if (!is.null(dim(x)) || !(is.numeric(x) || (name %in% ordered && is.factor(x)))) stop("観測変数は数値、または指定した順序因子にしてください：", name, call. = FALSE)
      if (is.numeric(x) && any(is.infinite(x))) stop("観測変数に無限大が含まれています：", name, call. = FALSE)
    }
    estimator <- settings$estimator
    if (estimator == "auto") estimator <- if (length(ordered)) "WLSMV" else "ML"
    if (length(ordered) && estimator != "WLSMV") stop("順序尺度にはestimator = \"WLSMV\"または\"auto\"を指定してください。", call. = FALSE)
    if (settings$missing_data == "fiml" && estimator == "WLSMV") stop("WLSMVではFIMLを利用できません。欠測のある行を除外する設定にしてください。", call. = FALSE)
    if (settings$missing_data == "fiml" && !settings$meanstructure) stop("FIMLでは平均構造を含めてください。", call. = FALSE)
    if (settings$se == "bootstrap" && estimator != "ML") stop("この版のブートストラップはML推定で利用してください。", call. = FALSE)
    list(data = data, model = paste(model, collapse = "\n"), estimator = estimator,
      missing_data = settings$missing_data, ordered = if (length(ordered)) ordered else NULL,
      std_lv = settings$std_lv, fixed_x = settings$fixed_x, meanstructure = settings$meanstructure,
      se = settings$se, bootstrap = as.integer(settings$bootstrap), seed = settings$seed, ci = settings$ci)
  }

  translate_sem_notice <- function(text) {
    translations <- c(
      "some estimated ov variances are negative" = "観測変数の推定分散に負の値が含まれています。",
      "some estimated lv variances are negative" = "潜在変数の推定分散に負の値が含まれています。",
      "covariance matrix of latent variables is not positive definite" = "潜在変数の共分散行列が正定値ではありません。",
      "Could not compute standard errors" = "標準誤差を計算できませんでした。モデルの識別と情報行列を確認してください。",
      "Model estimation FAILED" = "モデルの推定が収束しませんでした。",
      "not all bootstrap runs converged" = "一部のブートストラップ標本で推定が収束しませんでした。",
      "NaNs produced" = "計算過程で非数値（NaN）が発生しました。",
      "extreme order statistics used as endpoints" = "信頼区間の端点が極端な順位統計量に基づいています。ブートストラップ回数を確認してください。",
      "first-order delta method" = "定義した効果の標準誤差・信頼区間は一次のデルタ法に基づきます。",
      "due to missing values in exogenous" = "固定として扱う外生共変量に欠測があるため、一部の行が分析から除外されました。")
    vapply(text, function(x) {
      found <- names(translations)[vapply(names(translations), function(pattern) grepl(pattern, x, fixed = TRUE), logical(1))]
      if (length(found)) paste(unname(translations[found]), collapse = " ") else trimws(gsub("[[:space:]]+", " ", x))
    }, character(1), USE.NAMES = FALSE)
  }

  display_sem <- function(result, digits, p_digits) {
    cat("\n========================================\nSEM\n========================================\n")
    heading <- function(text) {title <- paste0("──【", text, "】"); cat("\n", title, strrep("─", max(0L, 40L - nchar(title, type = "width"))), "\n", sep = "")}
    heading("モデルと推定設定")
    se_label <- switch(result$options$se, standard = "通常", robust.huber.white = "頑健（Huber–White）",
      robust.sem = "頑健（推定法に対応）", bootstrap = "ブートストラップ", result$options$se)
    cat("推定法：", result$estimator, " / 欠測：", if (result$missing_data == "fiml") "FIML" else "行単位で除外",
      "\n標準誤差：", se_label, "\nデータ件数：", result$N_input, " / 分析件数N：", result$N,
      "\n信頼水準：", 100 * result$ci, "%（信頼区間は未補正）\np値の補正：", result$adjust,
      "（検定可能な全パラメータをまとめて補正）\n", sep = "")
    threshold <- function(x) sub("^0\\.", ".", format(x, scientific = FALSE, trim = TRUE))
    cat("有意 p < ", threshold(result$alpha), " / 有意傾向 ", threshold(result$alpha), " ≤ p < ", threshold(result$followup_alpha), "\n", sep = "")
    if (length(result$ordered)) cat("順序尺度：", paste(result$ordered, collapse = ", "), "\n", sep = "")
    if (result$se == "bootstrap") cat("ブートストラップ：", result$bootstrap, "回 / CI：パーセンタイル法\n", sep = "")
    if (!is.null(result$seed)) cat("乱数シード：", result$seed, "\n", sep = "")
    cat("\n【モデル式】\n", result$model_syntax, "\n", sep = "")
    if (result$converged) {
      heading("モデル適合度")
      if (!is.null(result$fit_table)) {
        shown <- result$fit_table
        for (name in setdiff(names(shown), "variant")) if (all(is.na(shown[[name]]))) shown[[name]] <- NULL
        names(shown)[names(shown) == "variant"] <- "種類"
        names(shown)[names(shown) == "chisq"] <- "χ²"
        names(shown)[names(shown) == "RMSEA.lower"] <- "RMSEA下限"
        names(shown)[names(shown) == "RMSEA.upper"] <- "RMSEA上限"
        print_table(shown, digits, p_digits)
        if ("df" %in% names(result$fit_measures) && isTRUE(unname(result$fit_measures["df"]) == 0)) cat("\n注記：自由度が0のため、モデル全体の適合度による評価はできません。\n")
      } else cat("適合度を取得できませんでした。\n")
      table <- result$parameter_table
      if (!is.null(table)) {
        groups <- list("回帰パス" = table$op == "~", "因子負荷量" = table$op == "=~",
          "定義した効果" = table$op == ":=", "分散・共分散" = table$op == "~~",
          "切片・平均" = table$op == "~1", "閾値" = table$op == "|")
        for (title in names(groups)) {
          own <- table[groups[[title]], , drop = FALSE]
          if (!nrow(own)) next
          heading(title)
          label <- if (title %in% c("回帰パス", "因子負荷量")) paste(own$rhs, own$lhs, sep = if (title == "回帰パス") " → " else " ← ") else
            if (title == "分散・共分散") ifelse(own$lhs == own$rhs, paste0(own$lhs, "（分散）"), paste(own$lhs, own$rhs, sep = " ↔ ")) else
              if (title == "閾値") paste(own$lhs, own$rhs, sep = " | ") else own$lhs
          if (title %in% c("回帰パス", "定義した効果")) summaries(own, label, result$alpha, result$followup_alpha)
          shown <- data.frame(関係 = label, estimate = own$est, std.all = own$std.all, SE = own$se, z.ratio = own$z,
            lower.CL = own$ci.lower, upper.CL = own$ci.upper, p.value = own$p.value, significance = own$significance, check.names = FALSE)
          names(shown)[names(shown) == "std.all"] <- if (title == "回帰パス") "標準化 β" else if (title == "因子負荷量") "標準化 λ" else "標準化推定値"
          print_table(shown, digits, p_digits, estimate_label = if (title == "回帰パス") "係数 b" else if (title == "因子負荷量") "負荷量 λ" else "推定値")
        }
        cat("\n注記：SE・z値・信頼区間・p値は非標準化推定値に対する値です。固定した係数の検定欄は「—」で表示します。\n")
      }
      if (!is.null(result$r2_table)) {heading("決定係数"); shown <- result$r2_table; names(shown)[names(shown) == "variable"] <- "変数"; print_table(shown, digits, p_digits)}
    }
    warnings <- unique(c(if (!result$converged) "モデルが収束していないため、係数の検定・適合度・R²は表示しません。",
      if (isFALSE(result$admissible)) "不適解の可能性があります。負の分散や共分散行列などを確認してください。",
      translate_sem_notice(result$warnings), unlist(result$diagnostics$errors, use.names = FALSE)))
    messages <- unique(translate_sem_notice(result$messages))
    if (length(warnings) || length(messages)) cat("\n")
    if (length(messages)) cat(paste0("通知：", messages, collapse = "\n"), "\n", sep = "")
    if (length(warnings)) cat(paste0("警告：", warnings, collapse = "\n"), "\n", sep = "")
    invisible(NULL)
  }

  run_sem <- function(data, model, settings, select, data_expression) {
    if (select == "gui") {selected <- choose_sem(data, model, settings); model <- selected$model; settings <- selected$settings}
    messages <- warnings <- character()
    result <- withCallingHandlers({arguments <- prepare_sem(data, model, settings); do.call(r_sem, arguments)},
      message = function(m) {messages <<- c(messages, conditionMessage(m)); invokeRestart("muffleMessage")},
      warning = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")})
    result$parameter_table <- result$parameter_estimates
    if (!is.null(result$parameter_table)) {
      table <- result$parameter_table
      names(table)[names(table) == "pvalue"] <- "p.value"
      table$p.uncorrected <- table$p.value
      testable <- is.finite(table$p.value) & table$se > 0
      table$p.value[testable] <- stats::p.adjust(table$p.value[testable], method = settings$adjust)
      table$significance <- significance(table$p.value)
      result$parameter_table <- table
    }
    if (!is.null(result$fit_table)) result$fit_table$significance <- significance(result$fit_table$p.value)
    result$adjust <- settings$adjust; result$alpha <- settings$alpha; result$followup_alpha <- settings$followup_alpha
    result$warnings <- unique(warnings); result$messages <- unique(messages)
    args <- c(list(data = data_expression, analysis = "sem", sem_model = result$model_syntax), settings, list(select = "arguments"))
    call <- as.call(c(list(as.name("relate")), args))
    code <- deparse(call, width.cutoff = 80L); code[1L] <- paste0("relation_result <- ", code[1L])
    result$selection <- list(method = select, call = call, code = paste(code, collapse = "\n"))
    class(result) <- c("analysis_relation_result", class(result))
    if (settings$show) {display_sem(result, settings$digits, settings$p_digits); if (select == "gui") cat("\n【再実行用コード】\n", result$selection$code, "\n", sep = "")}
    invisible(result)
  }

  choose_sem <- function(data, model, initial) {
    validate_sem_options(initial)
    if (!requireNamespace("lavaan", quietly = TRUE)) stop('SEMにはlavaanが必要です。install.packages("lavaan")で導入してください。', call. = FALSE)
    if (!is.data.frame(data) || !nrow(data) || anyDuplicated(names(data))) stop("列名が重複しない個票データを指定してください。", call. = FALSE)
    if (!isTRUE(capabilities("tcltk")) || !suppressWarnings(requireNamespace("tcltk", quietly = TRUE))) stop('GUIにはtcltkが必要です。sem_modelとselect = "arguments"を指定して実行できます。', call. = FALSE)
    window <- tryCatch(tcltk::tktoplevel(), error = function(e) stop("GUIを開けませんでした。\n", conditionMessage(e), call. = FALSE))
    on.exit({try(tcltk::tkgrab.release(window), silent = TRUE); try(tcltk::tkdestroy(window), silent = TRUE)}, add = TRUE)
    tcltk::tkwm.title(window, "関係性分析：SEMの設定")
    tcltk::tkwm.minsize(window, 850, 650)
    tcltk::tkgrid.columnconfigure(window, 0, weight = 1)
    tcltk::tkgrid.rowconfigure(window, 0, weight = 1)
    notebook <- tcltk::ttknotebook(window)
    tcltk::tkgrid(notebook, row = 0, column = 0, sticky = "nsew", padx = 12, pady = 10)
    panels <- list()
    for (name in c("モデル", "推定", "検定・表示")) {
      panels[[name]] <- tcltk::ttkframe(notebook, padding = 12)
      tcltk::tkadd(notebook, panels[[name]], text = name)
      tcltk::tkgrid.columnconfigure(panels[[name]], 1, weight = 1)
    }
    panel <- panels[["モデル"]]
    tcltk::tkgrid.rowconfigure(panel, 1, weight = 1)
    tcltk::tkgrid(tcltk::tklabel(panel, text = "データの列名"), row = 0, column = 0, sticky = "w")
    tcltk::tkgrid(tcltk::tklabel(panel, text = "モデル式（lavaan形式）"), row = 0, column = 1, sticky = "w")
    columns <- tcltk::tklistbox(panel, exportselection = FALSE, height = 16, width = 22)
    for (name in names(data)) tcltk::tkinsert(columns, "end", name)
    tcltk::tkgrid(columns, row = 1, column = 0, sticky = "nsw", padx = 4, pady = 6)
    model_text <- tcltk::tktext(panel, width = 65, height = 17, wrap = "none", undo = TRUE)
    tcltk::tkinsert(model_text, "1.0", if (is.null(model)) "# 回帰パス：Y ~ X\n# 潜在変数：F =~ x1 + x2 + x3\n# 共分散：X ~~ Z\n# 定義した効果：indirect := a*b\n" else paste(model, collapse = "\n"))
    tcltk::tkgrid(model_text, row = 1, column = 1, sticky = "nsew", padx = 4, pady = 6)
    scrollbar <- tcltk::ttkscrollbar(panel, orient = "vertical", command = function(...) tcltk::tkyview(model_text, ...))
    tcltk::tkconfigure(model_text, yscrollcommand = function(...) tcltk::tkset(scrollbar, ...))
    tcltk::tkgrid(scrollbar, row = 1, column = 2, sticky = "ns")
    tcltk::tkgrid(tcltk::tklabel(panel, text = "単一集団の個票データを使用します。パス・因子負荷・共分散・定義した効果をモデル式で指定できます。",
      wraplength = 740, justify = "left", anchor = "w"), row = 2, column = 0, columnspan = 3, sticky = "w", pady = 8)
    maps <- sem_choices();maps$adjust <- setNames(as.list(stats::p.adjust.methods), stats::p.adjust.methods)
    vars <- fields <- list()
    field <- function(panel, row, name, label, check = FALSE) {
      value <- initial[[name]]
      if (name %in% names(maps)) {
        match <- which(vapply(as.list(maps[[name]]), function(x) isTRUE(all.equal(x, value)), logical(1)))[1L]
        value <- names(maps[[name]])[match]
      } else value <- if (check) if (value) "1" else "0" else if (is.null(value)) "" else as.character(value)
      vars[[name]] <<- tcltk::tclVar(value)
      tcltk::tkgrid(tcltk::tklabel(panel, text = label, anchor = "w"), row = row, column = 0, sticky = "w", padx = 4, pady = 6)
      fields[[name]] <<- if (name %in% names(maps)) tcltk::ttkcombobox(panel, textvariable = vars[[name]], values = names(maps[[name]]), state = "readonly", width = 44) else
        if (check) tcltk::tkcheckbutton(panel, variable = vars[[name]], text = "有効にする") else tcltk::tkentry(panel, textvariable = vars[[name]])
      tcltk::tkgrid(fields[[name]], row = row, column = 1, sticky = "ew", padx = 4, pady = 6)
    }
    panel <- panels[["推定"]]
    field(panel, 0, "estimator", "推定法")
    field(panel, 1, "missing_data", "欠測の処理")
    field(panel, 2, "se", "標準誤差の方法")
    field(panel, 3, "bootstrap", "ブートストラップ回数")
    field(panel, 4, "seed", "乱数シード（空欄＝未指定）")
    tcltk::tkgrid(tcltk::tklabel(panel, text = "順序尺度の列（任意・複数可）", anchor = "w"), row = 5, column = 0, sticky = "nw", padx = 4, pady = 6)
    ordered_box <- tcltk::tklistbox(panel, selectmode = "extended", exportselection = FALSE, height = 5)
    for (name in names(data)) tcltk::tkinsert(ordered_box, "end", name)
    for (i in match(initial$ordered, names(data), nomatch = 0L)) if (i > 0L) tcltk::tkselection.set(ordered_box, i - 1L)
    tcltk::tkgrid(ordered_box, row = 5, column = 1, sticky = "ew", padx = 4, pady = 6)
    field(panel, 6, "std_lv", "潜在変数の分散を1に固定", TRUE)
    field(panel, 7, "fixed_x", "外生共変量を固定として扱う", TRUE)
    field(panel, 8, "meanstructure", "平均構造を含める", TRUE)
    refresh <- function() {
      est <- unname(maps$estimator[[tcltk::tclvalue(vars$estimator)]])
      has_ordered <- length(as.integer(tcltk::tkcurselection(ordered_box))) > 0L
      ml <- est == "ML" || (est == "auto" && !has_ordered)
      se_values <- if (ml) names(maps$se) else names(maps$se)[1L]
      tcltk::tkconfigure(fields$se, values = se_values)
      if (!tcltk::tclvalue(vars$se) %in% se_values) tcltk::tclvalue(vars$se) <- se_values[1L]
      boot <- tcltk::tclvalue(vars$se) == names(maps$se)[2L]
      for (name in c("bootstrap", "seed")) tcltk::tkconfigure(fields[[name]], state = if (boot) "normal" else "disabled")
      missing_values <- if (est == "WLSMV" || (est == "auto" && has_ordered)) names(maps$missing_data)[1L] else names(maps$missing_data)
      tcltk::tkconfigure(fields$missing_data, values = missing_values)
      if (!tcltk::tclvalue(vars$missing_data) %in% missing_values) tcltk::tclvalue(vars$missing_data) <- missing_values[1L]
    }
    for (name in c("estimator", "se")) tcltk::tkbind(fields[[name]], "<<ComboboxSelected>>", refresh)
    tcltk::tkbind(ordered_box, "<<ListboxSelect>>", refresh)
    refresh()
    panel <- panels[["検定・表示"]]
    labels <- c(alpha = "有意水準", followup_alpha = "有意傾向の基準", adjust = "パラメータのp値の補正法",
      ci = "信頼水準", digits = "数値の桁数", p_digits = "p値の桁数")
    for (i in seq_along(labels)) field(panel, i - 1L, names(labels)[i], labels[[i]])
    field(panel, 6, "show", "コンソールに結果を表示", TRUE)
    error_var <- tcltk::tclVar("")
    tcltk::tkgrid(tcltk::tklabel(window, textvariable = error_var, foreground = "#a32020", wraplength = 800, justify = "left"), row = 1, column = 0, sticky = "ew", padx = 12)
    result <- NULL
    confirm <- function() {
      chosen <- tryCatch({
        settings <- initial
        for (name in names(vars)) {
          value <- tcltk::tclvalue(vars[[name]])
          settings[name] <- list(if (name %in% names(maps)) as.list(maps[[name]])[[value]] else
            if (name %in% c("std_lv", "fixed_x", "meanstructure", "show")) value == "1" else
            if (name == "seed" && !nzchar(trimws(value))) NULL else suppressWarnings(as.numeric(value)))
        }
        settings$ordered <- names(data)[as.integer(tcltk::tkcurselection(ordered_box)) + 1L]
        syntax <- tcltk::tclvalue(tcltk::tkget(model_text, "1.0", "end-1c"))
        invisible(prepare_sem(data, syntax, settings))
        list(model = syntax, settings = settings)
      }, error = function(e) {tcltk::tclvalue(error_var) <- conditionMessage(e); NULL})
      if (is.null(chosen)) return(invisible(NULL))
      result <<- chosen;tcltk::tkdestroy(window)
    }
    buttons <- tcltk::tkframe(window)
    tcltk::tkgrid(buttons, row = 2, column = 0, sticky = "e", padx = 12, pady = 10)
    cancel <- function() tcltk::tkdestroy(window)
    tcltk::tkpack(tcltk::tkbutton(buttons, text = "キャンセル", command = cancel), side = "left", padx = 4)
    tcltk::tkpack(tcltk::tkbutton(buttons, text = "実行", command = confirm), side = "left", padx = 4)
    tcltk::tkwm.protocol(window, "WM_DELETE_WINDOW", cancel);tcltk::tkbind(window, "<Escape>", cancel)
    if (!as.logical(as.integer(tcltk::tkwinfo("viewable", window)))) tcltk::tkwait.visibility(window)
    tcltk::tkgrab.set(window);tcltk::tkwait.window(window)
    if (is.null(result)) stop("設定をキャンセルしました。SEMは実行していません。", call. = FALSE)
    result
  }

  choose <- function(data, formula, model, initial) {
    frame <- if (!is.null(model)) {
      if (inherits(model, "afex_aov")) model$data$long else if (inherits(model, "lm")) stats::model.frame(model) else NULL
    } else data
    if (!is.data.frame(frame)) stop("GUIにはデータフレーム、または対応する既存モデルが必要です。", call. = FALSE)
    if (!isTRUE(capabilities("tcltk")) || !suppressWarnings(requireNamespace("tcltk", quietly = TRUE))) {
      stop("GUIにはRのtcltkが必要です。引数を指定しselect = \"arguments\"で実行してください。", call. = FALSE)
    }
    window <- tryCatch(tcltk::tktoplevel(), error = function(e) {
      stop("GUIを開けませんでした。デスクトップ版のR/RStudioで実行してください。\n", conditionMessage(e), call. = FALSE)
    })
    on.exit({try(tcltk::tkgrab.release(window), silent = TRUE); try(tcltk::tkdestroy(window), silent = TRUE)}, add = TRUE)
    tcltk::tkwm.title(window, "関係性分析：回帰分析の設定")
    tcltk::tkwm.minsize(window, 740, 580)
    tcltk::tkgrid.columnconfigure(window, 0, weight = 1)
    tcltk::tkgrid.rowconfigure(window, 0, weight = 1)
    notebook <- tcltk::ttknotebook(window)
    tcltk::tkgrid(notebook, row = 0, column = 0, sticky = "nsew", padx = 12, pady = 10)
    panels <- lapply(c("変数・モデル", "検定・表示"), function(name) {
      panel <- tcltk::ttkframe(notebook, padding = 12)
      tcltk::tkadd(notebook, panel, text = name)
      tcltk::tkgrid.columnconfigure(panel, 1, weight = 1)
      panel
    })
    panel <- panels[[1L]]
    editable <- is.null(formula) && is.null(model)
    columns <- names(frame)
    numeric <- columns[vapply(frame, is.numeric, logical(1))]
    factors <- columns[vapply(frame, function(x) is.factor(x) || is.character(x), logical(1))]
    listbox <- function(label, row, values, multiple, selected) {
      tcltk::tkgrid(tcltk::tklabel(panel, text = label, anchor = "w"), row = row, column = 0, sticky = "nw", pady = 6)
      box <- tcltk::tklistbox(panel, selectmode = if (multiple) "extended" else "browse",
        exportselection = FALSE, height = if (multiple) 5 else 3, width = 35)
      for (value in values) tcltk::tkinsert(box, "end", value)
      for (index in match(selected, values, nomatch = 0L)) if (index > 0L) tcltk::tkselection.set(box, index - 1L)
      tcltk::tkgrid(box, row = row, column = 1, sticky = "ew", pady = 6)
      list(widget = box, values = values)
    }
    if (editable) {
      dv_box <- listbox("目的変数（数値・1つ）", 0, numeric, FALSE, initial$dv)
      predictor_box <- listbox("説明変数（Ctrl / Shiftで複数選択）", 1, columns, TRUE, initial$predictors)
    } else {
      tcltk::tkgrid(tcltk::tklabel(panel, text = if (!is.null(model)) "推定済みモデルを使用します。" else
        paste(deparse(formula), collapse = " "), anchor = "w", justify = "left", wraplength = 650),
        row = 0, column = 0, columnspan = 2, sticky = "w", pady = 8)
    }
    interaction_var <- tcltk::tclVar(if (initial$interactions) "1" else "0")
    if (editable) tcltk::tkgrid(tcltk::tkcheckbutton(panel, text = "説明変数同士の交互作用を含める", variable = interaction_var),
      row = 2, column = 1, sticky = "w", pady = 4)
    trend_var <- tcltk::tclVar(if (is.null(initial$trend)) "（指定なし）" else initial$trend)
    tcltk::tkgrid(tcltk::tklabel(panel, text = "傾きを調べる数値変数（任意）", anchor = "w"), row = 3, column = 0, sticky = "w", pady = 8)
    tcltk::tkgrid(tcltk::ttkcombobox(panel, textvariable = trend_var,
      values = c("（指定なし）", numeric), state = "readonly"), row = 3, column = 1, sticky = "ew", pady = 8)
    by_box <- listbox("条件別に調べる因子（任意・複数可）", 4, factors, TRUE, initial$by)
    compare_var <- tcltk::tclVar(if (initial$compare_slopes) "1" else "0")
    tcltk::tkgrid(tcltk::tkcheckbutton(panel, text = "条件間の傾きも比較する", variable = compare_var), row = 5, column = 1, sticky = "w")
    variables <- list()
    fields <- c(alpha = "有意水準", followup_alpha = "有意傾向・比較の表示基準", ci = "信頼水準", digits = "数値の小数点以下の桁数", p_digits = "p値の小数点以下の桁数")
    for (i in seq_along(fields)) {
      name <- names(fields)[i]
      variables[[name]] <- tcltk::tclVar(as.character(initial[[name]]))
      tcltk::tkgrid(tcltk::tklabel(panels[[2L]], text = fields[[name]], anchor = "w"), row = i - 1L, column = 0, sticky = "w", padx = 4, pady = 8)
      tcltk::tkgrid(tcltk::tkentry(panels[[2L]], textvariable = variables[[name]]), row = i - 1L, column = 1, sticky = "ew", pady = 8)
    }
    adjust_var <- tcltk::tclVar(initial$adjust)
    tcltk::tkgrid(tcltk::tklabel(panels[[2L]], text = "傾き・傾き差のp値の補正法"), row = 5, column = 0, sticky = "w", pady = 8)
    tcltk::tkgrid(tcltk::ttkcombobox(panels[[2L]], textvariable = adjust_var, values = stats::p.adjust.methods,
      state = "readonly"), row = 5, column = 1, sticky = "ew", pady = 8)
    show_var <- tcltk::tclVar(if (initial$show) "1" else "0")
    tcltk::tkgrid(tcltk::tkcheckbutton(panels[[2L]], text = "コンソールに結果を表示", variable = show_var), row = 6, column = 1, sticky = "w", pady = 8)
    error_var <- tcltk::tclVar("")
    tcltk::tkgrid(tcltk::tklabel(window, textvariable = error_var, foreground = "#a32020", wraplength = 700, justify = "left"),
      row = 1, column = 0, sticky = "w", padx = 12)
    result <- NULL
    selected <- function(box) box$values[as.integer(tcltk::tkcurselection(box$widget)) + 1L]
    confirm <- function() {
      value <- tryCatch({
        settings <- initial
        for (name in names(variables)) settings[[name]] <- suppressWarnings(as.numeric(tcltk::tclvalue(variables[[name]])))
        settings$adjust <- tcltk::tclvalue(adjust_var)
        settings$show <- tcltk::tclvalue(show_var) == "1"
        settings$compare_slopes <- tcltk::tclvalue(compare_var) == "1"
        settings$interactions <- tcltk::tclvalue(interaction_var) == "1"
        trend <- tcltk::tclvalue(trend_var)
        settings["trend"] <- list(if (trend == "（指定なし）") NULL else trend)
        settings$by <- selected(by_box)
        selected_formula <- if (editable) make_formula(selected(dv_box), selected(predictor_box), settings$interactions) else formula
        validate_options(settings$alpha, settings$followup_alpha, settings$adjust, settings$ci,
          settings$digits, settings$p_digits, settings$show, settings$compare_slopes)
        rhs <- if (!is.null(selected_formula)) all.vars(selected_formula[[3L]]) else if (inherits(model, "lm"))
          all.vars(stats::delete.response(stats::terms(model))) else names(frame)
        if (length(settings$by) && is.null(settings$trend)) stop("条件別の傾きを調べる数値変数も選択してください。")
        if (any(!c(settings$trend, settings$by) %in% rhs)) stop("傾き・条件にはモデルの説明変数を選択してください。")
        if (inherits(model, "afex_aov") && is.null(settings$trend)) stop("既存ANCOVAモデルでは、傾きを調べる数値変数を選択してください。")
        list(formula = selected_formula, settings = settings)
      }, error = function(e) {tcltk::tclvalue(error_var) <- conditionMessage(e); NULL})
      if (is.null(value)) return(invisible(NULL))
      result <<- value
      tcltk::tkdestroy(window)
    }
    cancel <- function() tcltk::tkdestroy(window)
    buttons <- tcltk::tkframe(window)
    tcltk::tkgrid(buttons, row = 2, column = 0, sticky = "e", padx = 12, pady = 10)
    tcltk::tkpack(tcltk::tkbutton(buttons, text = "キャンセル", command = cancel), side = "left", padx = 4)
    tcltk::tkpack(tcltk::tkbutton(buttons, text = "実行", command = confirm), side = "left", padx = 4)
    tcltk::tkwm.protocol(window, "WM_DELETE_WINDOW", cancel)
    tcltk::tkbind(window, "<Escape>", cancel)
    if (!as.logical(as.integer(tcltk::tkwinfo("viewable", window)))) tcltk::tkwait.visibility(window)
    tcltk::tkgrab.set(window)
    tcltk::tkwait.window(window)
    if (is.null(result)) stop("設定をキャンセルしました。回帰分析は実行していません。", call. = FALSE)
    result
  }
}, envir = .relation_private_environment)

relate <- function(data = NULL, formula = NULL, model = NULL, dv = NULL, predictors = NULL,
                   interactions = FALSE, trend = NULL, by = NULL, compare_slopes = TRUE,
                   alpha = .05, followup_alpha = .10, adjust = "bonferroni", ci = .95,
                   digits = 3, p_digits = 3, show = TRUE,
                   select = if (is.null(formula) && is.null(model) && (is.null(dv) || is.null(predictors))) "gui" else "arguments",
                   result_rows = NULL, analysis = NULL, variables = NULL,
                   method = "pearson", use = "pairwise.complete.obs",
                   alternative = "two.sided", exact = NULL,
                   sem_model = NULL, estimator = "auto", missing_data = "listwise", ordered = NULL,
                   std_lv = FALSE, fixed_x = TRUE, meanstructure = TRUE, se = "default",
                   bootstrap = 1000, seed = NULL) {
  internal <- .relation_private_environment
  if (is.null(analysis)) analysis <- if (!is.null(sem_model) || is.character(model)) "sem" else if (!is.null(variables) || !missing(method)) "correlation" else "regression"
  analysis <- match.arg(analysis, c("regression", "correlation", "sem"))
  if (analysis == "sem") {
    if (internal$is_result_input(data) || !is.null(formula) || length(dv) || length(predictors) || length(trend) ||
        length(by) || !identical(interactions, FALSE) || !missing(compare_slopes) || !is.null(result_rows) ||
        !is.null(variables) || !missing(method) || !missing(use) || !missing(alternative) || !missing(exact)) {
      stop("SEMにはデータとsem_modelを指定してください。相関・回帰用の引数やANCOVA結果は指定できません。", call. = FALSE)
    }
    if (!is.null(model)) {
      if (!is.null(sem_model) || !is.character(model)) stop("SEMのモデル式はsem_modelまたは文字列のmodelの一方で指定してください。", call. = FALSE)
      sem_model <- model
    }
    if (missing(adjust)) adjust <- "none"
    if (missing(select)) select <- if (is.null(sem_model)) "gui" else "arguments"
    select <- match.arg(select, c("gui", "arguments"))
    internal$validate_options(alpha, followup_alpha, adjust, ci, digits, p_digits, show, TRUE)
    settings <- list(estimator = estimator, missing_data = missing_data, ordered = ordered, std_lv = std_lv,
      fixed_x = fixed_x, meanstructure = meanstructure, se = se, bootstrap = bootstrap, seed = seed,
      alpha = alpha, followup_alpha = followup_alpha, adjust = adjust, ci = ci, digits = digits, p_digits = p_digits, show = show)
    return(invisible(internal$run_sem(data, sem_model, settings, select, substitute(data))))
  }
  if (!is.null(sem_model) || !missing(estimator) || !missing(missing_data) || !missing(ordered) ||
      !missing(std_lv) || !missing(fixed_x) || !missing(meanstructure) || !missing(se) || !missing(bootstrap) || !missing(seed)) {
    stop('SEM用の引数はanalysis = "sem"で指定してください。', call. = FALSE)
  }
  if (analysis == "correlation") {
    if (internal$is_result_input(data) || !is.null(formula) || !is.null(model) || length(dv) ||
        length(predictors) || length(trend) || length(by) || !identical(interactions, FALSE) ||
        !missing(compare_slopes) || !is.null(result_rows)) {
      stop("相関分析にはデータとvariablesを指定してください。回帰用の引数やANCOVA結果は指定できません。", call. = FALSE)
    }
    if (missing(select)) select <- if (is.null(variables)) "gui" else "arguments"
    select <- match.arg(select, c("gui", "arguments"))
    internal$validate_options(alpha, followup_alpha, adjust, ci, digits, p_digits, show, TRUE)
    method <- match.arg(method, c("pearson", "spearman", "kendall"))
    use <- match.arg(use, c("pairwise.complete.obs", "complete.obs"))
    alternative <- match.arg(alternative, c("two.sided", "greater", "less"))
    if (!is.null(exact) && (!is.logical(exact) || length(exact) != 1L || is.na(exact))) stop("exactにはNULL、TRUE、FALSEを指定してください。", call. = FALSE)
    settings <- list(method = method, use = use, alternative = alternative, exact = exact,
      alpha = alpha, followup_alpha = followup_alpha, adjust = adjust, ci = ci,
      digits = digits, p_digits = p_digits, show = show)
    return(invisible(internal$correlate(data, variables, settings, select, substitute(data))))
  }
  if (!is.null(variables) || !missing(method) || !missing(use) || !missing(alternative) || !missing(exact)) {
    stop('variables / method / use / alternative / exactはanalysis = "correlation"で指定してください。', call. = FALSE)
  }
  inherit <- c(alpha = missing(alpha), followup_alpha = missing(followup_alpha),
    adjust = missing(adjust), digits = missing(digits), p_digits = missing(p_digits))
  from_results <- internal$is_result_input(data)
  if (missing(select) && from_results) select <- "arguments"
  select <- match.arg(select, c("gui", "arguments"))
  internal$validate_options(alpha, followup_alpha, adjust, ci, digits, p_digits, show, compare_slopes)
  if (!is.logical(interactions) || length(interactions) != 1L || is.na(interactions)) stop("interactionsにはTRUEまたはFALSEを指定してください。", call. = FALSE)
  data_expression <- substitute(data)
  model_expression <- substitute(model)
  settings <- list(dv = dv, predictors = predictors, interactions = interactions, trend = trend, by = by,
    compare_slopes = compare_slopes, alpha = alpha, followup_alpha = followup_alpha,
    adjust = adjust, ci = ci, digits = digits, p_digits = p_digits, show = show)
  if (from_results) {
    if (!is.null(formula) || !is.null(model) || length(dv) || length(predictors) || length(trend) || length(by) || interactions) {
      stop("結果からの分析ではモデル・共変量・条件を自動取得します。変数を指定する場合はデータまたはmodelを渡してください。", call. = FALSE)
    }
    output <- internal$followup_results(data, settings, inherit, select, result_rows)
    replay <- output$settings
    gui_keys <- if (select == "gui") c("ci", "digits", "p_digits", "adjust", "compare_slopes", "show") else character()
    keys <- setdiff(names(replay), c("dv", "predictors", "interactions", "trend", "by"))
    keys <- keys[!keys %in% names(inherit)[inherit] | keys %in% gui_keys]
    args <- c(list(data = data_expression), replay[keys], list(select = "arguments"))
    if (select == "gui" || !is.null(result_rows)) args$result_rows <- output$selected_rows
    call <- as.call(c(list(as.name("relate")), args))
    code <- deparse(call, width.cutoff = 80L)
    code[1L] <- paste0("relation_result <- ", code[1L])
    output$selection <- list(input_type = "results", method = select, call = call, code = paste(code, collapse = "\n"))
    if (output$show) {
      internal$display_followups(output, output$settings$digits, output$settings$p_digits)
      if (select == "gui") cat("\n【再実行用コード】\n", output$selection$code, "\n", sep = "")
    }
    return(invisible(output))
  }
  if (!is.null(result_rows)) stop("result_rowsはANCOVA結果からの分析でのみ指定できます。", call. = FALSE)
  if (!is.null(model) && (!is.null(formula) || length(dv) || length(predictors))) stop("modelと、formula / dv / predictorsは同時に指定できません。", call. = FALSE)
  if (select == "gui") {
    chosen <- internal$choose(data, formula, model, settings)
    formula <- chosen$formula
    settings <- chosen$settings
  } else if (is.null(model) && is.null(formula)) formula <- internal$make_formula(dv, predictors, interactions)
  result <- internal$run_regression(data, formula, model, settings)
  args <- c(if (!is.null(data)) list(data = data_expression), if (!is.null(model)) list(model = model_expression) else list(formula = formula),
    settings[setdiff(names(settings), c("dv", "predictors", "interactions"))], list(select = "arguments"))
  call <- as.call(c(list(as.name("relate")), args))
  code <- deparse(call, width.cutoff = 80L)
  code[1L] <- paste0("result <- ", code[1L])
  result$selection <- list(method = select, call = call, code = paste(code, collapse = "\n"))
  class(result) <- unique(c("analysis_relation_result", class(result)))
  if (settings$show) {
    internal$display(result, settings$digits, settings$p_digits)
    if (select == "gui") cat("\n【再実行用コード】\n", result$selection$code, "\n", sep = "")
  }
  invisible(result)
}

if (!isTRUE(getOption("analysis.helper.quiet", FALSE)) && !isTRUE(getOption("analysis.load.quiet", FALSE))) {
  cat("\n========================================\nrelation helper：依存パッケージ確認\n========================================\n",
    "通常の線形回帰：R標準機能を使用\n",
    "相関分析（Pearson / Spearman / Kendall）：R標準機能を使用\n",
    "SEM：lavaan ", if (requireNamespace("lavaan", quietly = TRUE)) "OK\n" else "未導入（利用時に必要です）\n",
    "条件別の傾き・傾き差：emmeans ",
    if (requireNamespace("emmeans", quietly = TRUE)) "OK\n" else "未導入（利用時に必要です）\n",
    "置換ANCOVA後の傾き・傾き差：permuco ",
    if (requireNamespace("permuco", quietly = TRUE)) "OK\n" else "未導入（利用時に必要です）\n", sep = "")
}
rm(.relation_source, .relation_directory, .relation_child, .relation_frame, .relation_script)
