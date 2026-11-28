# ============================================================
# 検定・分析の統合インターフェース
#
# 公開関数：
#   analyze()
#
# 内部関数：
#   t_anova()
#   t_ancova()
#   t_art()
#   t_permuco()
#   t_rank()
#   t_brunner_munzel()
#   t_ats()
#   t_t()
#   t_welch()
#   t_rank_ancova()
#
# 下位関数は専用環境へ読み込み、RStudioのFunctionsには公開しません。
# ============================================================


# ------------------------------------------------------------
# このhelper自身の場所を取得
# ------------------------------------------------------------

.test_source_file <- NULL
.test_frames <- sys.frames()

for (.test_index in rev(seq_along(.test_frames))) {
  .test_candidate <- .test_frames[[.test_index]]$ofile

  if (!is.null(.test_candidate) && nzchar(.test_candidate)) {
    .test_source_file <- .test_candidate
    break
  }
}

if (is.null(.test_source_file)) {
  stop(
    "test_helper.Rの場所を取得できませんでした。\n",
    'コンソールからsource("script/test/test_helper.R", encoding = "UTF-8")を実行してください。',
    call. = FALSE
  )
}

.test_directory <- dirname(
  normalizePath(
    .test_source_file,
    winslash = "/",
    mustWork = TRUE
  )
)

.test_child_paths <- file.path(
  .test_directory,
  c(
    "t_anova.R",
    "t_ancova.R",
    "t_art.R",
    "permuco_followups.R",
    "t_permuco.R",
    "t_rank.R",
    "t_brunner_munzel.R",
    "t_ats.R",
    "t_t.R",
    "t_welch.R",
    "t_rank_ancova.R"
  )
)

.test_missing_scripts <- .test_child_paths[
  !file.exists(.test_child_paths)
]

if (length(.test_missing_scripts) > 0) {
  stop(
    "次の検定・分析スクリプトが見つかりません：\n",
    paste0("  - ", .test_missing_scripts, collapse = "\n"),
    "\n\ntest_helper.Rと10本のt_スクリプトを同じフォルダに配置してください。",
    call. = FALSE
  )
}


# ------------------------------------------------------------
# 下位関数を専用環境へ読み込む
# ------------------------------------------------------------

.test_private_environment <- new.env(parent = environment())

for (.test_path in .test_child_paths) {
  source(
    .test_path,
    local = .test_private_environment,
    encoding = "UTF-8",
    chdir = FALSE
  )
}

.test_required_functions <- c(
  "t_anova",
  "t_ancova",
  "t_art",
  "t_permuco",
  "t_rank",
  "t_brunner_munzel",
  "t_ats",
  "t_t",
  "t_welch",
  "t_rank_ancova"
)

.test_loaded <- vapply(
  .test_required_functions,
  exists,
  logical(1),
  envir = .test_private_environment,
  mode = "function",
  inherits = FALSE
)

if (any(!.test_loaded)) {
  stop(
    "次の検定・分析関数を読み込めませんでした：\n",
    paste0(
      "  - ",
      .test_required_functions[!.test_loaded],
      "()",
      collapse = "\n"
    ),
    call. = FALSE
  )
}


# ------------------------------------------------------------
# test系で使用する外部パッケージを読み込み時に確認
#
# ここでは不足を通知するだけでhelper自体は停止しない。
# 各t_*()は実行時にも必要パッケージを確認し、
# 実際に必要な分析を実行するときの最終安全確認を行う。
#
# options(analysis.helper.quiet = TRUE)
# またはoptions(analysis.load.quiet = TRUE)で表示を抑制できる。
# ------------------------------------------------------------

.test_package_groups <- list(
  "ANOVA / ANCOVA" = c("afex", "emmeans", "car"),
  "置換ANOVA / ANCOVA" = "permuco",
  "ART" = c("ARTool", "emmeans", "car"),
  "ランク（Dunn下位検定のみ）" = "PMCMRplus",
  "Brunner–Munzel" = "lawstat",
  "ATS＋MCTP" = "rankFD",
  "Welch（Games–Howell下位検定のみ）" = "PMCMRplus",
  "順位ANCOVA（Quade）" = "npANCOVA"
)

.test_package_status <- lapply(
  .test_package_groups,
  function(packages) {
    installed <- vapply(
      packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )

    list(
      installed = packages[installed],
      missing = packages[!installed]
    )
  }
)

.test_missing_packages <- unique(
  unlist(
    lapply(.test_package_status, `[[`, "missing"),
    use.names = FALSE
  )
)

.test_helper_quiet <- isTRUE(getOption("analysis.helper.quiet", FALSE)) ||
  isTRUE(getOption("analysis.load.quiet", FALSE))

.test_group_name <- NULL
.test_group_missing <- NULL

if (!.test_helper_quiet) {
  cat(
    "\n========================================\n",
    "test helper：依存パッケージ確認\n",
    "========================================\n",
    sep = ""
  )

  for (.test_group_name in names(.test_package_groups)) {
    .test_group_missing <- .test_package_status[[.test_group_name]]$missing

    if (length(.test_group_missing) == 0) {
      cat("  ", .test_group_name, "：OK\n", sep = "")
    } else {
      cat(
        "  ",
        .test_group_name,
        "：不足 -> ",
        paste(.test_group_missing, collapse = ", "),
        "\n",
        sep = ""
      )
    }
  }

  if (length(.test_missing_packages) > 0) {
    cat(
      "\n不足パッケージがあってもhelperは読み込みます。\n",
      "該当する分析の実行時に再確認し、不足していれば停止します。\n",
      sep = ""
    )
  } else {
    cat("\n必要な外部パッケージはすべてインストール済みです。\n")
  }
}


# 共通処理は下位関数と同じ専用環境に保持する。
evalq({
  significance <- function(p) {
    ifelse(is.na(p), "", ifelse(p < .001, "***", ifelse(p < .01, "**",
      ifelse(p < .05, "*", ifelse(p < .10, ".", "")))))
  }

  select_followups <- function(table, alpha, trend_alpha, simple_effects,
                               posthoc = TRUE, followup_scope = "followup") {
    effects <- table$Effect
    p <- table$p.value
    valid <- is.finite(p) & !grepl("Residual|Error", effects, ignore.case = TRUE)
    split_effect <- function(x) strsplit(gsub("`", "", x, fixed = TRUE), ":", fixed = TRUE)[[1]]
    cutoff <- if (followup_scope == "significant") alpha else if (followup_scope == "all") Inf else trend_alpha
    selected <- if (posthoc) effects[valid & p < cutoff] else character()
    main <- selected[lengths(lapply(selected, split_effect)) == 1L]
    forced <- character()
    if (posthoc && simple_effects && length(main)) {
      pairs <- effects[valid & lengths(lapply(effects, split_effect)) == 2L]
      forced <- pairs[!pairs %in% selected & vapply(pairs, function(x) any(split_effect(x) %in% main), logical(1))]
    }
    list(significant = effects[valid & p < alpha],
      trend = effects[valid & p >= alpha & p < trend_alpha],
      main = main, forced = forced, effects = effects[effects %in% c(selected, forced)])
  }

  followup_scope_label <- function(scope) switch(scope, significant = "有意な効果のみ", followup = "有意＋有意傾向の効果", all = "すべて")
  display_scope_label <- function(scope) switch(scope, significant = "有意な比較のみ", effect = "有意な効果の全比較", all = "計算済みの全結果")

  # 表示だけの判定。計算や補正済みの値を変更しない。
  visible_followup_rows <- function(result, effect, p) {
    scope <- if (is.null(result$display_scope)) "significant" else result$display_scope
    switch(scope, significant = is.finite(p) & p < result$alpha,
      effect = rep(effect %in% result$significant_effects, length.out = length(p)), all = rep(TRUE, length(p)))
  }

  apply_sphericity <- function(table, sphere, correction, alpha = .05) {
    table$p.uncorrected <- table$p.value
    table$df1.uncorrected <- table$df1
    table$df2.uncorrected <- table$df2
    table$correction <- "none"
    rows <- match(table$Effect, sphere$Effect)
    within <- !is.na(rows) & nzchar(sphere$Within[rows])
    within_df <- if ("Within.df" %in% names(sphere)) sphere$Within.df[rows] else table$df1
    applicable <- within & is.finite(within_df) & within_df > 1 + sqrt(.Machine$double.eps)
    decisions <- data.frame(Effect = table$Effect, requested = correction,
      Mauchly.p = sphere$Mauchly.p[rows], GG.estimate = sphere$GG[rows], HF.estimate = sphere$HF[rows],
      correction = "none", epsilon = 1, selection_rule = "none", epsilon_threshold = .75,
      status = ifelse(applicable, "not_requested", "not_applicable"), stringsAsFactors = FALSE)
    notices <- character()
    for (i in which(applicable)) {
      selected <- correction
      if (correction == "none") next
      if (correction == "auto") {
        p <- decisions$Mauchly.p[i]
        if (is.finite(p) && p >= alpha) {
          decisions$status[i] <- "sphericity_not_rejected"
          decisions$selection_rule[i] <- "mauchly_not_rejected"
          next
        }
        gg <- decisions$GG.estimate[i]
        valid_gg <- is.finite(gg) && gg > 0 && gg <= 1 + sqrt(.Machine$double.eps)
        selected <- if (!valid_gg) "unavailable" else if (gg <= .75) "GG" else "HF"
        decisions$selection_rule[i] <- if (!valid_gg) "GG_epsilon_unavailable" else
          if (selected == "GG") "GG_epsilon_le_0.75" else "GG_epsilon_gt_0.75"
        decisions$status[i] <- if (is.finite(p)) "sphericity_rejected" else "mauchly_unavailable_epsilon_selection"
        if (!is.finite(p)) notices <- c(notices,
          paste0("球面性の自動判定不可 [", table$Effect[i], "]：GGのεによるGG/HF選択を行います。"))
      } else {
        decisions$status[i] <- "explicit_correction"
        decisions$selection_rule[i] <- "explicit_correction"
      }
      epsilon <- if (selected == "unavailable") NA_real_ else sphere[[selected]][rows[i]]
      if (length(epsilon) != 1L || !is.finite(epsilon) || epsilon <= 0 || epsilon > 1 + sqrt(.Machine$double.eps)) {
        table$df1[i] <- table$df2[i] <- table$p.value[i] <- NA_real_
        table$correction[i] <- decisions$correction[i] <- "unavailable"
        decisions$epsilon[i] <- NA_real_
        decisions$status[i] <- "epsilon_unavailable"
        notices <- c(notices, paste0("球面性補正のεを推定できません [", table$Effect[i], "]：補正後の自由度・pは算出不可です。"))
        next
      }
      epsilon <- min(1, epsilon)
      if ("MSE" %in% names(table)) table$MSE[i] <- table$MSE[i] / epsilon
      table$df1[i] <- table$df1[i] * epsilon
      table$df2[i] <- table$df2[i] * epsilon
      table$p.value[i] <- stats::pf(table$F[i], table$df1[i], table$df2[i], lower.tail = FALSE)
      table$correction[i] <- decisions$correction[i] <- selected
      decisions$epsilon[i] <- epsilon
    }
    attr(table, "sphericity_decisions") <- decisions
    attr(table, "sphericity_messages") <- unique(notices)
    table
  }

  # 共通補正の対応表。固有補正・下位検定なしの方法は共通設定を適用しない。
  adjustment_policy <- function(methods) {
    p_adjust <- c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "none")
    emmeans_adjust <- c("bonferroni", "holm", "sidak", "tukey", "hochberg", "hommel", "BH", "BY", "none")
    configurable <- intersect(methods, c("anova", "ancova", "art", "rank", "rank_ancova", "permuco", "permuco_ancova"))
    choices <- emmeans_adjust
    if (any(configurable %in% c("rank", "rank_ancova", "permuco", "permuco_ancova"))) choices <- intersect(choices, p_adjust)
    list(configurable = configurable, choices = choices,
      emmeans = intersect(methods, c("anova", "ancova", "art")),
      p_adjust = intersect(methods, c("rank", "rank_ancova", "permuco", "permuco_ancova")),
      intrinsic = intersect(methods, c("welch", "ats")),
      inactive = intersect(methods, c("t", "brunner_munzel")))
  }

  normalize_adjustment <- function(adjust) {
    if (is.character(adjust) && length(adjust) == 1L && !is.na(adjust) && adjust == "fdr") "BH" else adjust
  }

  validate_adjustment <- function(methods, adjust) {
    adjust <- normalize_adjustment(adjust)
    policy <- adjustment_policy(methods)
    if (!is.character(adjust) || length(adjust) != 1L || is.na(adjust) || !adjust %in% policy$choices) {
      stop("選択した方法の共通補正には次から指定してください: ", paste(policy$choices, collapse = ", "),
        "。Games–Howell/MCTPは固有の調整を使います。", call. = FALSE)
    }
    invisible(adjust)
  }

  adjust_contrasts <- function(contrast, adjust, family, confidence_level = .95) {
    contrast$requested_adjust <- adjust
    if (!is.null(contrast$error)) return(contrast)
    if (isTRUE(contrast$kind %in% c("mctp", "games_howell"))) {
      # 全ペアを考慮したp値・同時信頼区間が渡される。二重補正しない。
      table <- contrast$formatted
      table$significance <- significance(table$p.value)
      table$adjust <- contrast$intrinsic_adjustment
      table$family <- family
      table$family_size <- nrow(table)
      contrast$formatted <- table
      contrast$adjust <- contrast$intrinsic_adjustment
      contrast$p_adjustment <- contrast$ci_adjustment <- contrast$intrinsic_adjustment
      contrast$additional_adjustment <- "none"
      table$p_adjustment <- table$ci_adjustment <- contrast$intrinsic_adjustment
      contrast$formatted <- table
      contrast$family <- family
      contrast$family_size <- nrow(table)
      contrast$adjustment_messages <- paste0("Confidence level used: ", contrast$confidence_level)
      return(contrast)
    }
    original <- contrast$unadjusted
    if (isTRUE(contrast$kind %in% c("rank", "rank_ancova", "permuco"))) {
      # 全ペアを1つのfamilyとし、抽出・表示より先に補正する。
      original$p.uncorrected <- original$p.value
      original$p.value <- stats::p.adjust(original$p.value, method = adjust, n = nrow(original))
      original$significance <- significance(original$p.value)
      original$adjust <- adjust
      original$family <- family
      original$family_size <- nrow(original)
      contrast$formatted <- original
      contrast$adjust <- adjust
      contrast$p_adjustment <- adjust
      contrast$ci_adjustment <- "not_computed"
      original$p_adjustment <- adjust
      original$ci_adjustment <- "not_computed"
      contrast$formatted <- original
      contrast$family <- family
      contrast$family_size <- nrow(original)
      return(contrast)
    }
    # emmGridは比較間の共分散を保持。Tukey等もemmeansの標準処理に委譲する。
    adjusted_summary <- summary(contrast$comparison_grid, adjust = adjust, infer = c(TRUE, TRUE), level = confidence_level)
    corrected <- as.data.frame(adjusted_summary)
    if (nrow(corrected) != nrow(original)) stop("補正前後の比較行が一致しません。")
    original$p.uncorrected <- original$p.value
    original$p.value <- corrected$p.value
    for (name in intersect(c("lower.CL", "upper.CL", "asymp.LCL", "asymp.UCL"), names(corrected))) {
      original[[name]] <- corrected[[name]]
    }
    original$significance <- significance(original$p.value)
    original$adjust <- adjust
    by_vars <- attr(adjusted_summary, "by.vars")
    family_labels <- rep(family, nrow(original))
    if (length(by_vars)) {
      conditions <- apply(original[, by_vars, drop = FALSE], 1, function(row) {
        paste(paste0(by_vars, "=", row), collapse = ", ")
      })
      family_labels <- paste(family, conditions, sep = " | ")
    }
    original$family <- family_labels
    counts <- table(family_labels)
    original$family_size <- as.integer(counts[family_labels])
    contrast$formatted <- original
    contrast$adjust <- adjust
    contrast$adjustment_messages <- attr(adjusted_summary, "mesg")
    contrast$confidence_level <- confidence_level
    contrast$p_adjustment <- sub(" method.*$", "", sub("^P value adjustment: *", "",
      contrast$adjustment_messages[startsWith(contrast$adjustment_messages, "P value adjustment:")]))
    contrast$ci_adjustment <- sub(" method.*$", "", sub("^Conf-level adjustment: *", "",
      contrast$adjustment_messages[startsWith(contrast$adjustment_messages, "Conf-level adjustment:")]))
    if (!length(contrast$p_adjustment)) contrast$p_adjustment <- if (all(original$family_size == 1L)) "none" else adjust
    if (!length(contrast$ci_adjustment)) contrast$ci_adjustment <- if (all(original$family_size == 1L)) "none" else
      if (adjust %in% c("holm", "hochberg", "hommel", "BH", "BY")) "bonferroni" else adjust
    contrast$p_adjustment <- unique(contrast$p_adjustment)
    contrast$ci_adjustment <- unique(contrast$ci_adjustment)
    original$p_adjustment <- paste(contrast$p_adjustment, collapse = "/")
    original$ci_adjustment <- paste(contrast$ci_adjustment, collapse = "/")
    original$adjust <- original$p_adjustment
    contrast$formatted <- original
    contrast$adjust <- original$p_adjustment[1L]
    contrast$family <- names(counts)
    contrast$family_size <- as.integer(counts)
    contrast
  }

  prepare_result <- function(result, method, alpha, trend_alpha,
                             simple_effects, sphericity, adjust, correction,
                             posthoc = TRUE, followup_scope = "followup", display_scope = "significant") {
    direct_method <- method %in% c("rank", "brunner_munzel", "ats", "t", "welch", "rank_ancova")
    table <- if (direct_method) result$test_table else result$anova_table
    if (method %in% c("anova", "ancova", "art")) {
      selected <- if (method == "art" && !is.null(sphericity)) sphericity else correction
      # 群間要因との交互作用でも、球面性の必要性は被験者内の対比次元で決まる。
      # 例：3群×2時点は交互作用df1=2でも、被験者内の対比は1つなので対象外。
      result$sphericity_table$Within.df <- vapply(result$sphericity_table$Effect, function(effect) {
        terms <- strsplit(gsub("`", "", effect, fixed = TRUE), ":", fixed = TRUE)[[1L]]
        factors <- intersect(terms, result$within)
        if (!length(factors)) return(0)
        prod(vapply(result$analysis_data[factors], function(x) nlevels(factor(x)) - 1L, integer(1)))
      }, numeric(1))
      table <- apply_sphericity(result$anova_table_uncorrected, result$sphericity_table, selected, alpha)
      result$sphericity <- selected
      result$legacy_sphericity_override <- method == "art" && !is.null(sphericity)
      result$sphericity_decisions <- attr(table, "sphericity_decisions")
      result$applied_within_correction <- unique(result$sphericity_decisions$correction)
      result$messages <- unique(c(result$messages, attr(table, "sphericity_messages")))
    } else {
      table$p.uncorrected <- table$p.value
      table$correction <- if (direct_method) "none" else "permutation"
    }
    table$significance <- significance(table$p.value)
    result$test_table <- table
    if (!direct_method || method == "rank_ancova" || (method == "welch" && result$n_levels > 2L)) result$anova_table <- table
    if (method == "rank_ancova") result$covariate_table$significance <- significance(result$covariate_table$p.value)
    if (!is.null(result$levene)) {
      result$levene$significance <- significance(result$levene$p.value)
      result$levene$judgement <- ifelse(is.na(result$levene$p.value), "判定不可",
        ifelse(result$levene$p.value < alpha, "等分散性に注意", "等分散性を棄却せず"))
    }
    result$requested_adjust <- adjust
    result$adjustment_engine <- if (!is.function(result$compute_followups)) "not_applied" else
      if (method %in% c("anova", "ancova", "art")) "emmeans" else
      if (method %in% c("rank", "rank_ancova", "permuco", "permuco_ancova")) "p.adjust" else if (method %in% c("welch", "ats")) "intrinsic" else "not_applied"
    result$adjust <- if (method == "ats") "MCTP (Fisher, Tukey)" else if (method == "welch" && result$n_levels > 2L)
      "Games–Howell" else if (method %in% c("brunner_munzel", "t", "welch") ||
        (method == "rank_ancova" && result$n_levels == 2L)) "none" else adjust
    result$trend_alpha <- trend_alpha
    result$followup_alpha <- trend_alpha # 旧結果との互換名。同じ共通基準を保存する。
    result$posthoc <- posthoc
    if (!posthoc) {
      result$adjustment_engine <- "not_applied"
      result$adjust <- "none"
    }
    result$followup_scope <- followup_scope
    result$display_scope <- display_scope
    selection <- select_followups(table, alpha, trend_alpha, FALSE, posthoc, followup_scope)
    result$significant_effects <- selection$significant
    result$trend_effects <- selection$trend
    # 共変量とその交互作用は、実験要因の水準間比較には含めない。
    eligible <- rep(TRUE, nrow(table))
    if (!direct_method && !is.null(result$followup_factors)) eligible <- vapply(table$Effect, function(effect) {
      all(strsplit(gsub("`", "", effect, fixed = TRUE), ":", fixed = TRUE)[[1L]] %in% result$followup_factors)
    }, logical(1))
    selection <- select_followups(table[eligible, , drop = FALSE], alpha, trend_alpha,
      simple_effects && !direct_method, posthoc, followup_scope)
    result$forced_simple_effects <- selection$forced
    if (method %in% c('permuco','permuco_ancova')) {
      selection$effects <- selection$effects[lengths(strsplit(selection$effects, ':', fixed=TRUE)) <= 2L]
      if (posthoc) result$messages <- unique(c(result$messages,
        '下位検定：参加者単位HC2の置換対比検定（独自拡張）。残差置換は近似的方法です。被験者内の符号反転には帰無下の対称性が必要です。信頼区間は算出しません。'))
    }
    result$followup_main_effects <- selection$main
    result$follow_up_effects <- if (is.function(result$compute_followups)) selection$effects else character()
    result
  }

  finalize_effect <- function(result, effect, method, adjust) {
    compute <- function() result$compute_followups(effect, result$forced_simple_effects, result$followup_main_effects)
    contrasts <- if (method %in% c("rank", "ats", "welch", "rank_ancova")) tryCatch(compute(), error = function(e) {
      setNames(list(list(kind = switch(method, ats = "mctp", welch = "games_howell", rank_ancova = "rank_ancova", "rank"), error = conditionMessage(e), source_effect = effect,
        error_details = list(stage = "下位検定", message = conditionMessage(e)))), effect)
    }) else compute()
    adjusted <- lapply(names(contrasts), function(name) {
      tryCatch(adjust_contrasts(contrasts[[name]], adjust, name, 1 - result$alpha), error = function(e) {
        failed <- contrasts[[name]]
        failed$error <- conditionMessage(e)
        failed$error_details <- list(stage = "多重比較補正", message = conditionMessage(e))
        failed
      })
    })
    names(adjusted) <- names(contrasts)
    result$contrasts <- c(result$contrasts, adjusted)
    if (method %in% c("anova", "ancova")) {
      comparisons <- attr(contrasts, "comparisons")
      for (comparison in adjusted) {
        if (!is.null(comparison$error)) next
        comparisons[[effect]][[comparison$source_slot]]$table <- comparison$formatted
      }
      result$comparisons <- c(result$comparisons, comparisons)
    }
    result
  }

  method_label <- function(method, followup = FALSE) {
    if (followup && method == "art") return("ART-C")
    if (followup && method == "ats") return("MCTP")
    if (followup && method == "welch") return("Games–Howell")
    switch(method, anova = "ANOVA", ancova = "ANCOVA", art = "ART", permuco = "permuco",
           permuco_ancova = "permuco ANCOVA", rank = "ランク",
           brunner_munzel = "Brunner–Munzel", ats = "ATS", t = "t検定", welch = "Welch", rank_ancova = "順位ANCOVA", method)
  }

  rank_test_label <- function(test_id) {
    switch(test_id, wilcoxon_rank_sum = "Wilcoxon順位和検定（Mann–Whitney U検定）",
      wilcoxon_signed_rank = "Wilcoxon符号付順位検定", kruskal_wallis = "Kruskal–Wallis検定",
      friedman = "Friedman検定", dunn = "Dunn検定",
      brunner_munzel = "Brunner–Munzel検定", ats = "ATS（ANOVA型統計量・F近似）", test_id)
  }

  type_label <- function(type) {
    switch(as.character(type), "1" = "I", "2" = "II", "3" = "III", as.character(type))
  }

  format_table <- function(table, digits, p_digits) {
    if (!nrow(table)) return(invisible(NULL))
    table <- as.data.frame(table, stringsAsFactors = FALSE)
    right <- intersect(c("pes", "ges", "p.value", "significance"), names(table))
    table <- table[, c(setdiff(names(table), right), right), drop = FALSE]
    for (name in names(table)) {
      values <- table[[name]]
      if (is.numeric(values)) {
        valid <- is.finite(values)
        rendered <- rep("—", length(values))
        if (grepl("(^p[._]|^p$|p.value|Mauchly.p|P\\(>|Pr\\()", name)) {
          threshold <- 10^(-p_digits)
          rendered[valid] <- formatC(values[valid], format = "f", digits = p_digits)
          rendered[valid & values < threshold] <- paste0("<", formatC(threshold, format = "f", digits = p_digits))
        } else {
          integer_df <- grepl("^df[12]?$", name) && all(values[valid] == floor(values[valid]))
          rendered[valid] <- formatC(values[valid], format = "f", digits = if (integer_df) 0 else digits)
        }
        if (grepl("(^p[._]|^p$|p.value|Mauchly.p|P\\(>|Pr\\()", name) || name %in% c("pes", "ges")) {
          rendered <- sub("^(<?)0\\.", "\\1.", rendered)
        }
        table[[name]] <- rendered
      } else {
        table[[name]] <- as.character(values)
        table[[name]][is.na(values)] <- "—"
        if (name == "correction") {
          table[[name]][table[[name]] == "none"] <- "未補正"
          table[[name]][table[[name]] == "unavailable"] <- "算出不可"
        }
      }
    }
    labels <- c(Effect = "効果", method = "方法", comparison = "比較", estimate = "推定差",
                significance = "", correction = "補正", p.value = "p", pes = "偏η²",
                ges = "一般化η²", Mauchly.p = "p", Condition = "条件",
                t.ratio = "t", z.ratio = "z", q.ratio = "q", CI = "信頼区間", mean_difference = "平均差",
                coefficient = "係数", Covariate = "共変量",
                statistic_name = "統計量", statistic = "値", relative_effect = "相対効果",
                relative_difference = "相対効果差")
    names(table) <- vapply(names(table), function(name) {
      if (name %in% names(labels)) labels[[name]] else name
    }, character(1))
    table
  }

  print_table <- function(table, digits, p_digits) {
    if (!nrow(table)) return(invisible(NULL))
    table <- format_table(table, digits, p_digits)
    widths <- vapply(seq_along(table), function(j) {
      max(nchar(c(names(table)[j], table[[j]]), type = "width"))
    }, numeric(1))
    old <- options(width = min(10000L, max(getOption("width"), sum(widths) + 3L * ncol(table) + 4L)))
    on.exit(options(old), add = TRUE)
    print(table, row.names = FALSE, right = TRUE, quote = FALSE)
    invisible(NULL)
  }

  target_heading <- function(target) {
    cat("\n========================================\n", target,
        "\n========================================\n", sep = "")
  }

  heading <- function(target, section, methods, followup = FALSE) {
    labels <- vapply(methods, method_label, character(1), followup = followup)
    title <- paste0("──【", section, "（", paste(labels, collapse = " | "), "）】")
    detail <- NULL
    if (nchar(title, type = "width") > 38L) {
      title <- paste0("──【", section, "】")
      detail <- paste0("  （", paste(labels, collapse = " | "), "）")
    }
    padding <- max(2L, 40L - nchar(title, type = "width"))
    cat("\n", title, strrep("─", padding), "\n", sep = "")
    if (!is.null(detail)) cat(detail, "\n", sep = "")
  }

  translate_notice <- function(text) {
    # 原文は結果オブジェクトに保持し、表示時だけ既知のメッセージを翻訳する。
    text <- sub("^NOTE: *", "", text)
    translations <- c(
      "HF eps > 1 treated as 1" = "HFのεが1を超えたため、1として扱いました。",
      "Results may be misleading due to involvement in interactions" =
        "比較対象の要因が交互作用に含まれるため、比較結果の解釈には注意が必要です。"
    )
    for (original in names(translations)) {
      text <- gsub(original, translations[[original]], text, fixed = TRUE)
    }
    text <- sub("^Degrees-of-freedom method: *", "自由度の算出法：", text)
    text <- sub("^Results are averaged over the levels of: *", "平均化した要因：", text)
    text <- sub("^Confidence level used: *", "信頼水準：", text)
    text <- sub("^Conf-level adjustment: *([^ ]+) method for ([0-9]+) estimates? *$",
      "信頼区間の補正法：\\1（推定値\\2個）", text)
    text <- sub("^P value adjustment: *([^ ]+) method for ([0-9]+) tests? *$",
      "p値の補正法：\\1（検定\\2個）", text)
    text
  }

  print_notes <- function(notes = character(), messages = character(), warnings = character(), label = NULL) {
    clean <- function(x) unique(translate_notice(sub("[[:space:]]+$", "", as.character(x))))
    scope <- if (is.null(label)) "" else paste0(" [", label, "]")
    lines <- c(if (length(notes)) paste0("注記", scope, "：", clean(notes)),
               if (length(messages)) paste0("通知", scope, "：", clean(messages)),
               if (length(warnings)) paste0("警告", scope, "：", clean(warnings)))
    if (length(lines)) cat("\n", paste(lines, collapse = "\n"), "\n", sep = "")
    invisible(NULL)
  }

  print_method_tables <- function(tables, key, digits, p_digits, row_keys = NULL, captions = NULL,
                                  before_table = NULL, after_table = NULL, followup = FALSE) {
    if (is.null(row_keys)) row_keys <- unique(unlist(lapply(tables, `[[`, key), use.names = FALSE))
    if (!length(row_keys)) return(invisible(NULL))
    for (method in names(tables)) {
      table <- tables[[method]]
      if (anyDuplicated(table[[key]])) stop("比較表の行キーが重複しています。", call. = FALSE)
      current_keys <- if (followup) row_keys[row_keys %in% table[[key]]] else row_keys
      aligned <- table[match(current_keys, table[[key]]), setdiff(names(table), key), drop = FALSE]
      displayed <- data.frame(current_keys, aligned, check.names = FALSE)
      names(displayed)[1L] <- key
      if (method %in% names(captions)) cat("\n", captions[[method]], "\n", sep = "")
      else if (length(tables) > 1L) cat("\n", method_label(method, followup), "\n", sep = "")
      if (is.function(before_table)) before_table(method, displayed)
      print_table(displayed, digits, p_digits)
      if (followup && !nrow(displayed)) cat("表示条件に該当する比較はありません。\n")
      if (is.function(after_table)) after_table(method, displayed)
    }
    invisible(NULL)
  }

  sphericity_caption <- function(result) {
    if (!length(result$within)) return("球面性補正：対象なし")
    labels <- c(none = "未補正", GG = "GG", HF = "HF", unavailable = "算出不可")
    actual <- unique(as.character(result$test_table$correction))
    actual <- actual[!is.na(actual)]
    used <- unname(labels[actual])
    used[is.na(used)] <- actual[is.na(used)]
    if (!length(used)) used <- "算出不可"
    requested <- if (identical(result$sphericity, "auto")) "自動選択" else
      if (result$sphericity %in% names(labels)) paste0(labels[[result$sphericity]], "指定") else as.character(result$sphericity)
    paste0("球面性補正：", requested, " → ", paste(used, collapse = " / "), "（効果別は補正列）")
  }

  display_sphericity <- function(values, target, digits, p_digits) {
    rows <- list()
    for (method in names(values)) {
      sphere <- values[[method]]$sphericity_table
      if (is.null(sphere)) next
      sphere <- sphere[nzchar(sphere$Within), , drop = FALSE]
      if (!nrow(sphere)) next
      sphere$method <- method_label(method)
      decisions <- values[[method]]$sphericity_decisions
      sphere$correction <- if (is.null(decisions)) "none" else decisions$correction[match(sphere$Effect, decisions$Effect)]
      rows[[method]] <- sphere[, c("Effect", "method", "W", "Mauchly.p", "GG", "HF", "correction")]
    }
    available <- union(names(rows), names(values)[lengths(lapply(values, `[[`, "sphericity_errors")) > 0L])
    if (!length(available)) return(invisible(NULL))
    heading(target, "球面性", available)
    rows <- lapply(rows, function(table) table[, setdiff(names(table), "method"), drop = FALSE])
    supplement <- function(method, table = NULL) {
      errors <- values[[method]]$sphericity_errors
      warnings <- vapply(names(errors), function(effect) paste0("球面性推定失敗 [", effect, "] ", errors[[effect]]), character(1))
      print_notes(if (method == "art" && method %in% names(rows))
        "効果別順位データによる評価。GG/HF補正は拡張処理。" else character(),
        warnings = warnings, label = method_label(method))
    }
    print_method_tables(rows, "Effect", digits, p_digits, after_table = supplement)
    for (method in setdiff(available, names(rows))) supplement(method)
    invisible(NULL)
  }

  display_followup_intro <- function(values, target, digits, p_digits, followup_alpha) {
    values <- values[vapply(seq_along(values), function(i) names(values)[i] %in% c("anova", "ancova", "art") ||
      names(values)[i] == "ats" ||
      (names(values)[i] == "welch" && values[[i]]$n_levels > 2L) ||
      (names(values)[i] == "rank_ancova" && values[[i]]$n_levels > 2L) ||
      (names(values)[i] == "rank" && values[[i]]$n_levels > 2L), logical(1))]
    if (!length(values)) return(invisible(NULL))
    heading(target, "下位検定", names(values), followup = TRUE)
    if (all(vapply(values, function(result) identical(result$posthoc, FALSE), logical(1)))) {
      cat("下位検定は実行しない設定です。\n")
      return(invisible(NULL))
    }
    effects_text <- function(x) if (length(x)) paste(x, collapse = ", ") else "なし"
    for (method in names(values)) {
      result <- values[[method]]
      table <- data.frame(method = method_label(method, TRUE), 有意な効果 = effects_text(result$significant_effects),
        有意傾向 = effects_text(result$trend_effects),
        推定差の尺度 = if (method == "rank_ancova") "各ペア内の順位尺度（推定差・CIなし）" else if (method == "ats") "群を等しく重み付けした相対効果の差" else if (method == "art") "順位尺度（ART-C）" else if (method == "rank")
          "対象外（順位に基づく検定）" else "元尺度", check.names = FALSE)
      if (length(values) > 1L) cat("\n")
      print_table(table, digits, p_digits)
      notes <- if (method == "ancova")
        "推定差は共変量の参加者平均で評価。共変量を含む効果は水準間比較の対象外です。" else character()
      if (method == "rank") notes <- c(notes,
        paste0("下位検定：", rank_test_label(result$followup_test_id), "。全ペアをまとめて補正します。"),
        "この下位検定では信頼区間を算出しません。",
        if (result$paired) "各ペアの差に順位を付ける検定で、差の分布の対称性を前提とします。" else
          "全群を通じた順位を用いるDunn検定です。")
      if (method == "ats") notes <- c(notes,
        "下位検定：全ペアのMCTP（Fisher近似）。p値と同時信頼区間は補正済みです。",
        "相対効果差は第1水準−第2水準。元尺度の平均差や2群だけの優越確率とは異なります。")
      if (method == "welch") notes <- c(notes,
        "下位検定：全ペアのGames–Howell。p値は全ペアを考慮済みで、追加補正は行いません。",
        "推定差は第1水準−第2水準の平均差。信頼区間はGames–Howellの近似同時信頼区間です。")
      if (method == "rank_ancova") notes <- c(notes,
        "下位検定：各2群内で順位化・共変量調整をやり直すQuade法。全ペアをまとめて補正します。",
        "元尺度の推定差・信頼区間は算出しません。")
      if (length(result$forced_simple_effects)) notes <- c(notes,
        paste0("追加の単純比較 ", effects_text(result$forced_simple_effects)))
      print_notes(notes, label = method_label(method, TRUE))
    }
    alpha <- values[[1L]]$alpha
    trend_alpha <- values[[1L]]$trend_alpha
    if (is.null(trend_alpha)) trend_alpha <- followup_alpha
    threshold_text <- function(x) sub("^0\\.", ".", format(x, scientific = FALSE, trim = TRUE))
    cat("\n有意 p < ", threshold_text(alpha), " / 有意傾向 ", threshold_text(alpha),
        " ≤ p < ", threshold_text(trend_alpha), "\n", sep = "")
    adjustments <- vapply(values, `[[`, character(1), "adjust")
    descriptions <- if (length(unique(adjustments)) == 1L) unique(adjustments) else
      vapply(names(values), function(method) {
        if (method %in% c("ats", "welch")) adjustments[[method]] else
          paste(method_label(method, TRUE), adjustments[[method]])
      }, character(1))
    config <- values[[1L]]
    scope <- if (is.null(config$followup_scope)) "followup" else config$followup_scope
    shown <- if (is.null(config$display_scope)) "significant" else config$display_scope
    cat("計算範囲 ", followup_scope_label(scope),
      if (scope == "followup") paste0("（p < ", threshold_text(followup_alpha), "）") else "",
      " / 表示範囲 ", display_scope_label(shown), "\n", sep = "")
    cat("多重比較補正（指定／固有） ", paste(descriptions, collapse = " / "), "\n", sep = "")
    p_adjustments <- lapply(values, function(result) unique(unlist(lapply(result$contrasts, `[[`, "p_adjustment"), use.names = FALSE)))
    p_adjustments <- p_adjustments[lengths(p_adjustments) > 0L]
    if (length(p_adjustments)) {
      descriptions <- vapply(names(p_adjustments), function(method) paste0(method_label(method, TRUE), " ",
        paste(p_adjustments[[method]], collapse = ", ")), character(1))
      cat("p値の補正（実際） ", paste(descriptions, collapse = " / "), "\n", sep = "")
    }
    ci_adjustments <- lapply(values, function(result) setdiff(unique(unlist(lapply(result$contrasts, `[[`, "ci_adjustment"), use.names = FALSE)), "not_computed"))
    ci_adjustments <- ci_adjustments[lengths(ci_adjustments) > 0L]
    if (length(ci_adjustments)) {
      descriptions <- vapply(names(ci_adjustments), function(method) paste0(method_label(method, TRUE), " ",
        paste(ci_adjustments[[method]], collapse = ", ")), character(1))
      cat("信頼区間の補正 ", paste(descriptions, collapse = " / "), "\n", sep = "")
    }
    confidence <- lapply(values, function(result) {
      messages <- as.character(unlist(lapply(result$contrasts, `[[`, "adjustment_messages"), use.names = FALSE))
      levels <- suppressWarnings(as.numeric(sub("^Confidence level used: *", "",
        messages[startsWith(messages, "Confidence level used:")])))
      if (result$analysis %in% c("ats", "welch")) levels <- c(levels, result$confidence_level)
      unique(levels[is.finite(levels) & levels > 0 & levels < 1])
    })
    confidence <- confidence[lengths(confidence) > 0L]
    if (length(confidence)) {
      percent <- function(x) paste0(format(100 * x, scientific = FALSE, trim = TRUE), "%")
      levels <- unique(unlist(confidence, use.names = FALSE))
      if (length(levels) == 1L && length(confidence) == length(values)) {
        cat("信頼水準 ", percent(levels), "\n", sep = "")
      } else {
        descriptions <- vapply(names(confidence), function(method) {
          paste0(method_label(method, TRUE), " ", paste(percent(confidence[[method]]), collapse = ", "))
        }, character(1))
        cat("信頼水準 ", paste(descriptions, collapse = " / "), "\n", sep = "")
      }
    }
    if (!any(lengths(lapply(values, `[[`, "follow_up_effects")))) cat("下位検定の対象となる効果はありません。\n")
    invisible(NULL)
  }

  display_followup <- function(values, target, digits, p_digits, followup_alpha) {
    original_levels <- function(result, method, factor_name, labels) {
      levels <- levels(factor(result$analysis_data[[factor_name]]))
      native <- levels
      if (method %in% c("anova", "ancova")) {
        model_levels <- levels(result$model$data$long[[factor_name]])
        if (length(model_levels) == length(levels)) native <- model_levels
      }
      indices <- match(labels, native)
      missing <- is.na(indices)
      indices[missing] <- match(labels[missing], paste0(factor_name, native))
      converted <- labels
      converted[!is.na(indices)] <- levels[indices[!is.na(indices)]]
      converted
    }
    rows <- list()
    effect_names <- character()
    failures <- setNames(lapply(values, function(x) character()), names(values))
    for (method in names(values)) {
      result <- values[[method]]
      for (name in names(result$contrasts)) {
        comparison <- result$contrasts[[name]]
        effect <- if (isTRUE(comparison$kind %in% c("rank", "mctp", "games_howell", "rank_ancova"))) comparison$source_effect else sub("/.*", "", name)
        effect_names <- union(effect_names, effect)
        if (!is.null(comparison$error)) {
          failures[[method]] <- c(failures[[method]], paste0("下位検定失敗 [", name, "] ", comparison$error))
          next
        }
        table <- comparison$formatted
        factors <- strsplit(effect, ":", fixed = TRUE)[[1L]]
        compared <- if ("compared_factor" %in% names(table)) as.character(table$compared_factor)
          else rep(if (is.null(comparison$compared_factor)) factors[1L] else comparison$compared_factor, nrow(table))
        fixed_columns <- intersect(c(result$within, result$between), names(table))
        conditions <- if ("fixed_conditions" %in% names(table)) as.character(table$fixed_conditions)
          else if (length(fixed_columns)) apply(table[, fixed_columns, drop = FALSE], 1, function(row) {
            paste(paste0(fixed_columns, "=", row), collapse = ", ")
          }) else rep("", nrow(table))
        if (length(fixed_columns)) {
          condition_table <- table[, fixed_columns, drop = FALSE]
          for (factor_name in fixed_columns) condition_table[[factor_name]] <-
            original_levels(result, method, factor_name, as.character(condition_table[[factor_name]]))
          conditions <- apply(condition_table, 1, function(row) paste(paste0(fixed_columns, "=", row), collapse = ", "))
        } else if (all(c("fixed_factor", "fixed_level") %in% names(table))) {
          for (j in seq_len(nrow(table))) if (table$fixed_factor[j] %in% c(result$within, result$between)) {
            conditions[j] <- paste0(table$fixed_factor[j], "=", original_levels(result, method,
              table$fixed_factor[j], as.character(table$fixed_level[j])))
          }
        }
        labels <- table$comparison
        for (column in fixed_columns) {
          for (j in seq_len(nrow(table))) {
            fixed <- paste0(column, "=", table[[column]][j])
            labels[j] <- gsub(paste0(" | ", fixed), "", labels[j], fixed = TRUE)
            labels[j] <- gsub(paste0(fixed, " | "), "", labels[j], fixed = TRUE)
          }
        }
        # emmeansの比較ラベルから、同じ向きの水準対を方法間で対応させる。
        raw_labels <- if (isTRUE(comparison$kind %in% c("rank", "mctp", "games_howell", "rank_ancova"))) character() else {
          grid <- methods::slot(comparison$comparison_grid, "grid")
          as.character(grid$contrast)
        }
        pair_keys <- if (isTRUE(comparison$kind %in% c("rank", "mctp", "games_howell", "rank_ancova"))) table$pair_key else
          paste(method, name, seq_len(nrow(table)), sep = "/")
        if (length(raw_labels) == nrow(table)) for (j in seq_len(nrow(table))) {
          pair <- trimws(strsplit(raw_labels[j], " - ", fixed = TRUE)[[1L]])
          if (length(pair) != 2L) next
          if (method == "art" && length(factors) > 1L) {
            sides <- strsplit(pair, ",", fixed = TRUE)
            position <- match(compared[j], factors)
            if (is.na(position) || any(lengths(sides) != length(factors))) next
            pair <- vapply(sides, function(side) trimws(side[position]), character(1))
          }
          pair <- original_levels(result, method, compared[j], pair)
          level_order <- levels(factor(result$analysis_data[[compared[j]]]))
          indices <- match(pair, level_order)
          if (anyNA(indices)) next
          pair_keys[j] <- paste(indices, collapse = ":")
          labels[j] <- paste(pair, collapse = " - ")
        }
        column <- function(name) if (name %in% names(table)) table[[name]] else rep(NA_real_, nrow(table))
        ci <- if (all(c("lower.CL", "upper.CL") %in% names(table))) {
          ifelse(is.finite(table$lower.CL) & is.finite(table$upper.CL),
            paste0("[", formatC(table$lower.CL, format = "f", digits = digits), ", ", formatC(table$upper.CL, format = "f", digits = digits), "]"), NA_character_)
        } else rep(NA_character_, nrow(table))
        rows[[length(rows) + 1L]] <- data.frame(
          group = ifelse(nzchar(conditions), paste(compared, conditions, sep = " / "), compared),
          method = method, source = name, visible = visible_followup_rows(result, effect, table$p.value), pair_key = pair_keys, comparison = labels, estimate = column("estimate"),
          SE = column("SE"), df = column("df"), df1 = column("df1"), df2 = column("df2"), t.ratio = column("t.ratio"), z.ratio = column("z.ratio"),
          q.ratio = column("q.ratio"),
          CI = ci, statistic_name = if ("statistic_name" %in% names(table)) table$statistic_name else NA_character_,
          statistic = column("statistic"), p.value = table$p.value, significance = table$significance,
          family_size = table$family_size, stringsAsFactors = FALSE)
      }
    }
    cat("\n【", paste(effect_names, collapse = ", "), "】\n", sep = "")
    supplements <- function(method, sources = names(values[[method]]$contrasts), include_runtime = FALSE) {
      result <- values[[method]]
      comparisons <- result$contrasts[sources]
      notes <- as.character(unique(unlist(lapply(comparisons, `[[`, "adjustment_messages"), use.names = FALSE)))
      # 信頼水準・補正法は冒頭へ集約。平均化要因と比較数は表示しない。
      shared <- startsWith(notes, "Confidence level used:") |
        startsWith(notes, "Results are averaged over the levels of:")
      adjustment <- grepl("^(Conf-level|P value) adjustment:", notes)
      actual_adjust <- sub(" method.*$", "", sub("^[^:]+: *", "", notes))
      redundant <- adjustment & actual_adjust %in% result$adjust
      interval_header <- startsWith(notes, "Conf-level adjustment:")
      notes <- notes[!shared & !redundant & !interval_header]
      for (comparison in comparisons) {
        if (!is.null(comparison$error)) next
        current_p <- comparison$p_adjustment
        current_ci <- setdiff(comparison$ci_adjustment, "not_computed")
        if (length(current_p) && !setequal(current_p, result$header_p_adjustment))
          notes <- c(notes, paste0("この比較のp値の補正：", paste(current_p, collapse = ", ")))
        if (length(current_ci) && !setequal(current_ci, result$header_ci_adjustment))
          notes <- c(notes, paste0("この比較の信頼区間の補正：", paste(current_ci, collapse = ", ")))
      }
      warnings <- unlist(lapply(comparisons, function(comparison) {
        if (is.null(comparison$warnings)) return(character())
        if (is.character(comparison$warnings)) return(comparison$warnings)
        vapply(comparison$warnings, function(w) w$message, character(1))
      }), use.names = FALSE)
      warnings <- c(warnings, failures[[method]])
      if (include_runtime) warnings <- c(warnings, result$warnings)
      print_notes(notes, messages = if (include_runtime) result$messages else character(),
        warnings = warnings, label = method_label(method, TRUE))
    }
    table <- if (length(rows)) do.call(rbind, rows) else NULL
    for (method in setdiff(names(values), unique(table$method))) {
      cat("\n", method_label(method, TRUE), "\n", sep = "")
      supplements(method, include_runtime = TRUE)
    }
    if (is.null(table)) return(invisible(NULL))
    runtime_printed <- character()
    for (group in unique(table$group)) {
      complete <- table[table$group == group, , drop = FALSE]
      if (!group %in% effect_names) cat("\n", group, "\n", sep = "")
      # 方法ごとに表示条件を判定。非表示の比較は補正・保存結果から削除しない。
      selected_keys <- unique(complete$pair_key[complete$visible])
      key_labels <- complete$comparison[match(selected_keys, complete$pair_key)]
      if (anyDuplicated(key_labels)) key_labels <- paste0(key_labels, " [", selected_keys, "]")
      tables <- list()
      for (method in names(values)[names(values) %in% unique(complete$method)]) {
        shown <- complete[complete$method == method & complete$visible & complete$pair_key %in% selected_keys, , drop = FALSE]
        shown$comparison <- key_labels[match(shown$pair_key, selected_keys)]
        shown <- shown[, if (method == "rank_ancova") c("comparison", "df1", "df2", "statistic", "p.value", "significance") else
          if (method == "rank") c("comparison", "statistic_name", "statistic", "p.value", "significance") else
          c("comparison", "estimate", "SE", "df", "t.ratio", "z.ratio", "q.ratio", "CI", "p.value", "significance")]
        if (method == "rank_ancova") names(shown)[names(shown) == "statistic"] <- "F"
        if (method == "ats") names(shown)[names(shown) == "estimate"] <- "relative_difference"
        for (name in intersect(c("t.ratio", "z.ratio", "q.ratio", "CI"), names(shown))) if (all(is.na(shown[[name]]))) shown[[name]] <- NULL
        tables[[method]] <- shown
      }
      before <- function(method, displayed) {
        own <- complete[complete$method == method & complete$visible, , drop = FALSE]
        alpha <- values[[method]]$alpha
        trend_alpha <- values[[method]]$trend_alpha
        if (is.null(trend_alpha)) trend_alpha <- followup_alpha
        pairs_text <- function(selected) {
          pairs <- unique(own$comparison[is.finite(own$p.value) & selected])
          if (length(pairs)) paste(pairs, collapse = ", ") else "なし"
        }
        cat("有意なペア：", pairs_text(own$p.value < alpha), "\n",
            "有意傾向のペア：", pairs_text(own$p.value >= alpha & own$p.value < trend_alpha), "\n\n", sep = "")
      }
      after <- function(method, displayed) {
        sources <- unique(complete$source[complete$method == method])
        supplements(method, sources, include_runtime = !method %in% runtime_printed)
        runtime_printed <<- union(runtime_printed, method)
        failures[[method]] <<- character()
      }
      if (length(selected_keys)) {
        print_method_tables(tables, "comparison", digits, p_digits, row_keys = key_labels,
          before_table = before, after_table = after, followup = TRUE)
      } else for (method in names(tables)) {
        if (length(tables) > 1L) cat("\n", method_label(method, TRUE), "\n", sep = "")
        before(method, NULL)
        cat("表示基準を満たす比較はありません。\n")
        after(method, NULL)
      }
    }
    invisible(NULL)
  }

  mean_t_methods <- function(values) {
    names(values)[vapply(values, function(result) isTRUE(result$test_id %in% c("student_t", "paired_t", "welch_t")), logical(1))]
  }

  anova_methods <- function(values) {
    names(values)[vapply(seq_along(values), function(i) names(values)[i] %in% c("anova", "art", "permuco") ||
      (names(values)[i] == "welch" && values[[i]]$n_levels > 2L), logical(1))]
  }

  display_t_tests <- function(values, target, digits, p_digits) {
    values <- values[mean_t_methods(values)]
    if (!length(values)) return(invisible(NULL))
    heading(target, "平均値の比較", names(values))
    tables <- lapply(values, function(result) {
      table <- result$test_table
      data.frame(Effect = table$Effect, mean_difference = table$estimate, SE = table$SE,
        df = table$df, t.ratio = table$statistic,
        CI = ifelse(is.finite(table$lower.CL) & is.finite(table$upper.CL),
          paste0("[", formatC(table$lower.CL, digits = digits, format = "f"), ", ",
            formatC(table$upper.CL, digits = digits, format = "f"), "]"), NA_character_),
        p.value = table$p.value, significance = table$significance)
    })
    captions <- lapply(values, function(result) switch(result$test_id,
      student_t = "Studentのt検定（対応なし・等分散を仮定）",
      paired_t = "対応のあるt検定", welch_t = "Welchのt検定（対応なし）"))
    supplement <- function(method, table) {
      result <- values[[method]]
      groups <- paste(paste0(names(result$group_sizes), " n=", as.integer(result$group_sizes)), collapse = " / ")
      notes <- c(paste0("平均差：", result$levels[1L], " − ", result$levels[2L], " / ", groups),
        paste0("両側検定 / 信頼水準 ", format(100 * result$confidence_level, trim = TRUE), "%"),
        if (result$paired) "同じIDで対応させた差を検定しています。差の正規性を前提とします。")
      print_notes(notes, messages = result$messages, warnings = result$warnings, label = method_label(method))
    }
    print_method_tables(tables, "Effect", digits, p_digits, captions = captions, after_table = supplement)
    invisible(NULL)
  }

  display_results <- function(collection, digits, p_digits, followup_alpha, stage = "all") {
    stage <- match.arg(stage, c("all", "sphericity", "anova", "followup_intro", "followup", "diagnostics", "ancova", "rank", "t"))
    covariate_stage <- stage == "ancova"
    if (covariate_stage) stage <- "anova"
    targets <- unique(unlist(lapply(collection, names), use.names = FALSE))
    for (target in targets) {
      values <- lapply(collection, function(results) results[[target]])
      values <- Filter(Negate(is.null), values)
      if (stage == "all") {
        target_heading(target)
        subset <- lapply(values, function(result) setNames(list(result), target))
        for (phase in c("sphericity", "diagnostics", "t", "anova", "ancova", "rank", "followup_intro", "followup")) {
          selected <- if (phase == "anova") anova_methods(values) else
            if (phase == "t") mean_t_methods(values) else
            if (phase == "ancova") intersect(names(subset), c("ancova", "permuco_ancova", "rank_ancova")) else
            if (phase == "rank") intersect(names(subset), c("rank", "brunner_munzel", "ats")) else names(subset)
          if (length(selected)) display_results(subset[selected], digits, p_digits, followup_alpha, phase)
        }
        next
      }
      if (stage %in% c("sphericity", "followup_intro", "followup")) {
        switch(stage,
          sphericity = display_sphericity(values, target, digits, p_digits),
          followup_intro = display_followup_intro(values, target, digits, p_digits, followup_alpha),
          followup = display_followup(values, target, digits, p_digits, followup_alpha))
        next
      }
      if (stage == "t") {
        display_t_tests(values, target, digits, p_digits)
        next
      }
      if (stage == "rank") {
        direct_methods <- intersect(names(values), c("rank", "brunner_munzel", "ats"))
        if (!length(direct_methods)) next
        heading(target, "順位検定", direct_methods)
        for (method in direct_methods) {
          result <- values[[method]]
          cat("\n", rank_test_label(result$test_id), "\n", sep = "")
          table <- result$test_table
          shown <- table[, intersect(c("Effect", "statistic_name", "statistic", "df", "df1", "df2",
            "relative_effect", "p.value", "significance"), names(table)), drop = FALSE]
          if ("df" %in% names(shown) && all(is.na(shown$df))) shown$df <- NULL
          if (method == "brunner_munzel") {
            shown$CI <- paste0("[", formatC(table$lower.CL, format = "f", digits = digits), ", ",
              formatC(table$upper.CL, format = "f", digits = digits), "]")
            shown$CI[!is.finite(table$lower.CL) | !is.finite(table$upper.CL)] <- NA_character_
          }
          print_table(shown, digits, p_digits)
          notes <- c(paste0("参加者数 ", result$n_participants,
            " / ", if (result$paired) "対応あり" else "対応なし", " / ", result$n_levels, "水準"),
            if (result$test_id == "wilcoxon_signed_rank") "参加者ごとの差の分布の対称性を前提とします。",
            if (result$test_id %in% c("wilcoxon_rank_sum", "kruskal_wallis"))
              "群間で分布の形が異なる場合、中央値だけの差とは解釈できません。",
            if (result$n_levels == 2L) "2水準のため、同じ比較の下位検定は行いません。")
          if (method == "brunner_munzel") notes <- c(notes,
            paste0("相対効果：P(X < Y) + .5 P(X = Y)。X=", result$levels[1L], " / Y=", result$levels[2L], "。帰無値 .5。"),
            paste0("両側t近似 / 信頼水準 ", format(100 * result$confidence_level, trim = TRUE),
              "%（非制約のt信頼区間）。"))
          if (method == "ats") notes <- c(notes,
            "H0p：群を等しく重み付けした相対効果が等しいかを検定します（unweighted pseudo-ranks）。")
          print_notes(notes, messages = result$messages, warnings = result$warnings, label = method_label(method))
        }
        next
      }
      if (stage %in% c("all", "anova")) {
        tables <- list()
        captions <- list()
        section <- if (covariate_stage) "共分散分析" else "分散分析"
        heading(target, section, names(values))
        for (method in names(values)) {
          result <- values[[method]]
          table <- result$test_table
          if (is.null(table)) next
          if (!"Effect" %in% names(table)) table$Effect <- rownames(table)
          if (method %in% c("permuco", "permuco_ancova") && "Residuals" %in% table$Effect) {
            table$df2 <- table$df[match("Residuals", table$Effect)]
            table <- table[table$Effect != "Residuals", , drop = FALSE]
          }
          column <- function(candidates) {
            found <- intersect(candidates, names(table))
            if (length(found)) table[[found[1L]]] else rep(NA_real_, nrow(table))
          }
          tables[[method]] <- data.frame(Effect = table$Effect,
            df1 = column(c("df1", "dfn", "df")), df2 = column(c("df2", "dfd")), F = column(c("F", "F value")),
            pes = column("pes"), ges = column("ges"), p.value = column(c("p.value", "resampled P(>F)")),
            significance = significance(column(c("p.value", "resampled P(>F)"))))
          for (name in c("pes", "ges")) if (all(is.na(tables[[method]][[name]]))) tables[[method]][[name]] <- NULL
          if (method %in% c("anova", "ancova", "art") && length(result$within)) {
            tables[[method]]$correction <- table$correction
          }
          caption <- if (method == "rank_ancova") "Quadeの順位ANCOVA（対応なし・1要因）" else if (method == "welch") "Welchの一元配置分散分析（対応なし・不等分散を許容）" else
            paste0(method_label(method), "  Type ", type_label(result$anova_type))
          if (!is.null(result$permutation_method)) {
            caption <- paste0(caption, " / 置換法 ", result$permutation_method, " / ", result$np, "回")
          } else if (!method %in% c("welch", "rank_ancova")) {
            caption <- paste0(caption, " / ", sphericity_caption(result))
            if (isTRUE(result$legacy_sphericity_override)) caption <- paste0(caption, "（旧sphericity指定）")
          }
          if (method == "rank_ancova") caption <- paste0(caption,
            "\n共変量（順位尺度で調整） ", paste(result$covariate, collapse = ", "))
          if (isTRUE(result$covariate_adjusted) && method != "rank_ancova") {
            centers <- paste(paste0(names(result$covariate_centers), " = ",
              formatC(result$covariate_centers, digits = digits, format = "f")), collapse = ", ")
            within_covariate <- if (method == "ancova") result$within_covariate_interactions else
              length(result$within) > 0L && isTRUE(result$covariate_interactions)
            caption <- paste0(caption, "\n共変量の中心値（参加者平均） ", centers,
              "\n群間要因×共変量 ", if (!length(result$between)) "対象なし" else if (result$covariate_interactions) "あり" else "なし（共通傾斜）",
              " / 被験者内要因×共変量 ", if (!length(result$within)) "対象なし" else if (within_covariate) "あり" else "なし")
          }
          captions[[method]] <- caption
        }
        supplement <- function(method, table) {
          result <- values[[method]]
          if (method == "rank_ancova") {
            cat("\n【共変量の補助回帰（順位尺度）】\n")
            cov_table <- result$covariate_table
            names(cov_table)[match(c("Effect", "estimate"), names(cov_table))] <- c("Covariate", "coefficient")
            print_table(cov_table, digits, p_digits)
          }
          notes <- if (method == "art" && "pes" %in% names(table))
            "偏η²は順位変換後の尺度に基づきます。" else character()
          if (method == "rank_ancova") notes <- c(notes,
            paste0("参加者数 ", result$n_participants, " / 対応なし / ", result$n_levels, "群"),
            "補助回帰は群を含めない順位応答と順位共変量の回帰です。共変量のpはANCOVAの同時検定とは異なります。",
            "Quade法は、群間で共変量の周辺分布が共通であることを前提とします。")
          print_notes(notes, messages = result$messages, warnings = result$warnings, label = method_label(method))
        }
        print_method_tables(tables, "Effect", digits, p_digits, captions = captions, after_table = supplement)
      }
      for (method in names(collection)) {
        result <- collection[[method]][[target]]
        if (is.null(result)) next
        if (stage %in% c("all", "diagnostics") && !is.null(result$diagnostics)) {
          diagnostic <- result$diagnostics
          heading(target, "整列診断", method)
          if (!is.null(diagnostic$aligned.col.sums) && !is.null(diagnostic$aligned.anova)) {
            f_column <- intersect(c("F", "F value", "F.value"), names(diagnostic$aligned.anova))
            cat("整列値の列和（最大絶対値） ", format(max(abs(diagnostic$aligned.col.sums)), digits = digits), "\n", sep = "")
            if (length(f_column)) {
              f_values <- diagnostic$aligned.anova[[f_column[1]]]
              f_values <- f_values[is.finite(f_values)]
              cat("対象外効果のF（最大値） ", if (length(f_values)) format(max(f_values), digits = digits) else "対象なし", "\n", sep = "")
            }
          } else print(diagnostic)
        }
        if (stage %in% c("all", "diagnostics") && !is.null(result$levene)) {
          heading(target, "Levene検定", method)
          table <- result$levene
          print_table(table[, intersect(c("Condition", "df1", "df2", "F", "p.value", "significance"), names(table)), drop = FALSE], digits, p_digits)
        }
      }
    }
    invisible(NULL)
  }
}, envir = .test_private_environment)


# GUIは引数の受付のみを担当し、確定後は通常のanalyze()へ渡す。
evalq({
  validate_rank_settings <- function(experiment, adjust) {
    if (length(experiment$iv) != 1L ||
        !((identical(as.character(experiment$within), experiment$iv) && !length(experiment$between)) ||
          (identical(as.character(experiment$between), experiment$iv) && !length(experiment$within)))) {
      stop("ランクは1要因のwithinまたはbetween計画に対応しています。", call. = FALSE)
    }
    if (!is.character(adjust) || length(adjust) != 1L || is.na(adjust) || !adjust %in% stats::p.adjust.methods) {
      stop("ランクのadjustは次から指定してください: ", paste(stats::p.adjust.methods, collapse = ", "), call. = FALSE)
    }
    invisible(NULL)
  }

  validate_rank_ancova_settings <- function(experiment, adjust, covariate_interactions = NULL) {
    .t_rank_ancova_design(experiment)
    if (!is.character(adjust) || length(adjust) != 1L || is.na(adjust) || !adjust %in% stats::p.adjust.methods) {
      stop("順位ANCOVAのadjustは次から指定してください: ", paste(stats::p.adjust.methods, collapse = ", "), call. = FALSE)
    }
    if (!is.null(covariate_interactions) && !identical(covariate_interactions, FALSE)) {
      stop("順位ANCOVA（Quade）は共変量との交互作用を扱いません。", call. = FALSE)
    }
    invisible(NULL)
  }

  validate_common_settings <- function(alpha, followup_alpha, digits, p_digits, correction, sphericity) {
    scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
    if (!scalar(alpha) || !scalar(followup_alpha) || alpha <= 0 || alpha >= 1 || followup_alpha <= 0 || followup_alpha >= 1)
      stop("alphaとfollowup_alphaは0より大きく1未満にしてください。", call. = FALSE)
    if (!scalar(digits) || digits < 0 || digits != floor(digits) ||
        !scalar(p_digits) || p_digits < 1 || p_digits != floor(p_digits))
      stop("digitsは0以上、p_digitsは1以上の整数にしてください。", call. = FALSE)
    valid <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && x %in% c("auto", "none", "GG", "HF")
    if (!valid(correction) || (!is.null(sphericity) && !valid(sphericity)))
      stop("球面性補正はauto / GG / HF / noneから指定してください。", call. = FALSE)
    invisible(NULL)
  }

  resolve_shared_threshold <- function(trend_alpha = NULL, followup_alpha = NULL) {
    if (!is.null(trend_alpha) && !is.null(followup_alpha) && !isTRUE(all.equal(trend_alpha, followup_alpha)))
      stop("有意傾向と下位検定の基準は共用です。trend_alphaと旧followup_alphaに別の値は指定できません。trend_alphaだけを指定してください。", call. = FALSE)
    if (!is.null(trend_alpha)) trend_alpha else if (!is.null(followup_alpha)) followup_alpha else .10
  }

  validate_scope_settings <- function(alpha, trend_alpha, followup_alpha, posthoc, followup_scope, display_scope) {
    resolve_shared_threshold(trend_alpha, followup_alpha)
    if (!is.numeric(trend_alpha) || length(trend_alpha) != 1L || !is.finite(trend_alpha) || trend_alpha < alpha || trend_alpha >= 1)
      stop("有意傾向の上限trend_alphaはalpha以上、1未満にしてください。", call. = FALSE)
    if (!is.logical(posthoc) || length(posthoc) != 1L || is.na(posthoc)) stop("posthocにはTRUE/FALSEを指定してください。", call. = FALSE)
    valid <- function(x, choices) is.character(x) && length(x) == 1L && !is.na(x) && x %in% choices
    if (!valid(followup_scope, c("significant", "followup", "all"))) stop("followup_scopeはsignificant / followup / allから指定してください。", call. = FALSE)
    if (!valid(display_scope, c("significant", "effect", "all"))) stop("display_scopeはsignificant / effect / allから指定してください。", call. = FALSE)
    invisible(NULL)
  }

  aggregate_observation_units <- function(experiment, method) {
    dat <- experiment$data
    keys <- c(experiment$id, experiment$iv)
    if (!is.data.frame(dat) || !all(keys %in% names(dat)) || !anyDuplicated(dat[keys]))
      return(list(experiment = experiment, info = NULL))
    if (anyNA(dat[keys])) stop("ID・要因列に欠損があります。", call. = FALSE)
    for (name in experiment$between) {
      changed <- vapply(split(dat[[name]], dat[[experiment$id]], drop = TRUE),
        function(x) length(unique(x)) != 1L, logical(1))
      if (any(changed)) stop("被験者間要因がID内で変化しています: ", name, call. = FALSE)
    }
    target <- experiment$dv
    if (length(target) != 1L || !is.numeric(dat[[target]]) || any(!is.finite(dat[[target]])))
      stop("セル内反復の集約には欠損・無限値のない数値の従属変数が必要です。", call. = FALSE)
    if (method == "permuco_ancova") for (name in experiment$covariate) {
      if (anyNA(dat[[name]]) || any(vapply(split(dat[[name]], dat[[experiment$id]], drop = TRUE),
          function(x) length(unique(x)) != 1L, logical(1))))
        stop("共変量は欠損のない参加者内で一定の値が必要です: ", name, call. = FALSE)
    }
    group_codes <- lapply(dat[keys], function(x) as.integer(factor(x)))
    groups <- split(seq_len(nrow(dat)), interaction(group_codes, drop = TRUE, lex.order = TRUE))
    counts <- lengths(groups)
    result <- dat[vapply(groups, `[`, integer(1), 1L), , drop = FALSE]
    result[[target]] <- vapply(groups, function(i) mean(dat[[target]][i]), numeric(1))
    rownames(result) <- NULL
    experiment$data <- result
    list(experiment = experiment, info = list(aggregation = "mean", n_observations = nrow(dat),
      n_participant_cells = nrow(result), repetitions_min = min(counts), repetitions_max = max(counts)))
  }

  expand_methods <- function(method) {
    choices <- c("anova", "ancova", "art", "permuco", "permuco_ancova", "nonparam_anova", "rank", "brunner_munzel", "ats", "t", "welch", "rank_ancova")
    if (!is.character(method) || !length(method) || anyNA(method) || any(!method %in% choices)) {
      stop("methodには次の方法を指定してください: ", paste(choices, collapse = ", "), call. = FALSE)
    }
    expanded <- unique(unlist(lapply(method, function(x) {
      if (x == "nonparam_anova") c("art", "permuco") else x
    }), use.names = FALSE))
    if (sum(c("ancova", "permuco_ancova", "rank_ancova") %in% expanded) > 1L) {
      stop('共変量の分析はANCOVA・置換ANCOVA・順位ANCOVAから1つを選択してください（ancova / permuco_ancova / rank_ancova）。',
        call. = FALSE)
    }
    expanded
  }

  gui_enabled <- function(experiment, methods, settings = list(posthoc = TRUE, followup_scope = "followup")) {
    normal <- any(methods %in% c("anova", "ancova"))
    can_followup <- normal || "art" %in% methods || any(methods %in% c("ats", "rank", "welch", "rank_ancova")) &&
      length(experiment$iv) == 1L && nlevels(factor(experiment$data[[experiment$iv]])) > 2L
    enabled <- isTRUE(settings$posthoc) && can_followup
    followup <- enabled && length(adjustment_policy(methods)$configurable) > 0L
    permutation <- any(methods %in% c("permuco", "permuco_ancova"))
    within <- length(experiment$within) > 0L
    between <- length(experiment$between) > 0L
    list(alpha = TRUE, trend_alpha = TRUE, posthoc = can_followup,
      followup_scope = enabled, display_scope = enabled, adjust = followup,
      simple_effects = enabled && (normal || "art" %in% methods) && length(experiment$iv) > 1L && settings$followup_scope != "all",
      correction = (normal || "art" %in% methods) && within,
      type = normal, effect_size = normal, np = permutation, seed = permutation || "ats" %in% methods,
      permutation_method = permutation,
      covariate_interactions = ("ancova" %in% methods && between) ||
        ("permuco_ancova" %in% methods && (within || between)),
      digits = TRUE, p_digits = TRUE, show = TRUE)
  }

  gui_choices <- function(experiment, initial) {
    adjust <- adjustment_policy(initial$method)$choices
    permutation <- if (length(experiment$within)) {
      c("Rd_kheradPajouh_renaud", "Rde_kheradPajouh_renaud")
    } else {
      c("freedman_lane", "draper_stoneman", "manly", "terBraak", "kennedy", "huh_jhun", "dekker")
    }
    list(followup_scope = list("有意な効果のみ" = "significant", "有意＋有意傾向の効果" = "followup", "すべて" = "all"),
      display_scope = list("有意な比較のみ" = "significant", "有意な効果の全比較" = "effect", "計算済みの全結果" = "all"),
      adjust = setNames(as.list(adjust), adjust),
      correction = list("自動（Mauchly・ε → GG/HF）" = "auto", "GG" = "GG", "HF" = "HF", "なし" = "none"),
      type = list("II" = 2, "III" = 3),
      effect_size = list("偏η²" = "pes", "偏η²＋一般化η²" = "ges"),
      permutation_method = c(list("自動選択" = NULL), setNames(as.list(permutation), permutation)),
      covariate_interactions = list("方法ごとの既定" = NULL, "含める" = TRUE, "含めない" = FALSE))
  }

  validate_gui_settings <- function(experiment, settings) {
    methods <- expand_methods(settings$method)
    if (isTRUE(settings$posthoc)) validate_adjustment(methods, settings$adjust)
    validate_common_settings(settings$alpha, settings$followup_alpha, settings$digits, settings$p_digits,
      settings$correction, settings$sphericity)
    validate_scope_settings(settings$alpha, settings$trend_alpha, settings$followup_alpha, settings$posthoc, settings$followup_scope, settings$display_scope)
    checked_adjust <- if (settings$posthoc) settings$adjust else "bonferroni"
    if ("rank" %in% methods) validate_rank_settings(experiment, checked_adjust)
    for (method in intersect(methods, c("brunner_munzel", "ats"))) .t_relative_design(experiment, method)
    for (method in intersect(methods, c("t", "welch"))) .t_mean_design(experiment, method)
    if ("rank_ancova" %in% methods) validate_rank_ancova_settings(experiment, checked_adjust, settings$covariate_interactions)
    if (any(methods %in% c("ancova", "permuco_ancova")) && !length(experiment$covariate)) {
      stop("ANCOVAには、前処理で指定した共変量が必要です。", call. = FALSE)
    }
    scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
    for (name in c("digits", "p_digits")) {
      value <- settings[[name]]
      if (!scalar(value) || value != floor(value) || value < if (name == "digits") 0 else 1) {
        stop(if (name == "digits") "数値の桁数は0以上の整数にしてください。" else
          "p値の桁数は1以上の整数にしてください。", call. = FALSE)
      }
    }
    for (name in c("show", "simple_effects", "posthoc")) {
      if (!is.logical(settings[[name]]) || length(settings[[name]]) != 1L || is.na(settings[[name]])) {
        stop(name, "にはTRUEまたはFALSEを指定してください。", call. = FALSE)
      }
    }
    active <- gui_enabled(experiment, methods, settings)
    if (active$np && (!scalar(settings$np) || settings$np < 2 || settings$np != floor(settings$np))) {
      stop("置換回数は2以上の整数にしてください。", call. = FALSE)
    }
    if (active$seed && !is.null(settings$seed) && (!scalar(settings$seed) ||
        settings$seed < 0 || settings$seed > .Machine$integer.max || settings$seed != floor(settings$seed))) {
      stop("乱数シードは空欄、または0以上の整数にしてください。", call. = FALSE)
    }
    maps <- gui_choices(experiment, settings)
    for (name in names(maps)) {
      if (!isTRUE(active[[name]])) next
      if (!any(vapply(maps[[name]], function(x) isTRUE(all.equal(x, settings[[name]])), logical(1)))) {
        stop(name, "の選択値が利用できません。", call. = FALSE)
      }
    }
    invisible(settings)
  }

  parse_gui_settings <- function(experiment, initial, methods, text) {
    settings <- initial
    settings$method <- methods
    if ("rank_ancova" %in% methods) settings$covariate_interactions <- FALSE
    settings$posthoc <- identical(text$posthoc, "1")
    if (settings$posthoc) settings$followup_scope <- gui_choices(experiment, settings)$followup_scope[[text$followup_scope]]
    active <- gui_enabled(experiment, methods, settings)
    maps <- gui_choices(experiment, settings)
    for (name in names(active)) {
      if (!active[[name]]) next
      value <- text[[name]]
      if (name %in% names(maps)) {
        if (!value %in% names(maps[[name]])) stop(name, "を選択してください。", call. = FALSE)
        value <- maps[[name]][[value]]
      } else if (name %in% c("show", "simple_effects", "posthoc")) {
        if (!value %in% c("0", "1")) stop(name, "の設定を確認してください。", call. = FALSE)
        value <- value == "1"
      } else {
        value <- trimws(value)
        value <- if (name == "seed" && !nzchar(value)) NULL else suppressWarnings(as.numeric(value))
      }
      settings[name] <- list(value)
    }
    settings$followup_alpha <- settings$trend_alpha
    validate_gui_settings(experiment, settings)
    settings
  }

  choose_analysis <- function(experiment, initial) {
    if (!is.list(experiment) || !length(experiment$dv) || !length(experiment$iv)) {
      stop("前処理で作成したexperimentを指定してください。", call. = FALSE)
    }
    initial$trend_alpha <- resolve_shared_threshold(initial$trend_alpha, initial$followup_alpha)
    initial$followup_alpha <- initial$trend_alpha
    defaults <- list(posthoc = TRUE, followup_scope = "followup", display_scope = "significant")
    for (name in names(defaults)) if (is.null(initial[[name]])) initial[[name]] <- defaults[[name]]
    initial$method <- expand_methods(initial$method)
    initial$adjust <- normalize_adjustment(initial$adjust)
    if (!is.null(initial$sphericity)) {
      initial$correction <- initial$sphericity
      initial$sphericity <- NULL
    }
    validate_gui_settings(experiment, initial)
    if (!isTRUE(capabilities("tcltk")) ||
        !suppressWarnings(requireNamespace("tcltk", quietly = TRUE))) {
      stop("GUIにはRのtcltkが必要です。利用できない環境ではselect = \"arguments\"で実行してください。", call. = FALSE)
    }
    window <- tryCatch(tcltk::tktoplevel(), error = function(e) {
      stop("GUIを開けませんでした。デスクトップ版のR/RStudioで実行するか、select = \"arguments\"で実行してください。\n",
        conditionMessage(e), call. = FALSE)
    })
    on.exit({
      try(tcltk::tkgrab.release(window), silent = TRUE)
      try(tcltk::tkdestroy(window), silent = TRUE)
    }, add = TRUE)
    tcltk::tkwm.title(window, "検定・分析の設定")
    tcltk::tkwm.minsize(window, 700, 560)
    tcltk::tkgrid.columnconfigure(window, 0, weight = 1)
    tcltk::tkgrid.rowconfigure(window, 2, weight = 1)
    show_names <- function(x) if (length(x)) paste(x, collapse = ", ") else "なし"
    summary <- paste0("従属変数：", show_names(experiment$dv),
      "\n被験者内：", show_names(experiment$within), " / 被験者間：", show_names(experiment$between),
      " / 共変量：", show_names(experiment$covariate))
    tcltk::tkgrid(tcltk::tklabel(window, text = summary, anchor = "w", justify = "left", wraplength = 740),
      row = 0, column = 0, sticky = "ew", padx = 12, pady = 8)
    methods_panel <- tcltk::tkframe(window)
    tcltk::tkgrid(methods_panel, row = 1, column = 0, sticky = "ew", padx = 12, pady = 4)
    method_names <- c(t = "t検定（対応あり／なし・2水準）", anova = "ANOVA",
      welch = "Welch（対応なし・2群以上）", art = "ART（下位検定：ART-C）", permuco = "permuco",
      rank = "ランク", brunner_munzel = "Brunner–Munzel（対応なし2群）", ats = "ATS＋MCTP（対応なし3群以上）")
    method_positions <- c(0L, 1L, 3L, 4L, 5L, 6L, 7L, 8L)
    available <- vapply(names(method_names), function(name) {
      if (name == "rank") return(length(experiment$iv) == 1L)
      if (name %in% c("t", "welch")) return(tryCatch({ .t_mean_design(experiment, name); TRUE }, error = function(e) FALSE))
      if (!name %in% c("brunner_munzel", "ats")) return(TRUE)
      tryCatch({ .t_relative_design(experiment, name); TRUE }, error = function(e) FALSE)
    }, logical(1))
    method_vars <- list()
    for (i in seq_along(method_names)) {
      name <- names(method_names)[i]
      method_vars[[name]] <- tcltk::tclVar(if (name %in% initial$method &&
        available[[name]]) "1" else "0")
      checkbox <- tcltk::tkcheckbutton(methods_panel, text = method_names[[name]],
        variable = method_vars[[name]], command = function() refresh())
      if (!available[[name]]) {
        tcltk::tkconfigure(checkbox, state = "disabled")
      }
      tcltk::tkgrid(checkbox, row = method_positions[i] %/% 3L, column = method_positions[i] %% 3L, sticky = "w", padx = 4)
    }
    # 共変量の分析は独立した単一選択。共変量なしの各分析は従来どおり併用できる。
    covariate_method <- intersect(initial$method, c("ancova", "permuco_ancova", "rank_ancova"))
    covariate_var <- tcltk::tclVar(if (length(covariate_method)) covariate_method else "none")
    covariate_panel <- tcltk::tkframe(methods_panel)
    tcltk::tkgrid(covariate_panel, row = ceiling(length(method_names) / 3), column = 0, columnspan = 3, sticky = "ew", pady = 6)
    tcltk::tkgrid(tcltk::tklabel(covariate_panel, text = "共変量の分析（単一選択）", anchor = "w"),
      row = 0, column = 0, columnspan = 3, sticky = "w")
    covariate_names <- c(none = "実行しない", ancova = "ANCOVA", permuco_ancova = "置換ANCOVA（permuco）",
      rank_ancova = "順位ANCOVA（Quade・対応なし1要因）")
    for (i in seq_along(covariate_names)) {
      name <- names(covariate_names)[i]
      radio <- tcltk::tkradiobutton(covariate_panel, text = covariate_names[[name]],
        variable = covariate_var, value = name, command = function() refresh())
      if (name != "none" && !length(experiment$covariate)) tcltk::tkconfigure(radio, state = "disabled")
      if (name == "rank_ancova" && !tryCatch({ .t_rank_ancova_design(experiment); TRUE }, error = function(e) FALSE)) {
        tcltk::tkconfigure(radio, state = "disabled")
      }
      tcltk::tkgrid(radio, row = i, column = 0, sticky = "w", padx = 4)
    }
    notebook <- tcltk::ttknotebook(window)
    tcltk::tkgrid(notebook, row = 2, column = 0, sticky = "nsew", padx = 12, pady = 6)
    panels <- list()
    for (name in c("基本設定", "モデル・置換検定", "表示設定")) {
      panel <- tcltk::ttkframe(notebook, padding = 12)
      tcltk::tkadd(notebook, panel, text = name)
      tcltk::tkgrid.columnconfigure(panel, 1, weight = 1)
      panels[[name]] <- panel
    }
    maps <- gui_choices(experiment, initial)
    fields <- variables <- list()
    add_field <- function(panel, row, name, label, kind = "entry") {
      value <- initial[[name]]
      if (name %in% names(maps)) {
        matches <- vapply(maps[[name]], function(x) isTRUE(all.equal(x, value)), logical(1))
        value <- if (any(matches)) names(maps[[name]])[which(matches)[1L]] else names(maps[[name]])[1L]
      } else if (kind == "check") {
        value <- if (isTRUE(value)) "1" else "0"
      } else value <- if (is.null(value)) "" else as.character(value)
      variable <- tcltk::tclVar(value)
      variables[[name]] <<- variable
      tcltk::tkgrid(tcltk::tklabel(panel, text = label, anchor = "w"),
        row = row, column = 0, sticky = "w", padx = 4, pady = 3)
      widget <- if (name %in% names(maps)) {
        tcltk::ttkcombobox(panel, textvariable = variable, values = names(maps[[name]]), state = "readonly", width = 34)
      } else if (kind == "check") {
        tcltk::tkcheckbutton(panel, variable = variable, text = "有効にする", command = function() refresh())
      } else tcltk::tkentry(panel, textvariable = variable, width = 34)
      tcltk::tkgrid(widget, row = row, column = 1, sticky = "ew", padx = 4, pady = 3)
      fields[[name]] <<- list(widget = widget, kind = if (name %in% names(maps)) "choice" else kind)
    }
    add_field(panels[["基本設定"]], 0, "alpha", "有意水準")
    add_field(panels[["基本設定"]], 1, "trend_alpha", "有意傾向の上限")
    add_field(panels[["基本設定"]], 4, "posthoc", "下位検定を実行する", "check")
    add_field(panels[["基本設定"]], 5, "adjust", "多重比較補正法")
    add_field(panels[["基本設定"]], 6, "followup_scope", "下位検定の計算範囲")
    add_field(panels[["基本設定"]], 2, "correction", "球面性補正")
    tcltk::tkgrid(tcltk::tklabel(panels[["基本設定"]], anchor = "w", justify = "left", wraplength = 650,
      text = "自動：Mauchlyで棄却 → GGのε ≤ .75ならGG、> .75ならHF。\nARTへのGG/HF適用は、このスクリプト独自の拡張です。"),
      row = 3, column = 0, columnspan = 2, sticky = "w", padx = 4, pady = 6)
    add_field(panels[["モデル・置換検定"]], 0, "type", "検定タイプ（ANOVA / ANCOVA）")
    add_field(panels[["モデル・置換検定"]], 1, "effect_size", "効果量（ANOVA / ANCOVA）")
    add_field(panels[["モデル・置換検定"]], 2, "covariate_interactions", "共変量との交互作用")
    add_field(panels[["モデル・置換検定"]], 3, "np", "置換回数")
    add_field(panels[["モデル・置換検定"]], 4, "seed", "乱数シード（空欄＝未指定）")
    add_field(panels[["モデル・置換検定"]], 5, "permutation_method", "置換法")
    tcltk::tkgrid(tcltk::tklabel(panels[["モデル・置換検定"]],
      text = "通常ANCOVAの被験者内要因×共変量は常に含みます。\nART・permucoの検定タイプは各方法の仕様に従います。\nランクは1要因の比較に対応しています。\nBrunner–Munzel・ATS＋MCTPは対応なし・1要因のみです。\nATS＋MCTPの補正は固定で、乱数シードは数値積分に使用します。\nWelchの3群以上の下位検定はGames–Howellに固定です。\n順位ANCOVAはQuade法（対応なし1要因・共変量との交互作用なし）です。",
      anchor = "w", justify = "left", wraplength = 650),
      row = 6, column = 0, columnspan = 2, sticky = "w", padx = 4, pady = 8)
    add_field(panels[["表示設定"]], 0, "digits", "数値の小数点以下の桁数")
    add_field(panels[["表示設定"]], 1, "p_digits", "p値の小数点以下の桁数")
    add_field(panels[["表示設定"]], 2, "show", "コンソールに結果を表示", "check")
    add_field(panels[["表示設定"]], 3, "display_scope", "下位検定の表示範囲")
    add_field(panels[["表示設定"]], 4, "simple_effects", "主効果に基づく単純比較の追加", "check")
    tcltk::tkgrid(tcltk::tklabel(panels[["表示設定"]], text = "単純比較の追加は、表示だけでなく計算対象も増やします。",
      anchor = "w", justify = "left", wraplength = 650), row = 5, column = 0, columnspan = 2, sticky = "w", padx = 4, pady = 6)
    tcltk::tkgrid(tcltk::tklabel(panels[["基本設定"]], text = "共通補正の選択肢は選択した方法の共通部分です。\nGames–Howell/MCTPは固有調整。BH/BYはFDR制御です。",
      anchor = "w", justify = "left", wraplength = 650), row = 7, column = 0, columnspan = 2, sticky = "w", padx = 4, pady = 6)
    get_methods <- function() {
      selected <- names(method_vars)[vapply(method_vars, function(v) tcltk::tclvalue(v) == "1", logical(1))]
      covariate <- tcltk::tclvalue(covariate_var)
      c(selected, if (covariate != "none") covariate)
    }
    refresh <- function() {
      selected <- get_methods()
      settings <- initial
      settings$method <- selected
      settings$posthoc <- tcltk::tclvalue(variables$posthoc) == "1"
      settings$followup_scope <- maps$followup_scope[[tcltk::tclvalue(variables$followup_scope)]]
      active <- gui_enabled(experiment, selected, settings)
      settings$adjust <- tcltk::tclvalue(variables$adjust)
      maps$adjust <<- gui_choices(experiment, settings)$adjust
      tcltk::tkconfigure(fields$adjust$widget, values = names(maps$adjust))
      if (!tcltk::tclvalue(variables$adjust) %in% names(maps$adjust)) tcltk::tclvalue(variables$adjust) <- "bonferroni"
      for (name in names(fields)) {
        field <- fields[[name]]
        state <- if (!active[[name]]) "disabled" else if (field$kind == "choice") "readonly" else "normal"
        tcltk::tkconfigure(field$widget, state = state)
      }
      invisible(NULL)
    }
    tcltk::tkbind(fields$followup_scope$widget, "<<ComboboxSelected>>", function() refresh())
    refresh()
    error_text <- tcltk::tclVar("")
    tcltk::tkgrid(tcltk::tklabel(window, textvariable = error_text, foreground = "#a32020",
      anchor = "w", justify = "left", wraplength = 740),
      row = 3, column = 0, sticky = "ew", padx = 12, pady = 4)
    result <- NULL
    confirm <- function() {
      settings <- tryCatch(parse_gui_settings(experiment, initial, get_methods(),
        lapply(variables, tcltk::tclvalue)), error = function(e) {
          tcltk::tclvalue(error_text) <- conditionMessage(e)
          NULL
        })
      if (is.null(settings)) return(invisible(NULL))
      result <<- settings
      tcltk::tkdestroy(window)
      invisible(NULL)
    }
    cancel <- function() tcltk::tkdestroy(window)
    buttons <- tcltk::tkframe(window)
    tcltk::tkgrid(buttons, row = 4, column = 0, sticky = "e", padx = 12, pady = 10)
    tcltk::tkpack(tcltk::tkbutton(buttons, text = "キャンセル", command = cancel), side = "left", padx = 4)
    tcltk::tkpack(tcltk::tkbutton(buttons, text = "実行", command = confirm, default = "active"), side = "left", padx = 4)
    tcltk::tkwm.protocol(window, "WM_DELETE_WINDOW", cancel)
    tcltk::tkbind(window, "<Escape>", cancel)
    if (!as.logical(as.integer(tcltk::tkwinfo("viewable", window)))) tcltk::tkwait.visibility(window)
    tcltk::tkgrab.set(window)
    tcltk::tkwait.window(window)
    if (is.null(result)) stop("設定をキャンセルしました。分析は実行していません。", call. = FALSE)
    result
  }
}, envir = .test_private_environment)


# ------------------------------------------------------------
# 公開関数：分析方法を選択して実行
# analyze(experiment) または select = "gui" で設定画面を開く。
# methodを指定した呼出し、またはselect = "arguments"は引数で直接実行する。
# method = c("art", "permuco") で両方、単一指定で従来の結果形式。
# nonparam_anova は c("art", "permuco") の別名。
# ancova / permuco_ancova / rank_ancova は experiment$covariate を使用する。
# この3方法は排他選択。GUIは単一選択、引数での同時指定は実行前に停止する。
# covariate_interactions = NULL は既定を維持（ancova: FALSE, permuco: TRUE）。
# TRUE/FALSEの明示で群間要因×共変量を指定できる。
# ancovaの被験者内要因×共変量はafexの標準モデルに従い常に含む。
#
# posthocで実行可否、followup_scopeで計算範囲、display_scopeで表示範囲を指定。
# trend_alphaを有意傾向の判定と下位検定の共通基準に使用。followup_alphaは旧互換名。
# simple_effects = TRUE：
#   上記に加え、2要因交互作用が基準外でも、構成要因の主効果が
#   選択した計算基準内なら、その方向の単純比較を追加する。
# ------------------------------------------------------------

analyze <- function(
  experiment,
  method = "anova",
  alpha = 0.05,
  adjust = "bonferroni",
  digits = 3,
  p_digits = 3,
  followup_alpha = NULL, # 旧名。trend_alphaと同じ値を指す互換引数。
  simple_effects = FALSE,
  correction = "auto",
  effect_size = "pes",
  type = 3,
  show = TRUE,
  sphericity = NULL, # 旧ART専用指定。通常は共通のcorrectionを使用。
  np = 5000,
  seed = NULL,
  permutation_method = NULL,
  covariate_interactions = NULL,
  select = if (missing(method)) "gui" else "arguments",
  trend_alpha = .10,
  posthoc = TRUE,
  followup_scope = "followup",
  display_scope = "significant"
) {
  trend_alpha <- .test_private_environment$resolve_shared_threshold(if (missing(trend_alpha)) NULL else trend_alpha, followup_alpha)
  followup_alpha <- trend_alpha
  select <- match.arg(select, c("gui", "arguments"))
  if (select == "gui") {
    experiment_expression <- substitute(experiment)
    initial <- list(method = method, alpha = alpha, adjust = adjust, digits = digits,
      p_digits = p_digits, followup_alpha = followup_alpha, simple_effects = simple_effects,
      correction = correction, effect_size = effect_size, type = type, show = show,
      sphericity = sphericity, np = np, seed = seed, permutation_method = permutation_method,
      covariate_interactions = covariate_interactions, trend_alpha = trend_alpha,
      posthoc = posthoc, followup_scope = followup_scope, display_scope = display_scope)
    settings <- .test_private_environment$choose_analysis(experiment, initial)
    settings$followup_alpha <- NULL # 再実行コードには共通のtrend_alphaだけを出す。
    replay_call <- as.call(c(list(as.name("analyze")), list(experiment = experiment_expression),
      settings, list(select = "arguments")))
    code <- deparse(replay_call, width.cutoff = 80L)
    code[1L] <- paste0("result <- ", code[1L])
    code <- paste(code, collapse = "\n")
    if (settings$show) cat("\n【再実行用コード】\n", code, "\n", sep = "")
    result <- do.call(analyze, c(list(experiment = experiment), settings, list(select = "arguments")))
    attr(result, "selection") <- list(method = "gui", call = replay_call, code = code)
    return(invisible(result))
  }
  methods <- .test_private_environment$expand_methods(method)
  adjust <- .test_private_environment$normalize_adjustment(adjust)
  if (isTRUE(posthoc)) .test_private_environment$validate_adjustment(methods, adjust)
  .test_private_environment$validate_common_settings(alpha, followup_alpha, digits, p_digits, correction, sphericity)
  checked_adjust <- if (isTRUE(posthoc)) adjust else "bonferroni"
  if ("rank" %in% methods) .test_private_environment$validate_rank_settings(experiment, checked_adjust)
  for (selected in intersect(methods, c("brunner_munzel", "ats"))) {
    .test_private_environment$.t_relative_design(experiment, selected)
  }
  for (selected in intersect(methods, c("t", "welch"))) .test_private_environment$.t_mean_design(experiment, selected)
  if ("rank_ancova" %in% methods) .test_private_environment$validate_rank_ancova_settings(experiment, checked_adjust, covariate_interactions)
  for (flag in list(show, simple_effects)) {
    if (!is.logical(flag) || length(flag) != 1L || is.na(flag)) stop("showとsimple_effectsにはTRUEまたはFALSEを指定してください。", call. = FALSE)
  }
  if (!is.null(covariate_interactions) && (!is.logical(covariate_interactions) ||
      length(covariate_interactions) != 1L || is.na(covariate_interactions))) {
    stop("covariate_interactionsにはNULL、TRUE、FALSEを指定してください。", call. = FALSE)
  }
  if (!length(experiment$dv)) stop("従属変数を指定してください。", call. = FALSE)
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
                        seed < 0 || seed > .Machine$integer.max || seed != floor(seed))) {
    stop("seedには0以上の整数を指定してください。", call. = FALSE)
  }
  .test_private_environment$validate_scope_settings(alpha, trend_alpha, followup_alpha, posthoc, followup_scope, display_scope)
  engine_followup_alpha <- trend_alpha
  prepare_current <- function(raw, method) {
    result <- .test_private_environment$prepare_result(raw, method,
      alpha, trend_alpha, simple_effects, sphericity, adjust, correction, posthoc, followup_scope, display_scope)
    result
  }
  collected <- setNames(lapply(methods, function(x) list()), methods)
  rng_states <- list()
  warnings_by_method <- notices_by_method <- setNames(lapply(methods, function(x) character()), methods)
  shown_warnings <- shown_notices <- warnings_by_method

  # 同じメソッドの乱数列を従属変数間で引き継ぎ、明示seed時は呼出元の乱数状態を復元する。
  run_engine <- function(method, experiment) {
    aggregation <- NULL
    if (method %in% c("art", "permuco", "permuco_ancova")) {
      prepared <- .test_private_environment$aggregate_observation_units(experiment, method)
      experiment <- prepared$experiment
      aggregation <- prepared$info
    }
    analysis_arguments <- switch(
    method,
    anova = list(
      experiment = experiment,
      alpha = alpha,
      correction = correction,
      adjust = adjust,
      digits = digits,
      p_digits = p_digits,
      effect_size = effect_size,
      type = type,
      followup_alpha = engine_followup_alpha,
      simple_effects = simple_effects
    ),
    ancova = list(
      experiment = experiment, alpha = alpha, correction = correction,
      adjust = adjust, digits = digits, p_digits = p_digits, effect_size = effect_size,
      type = type, followup_alpha = engine_followup_alpha, simple_effects = simple_effects,
      covariate_interactions = if (is.null(covariate_interactions)) FALSE else covariate_interactions
    ),
    art = list(
      experiment = experiment,
      alpha = alpha,
      sphericity = "none", # 補正の適用はhelperのみで行う。
      adjust = adjust,
      digits = digits,
      p_digits = p_digits,
      followup_alpha = engine_followup_alpha,
      simple_effects = simple_effects
    ),
    rank = list(experiment = experiment, alpha = alpha, followup_alpha = engine_followup_alpha),
    rank_ancova = list(experiment = experiment, alpha = alpha, followup_alpha = engine_followup_alpha,
      covariate_interactions = if (is.null(covariate_interactions)) FALSE else covariate_interactions),
    t = list(experiment = experiment, alpha = alpha, followup_alpha = engine_followup_alpha),
    welch = list(experiment = experiment, alpha = alpha, followup_alpha = engine_followup_alpha),
    brunner_munzel = list(experiment = experiment, alpha = alpha, followup_alpha = engine_followup_alpha),
    ats = list(experiment = experiment, alpha = alpha, followup_alpha = engine_followup_alpha, seed = seed),
    permuco = list(
      experiment = experiment, alpha = alpha, followup_alpha = engine_followup_alpha,
      digits = digits, p_digits = p_digits, np = np, seed = seed,
      permutation_method = permutation_method, use_covariates = FALSE
    ),
    permuco_ancova = list(
      experiment = experiment, alpha = alpha, followup_alpha = engine_followup_alpha,
      digits = digits, p_digits = p_digits, np = np, seed = seed,
      permutation_method = permutation_method, use_covariates = TRUE,
      covariate_interactions = if (is.null(covariate_interactions)) TRUE else covariate_interactions
    )
  )
    analysis_function <- get(switch(method, anova = "t_anova", ancova = "t_ancova", art = "t_art",
      permuco = "t_permuco", permuco_ancova = "t_permuco", rank = "t_rank",
      brunner_munzel = "t_brunner_munzel", ats = "t_ats",
      t = "t_t", welch = "t_welch", rank_ancova = "t_rank_ancova"), envir = .test_private_environment, inherits = FALSE)
    if (method %in% c("permuco", "permuco_ancova") && !is.null(seed)) {
      had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      if (had_seed) previous_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit({
        if (had_seed) assign(".Random.seed", previous_seed, envir = .GlobalEnv)
        else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(list = ".Random.seed", envir = .GlobalEnv)
      }, add = TRUE)
      if (is.null(rng_states[[method]])) set.seed(seed)
      else assign(".Random.seed", rng_states[[method]], envir = .GlobalEnv)
      analysis_arguments$seed <- NULL
      analysis_arguments$.followup_seed <- seed
    }
    raw <- do.call(analysis_function, analysis_arguments)
    if (!is.null(aggregation)) raw <- lapply(raw, function(result) {
      result$cell_replication <- aggregation
      result$messages <- unique(c(result$messages, sprintf(
        "同一ID×条件の従属変数を平均しました（%d行→%dセル、セルあたり%d〜%d行）。",
        aggregation$n_observations, aggregation$n_participant_cells,
        aggregation$repetitions_min, aggregation$repetitions_max)))
      result
    })
    if (method %in% c("permuco", "permuco_ancova") && !is.null(seed)) {
      rng_states[[method]] <<- get(".Random.seed", envir = .GlobalEnv)
      raw <- lapply(raw, function(result) { result$seed <- seed; result })
    }
    raw
  }
  invoke <- function(method, expression, followup = FALSE) {
    # 結果が得られない場合にも、保留中の通知・警告は失わない。
    unavailable <- function() {
      if (!show) return(invisible(NULL))
      .test_private_environment$print_notes(
        messages = setdiff(notices_by_method[[method]], shown_notices[[method]]),
        warnings = setdiff(warnings_by_method[[method]], shown_warnings[[method]]),
        label = .test_private_environment$method_label(method, followup))
      shown_warnings[[method]] <<- union(shown_warnings[[method]], warnings_by_method[[method]])
      shown_notices[[method]] <<- union(shown_notices[[method]], notices_by_method[[method]])
    }
    value <- withCallingHandlers(force(expression), warning = function(w) {
      message <- conditionMessage(w)
      warnings_by_method[[method]] <<- unique(c(warnings_by_method[[method]], message))
      invokeRestart("muffleWarning")
    }, message = function(m) {
      message <- conditionMessage(m)
      notices_by_method[[method]] <<- unique(c(notices_by_method[[method]], message))
      invokeRestart("muffleMessage")
    }, error = function(e) unavailable())
    if (!length(value)) unavailable()
    value
  }
  emit <- function(values, stage) {
    if (!show || !length(values)) return(invisible(NULL))
    if (stage %in% c("anova", "ancova", "rank", "t", "followup")) for (method in names(values)) {
      warnings <- unique(c(values[[method]]$warnings, warnings_by_method[[method]]))
      notices <- unique(c(values[[method]]$messages, notices_by_method[[method]]))
      values[[method]]$warnings <- setdiff(warnings, shown_warnings[[method]])
      values[[method]]$messages <- setdiff(notices, shown_notices[[method]])
      shown_warnings[[method]] <<- union(shown_warnings[[method]], warnings)
      shown_notices[[method]] <<- union(shown_notices[[method]], notices)
    }
    collection <- lapply(values, function(result) setNames(list(result), result$dependent_variable))
    .test_private_environment$display_results(collection, digits, p_digits, followup_alpha, stage = stage)
    utils::flush.console()
    invisible(NULL)
  }
  for (target in experiment$dv) {
    if (show) .test_private_environment$target_heading(target)
    current_experiment <- experiment
    current_experiment$dv <- target
    warnings_by_method <- notices_by_method <- setNames(lapply(methods, function(x) character()), methods)
    shown_warnings <- shown_notices <- warnings_by_method
    current <- list()
    # 球面性評価に対応するモデルを先に計算する。
    for (method in intersect(methods, c("anova", "ancova", "art"))) {
      raw <- invoke(method, run_engine(method, current_experiment))
      if (length(raw)) current[[method]] <- prepare_current(raw[[target]], method)
    }
    emit(current, "sphericity")
    for (method in intersect(methods, c("permuco", "rank", "brunner_munzel", "ats", "t", "welch"))) {
      raw <- invoke(method, run_engine(method, current_experiment))
      if (length(raw)) current[[method]] <- raw[[target]]
    }
    for (method in setdiff(names(current), c("anova", "ancova", "art"))) {
      current[[method]] <- prepare_current(current[[method]], method)
    }
    current <- current[intersect(methods, names(current))]
    emit(current, "diagnostics")
    emit(current[.test_private_environment$mean_t_methods(current)], "t")
    emit(current[.test_private_environment$anova_methods(current)], "anova")
    # 共変量込みの各方法をまとめ、方法ごとの縦ブロックで表示する。
    if ("permuco_ancova" %in% methods) {
      raw <- invoke("permuco_ancova", run_engine("permuco_ancova", current_experiment))
      if (length(raw)) {
        current$permuco_ancova <- prepare_current(raw[[target]], "permuco_ancova")
      }
    }
    current <- current[intersect(methods, names(current))]
    if ("rank_ancova" %in% methods) {
      raw <- invoke("rank_ancova", run_engine("rank_ancova", current_experiment))
      if (length(raw)) current$rank_ancova <- prepare_current(raw[[target]], "rank_ancova")
    }
    current <- current[intersect(methods, names(current))]
    emit(current[intersect(names(current), c("ancova", "permuco_ancova", "rank_ancova"))], "ancova")
    emit(current[intersect(names(current), c("rank", "brunner_munzel", "ats"))], "rank")
    effects <- unique(unlist(lapply(current, `[[`, "follow_up_effects"), use.names = FALSE))
    intro_printed <- FALSE
    for (effect in effects) {
      completed <- list()
      for (method in names(current)) {
        result <- current[[method]]
        if (!effect %in% result$follow_up_effects) next
        before <- names(result$contrasts)
        current[[method]] <- invoke(method, .test_private_environment$finalize_effect(result, effect, method, adjust), followup = TRUE)
        displayed <- current[[method]]
        displayed$contrasts <- displayed$contrasts[setdiff(names(displayed$contrasts), before)]
        completed[[method]] <- displayed
      }
      if (!intro_printed) {
        emit(current, "followup_intro")
        for (method in names(current)) {
          current[[method]]$header_p_adjustment <- unique(unlist(lapply(current[[method]]$contrasts, `[[`, "p_adjustment"), use.names = FALSE))
          current[[method]]$header_ci_adjustment <- setdiff(unique(unlist(lapply(current[[method]]$contrasts, `[[`, "ci_adjustment"), use.names = FALSE)), "not_computed")
        }
        intro_printed <- TRUE
      }
      for (method in names(completed)) {
        completed[[method]]$header_p_adjustment <- current[[method]]$header_p_adjustment
        completed[[method]]$header_ci_adjustment <- current[[method]]$header_ci_adjustment
      }
      emit(completed, "followup")
    }
    if (!intro_printed && (any(names(current) %in% c("anova", "ancova", "art", "ats")) ||
        ("welch" %in% names(current) && current$welch$n_levels > 2L) ||
        ("rank_ancova" %in% names(current) && current$rank_ancova$n_levels > 2L) ||
        ("rank" %in% names(current) && current$rank$n_levels > 2L))) emit(current, "followup_intro")
    for (method in names(current)) {
      result <- current[[method]]
      result$compute_followups <- NULL
      result$header_p_adjustment <- result$header_ci_adjustment <- NULL
      result$analysis_method <- method
      result$warnings <- unique(c(result$warnings, warnings_by_method[[method]]))
      result$messages <- unique(c(result$messages, notices_by_method[[method]]))
      class(result) <- unique(c("analysis_result", "analysis_test_result", class(result)))
      collected[[method]][[target]] <- result
    }
  }
  for (method in methods) {
    results <- collected[[method]]
    attr(results, "stage") <- "test"
    attr(results, "type") <- "analysis"
    attr(results, "method") <- method
    attr(results, "experiment") <- experiment
    attr(results, "interface") <- "analyze"
    class(results) <- c("analysis_results", "analysis_test", "list")
    collected[[method]] <- results
  }
  if (length(methods) == 1L) return(invisible(collected[[1L]]))
  attr(collected, "method") <- methods
  attr(collected, "experiment") <- experiment
  class(collected) <- c("analysis_multi_results", "list")
  invisible(collected)

}


# ------------------------------------------------------------
# helperの読み込みにだけ使用した変数を削除
# .test_private_environment は公開関数が参照するため保持する。
# ------------------------------------------------------------

rm(
  list = c(
    ".test_source_file",
    ".test_frames",
    ".test_index",
    ".test_candidate",
    ".test_directory",
    ".test_child_paths",
    ".test_missing_scripts",
    ".test_path",
    ".test_required_functions",
    ".test_loaded",
    ".test_package_groups",
    ".test_package_status",
    ".test_missing_packages",
    ".test_helper_quiet",
    ".test_group_name",
    ".test_group_missing"
  ),
  envir = environment()
)
