# ============================================================
# 通常ANOVAの推定・効果量・診断・未補正対比。制御・補正・表示はhelperが担当。
# t_ancovaからも共通計算部として使用する。
#
# 公開関数:
#   t_anova()
#
# 必要パッケージ:
#   afex
#   emmeans
#   car
#
# 効果量：
#   偏イータ二乗（pes）を常に分散分析表へ保存します。
#   effect_size = "ges"の場合は一般化イータ二乗（ges）も追加します。
#
# 下位検定：
#   p < followup_alpha の効果（有意・有意傾向）を対象にします。
#   simple_effects = TRUE の場合は、基準外の2要因交互作用でも、
#   構成要因の主効果が有意・有意傾向なら、その方向だけ追加します。
#
# experimentに必要な要素:
#   data, id, iv, dv, within, between, design
# ============================================================

t_anova <- function(
  experiment,
  alpha = 0.05,
  correction = "auto",
  adjust = "bonferroni",
  digits = 3,
  p_digits = 3,
  effect_size = "pes",
  type = 3,
  followup_alpha = 0.10,
  simple_effects = FALSE,
  .on_result = NULL,
  .on_sphericity = NULL,
  .use_covariates = FALSE,
  .covariate_interactions = FALSE
) {

  # ----------------------------------------------------------
  # 補助関数
  # ----------------------------------------------------------

  clean_effect_name <- function(effect) {
    effect <- gsub("`", "", effect, fixed = TRUE)
    effect <- gsub("^\\(Intercept\\)$", "Intercept", effect)
    trimws(effect)
  }

  split_effect_factors <- function(effect) {
    effect <- clean_effect_name(effect)
    trimws(strsplit(effect, ":", fixed = TRUE)[[1]])
  }

  format_level_label <- function(factor_name, level) {
    text <- as.character(level)

    if (startsWith(text, factor_name)) {
      return(text)
    }

    if (grepl("^X?[0-9]+$", text)) {
      return(paste0(factor_name, sub("^X", "", text)))
    }

    paste0(factor_name, "=", text)
  }

  make_cell_label <- function(factors, levels) {
    labels <- mapply(
      format_level_label,
      factor_name = factors,
      level = levels,
      USE.NAMES = FALSE
    )

    separator <- if (any(grepl("=", labels, fixed = TRUE))) " | " else ""
    paste(labels, collapse = separator)
  }

  add_comparison_label <- function(
    data,
    effect_factors,
    compared_factor,
    fixed_factor = NULL
  ) {
    if (!"contrast" %in% names(data)) {
      return(data)
    }

    pairs <- strsplit(as.character(data$contrast), " - ", fixed = TRUE)

    data$comparison <- vapply(
      seq_len(nrow(data)),
      function(index) {
        pair <- pairs[[index]]
        if (length(pair) != 2) {
          return(as.character(data$contrast[index]))
        }

        left <- setNames(rep(NA_character_, length(effect_factors)), effect_factors)
        right <- left
        left[[compared_factor]] <- pair[1]
        right[[compared_factor]] <- pair[2]

        if (!is.null(fixed_factor)) {
          if (!fixed_factor %in% names(data)) {
            return(as.character(data$contrast[index]))
          }
          fixed_level <- as.character(data[[fixed_factor]][index])
          left[[fixed_factor]] <- fixed_level
          right[[fixed_factor]] <- fixed_level
        }

        if (anyNA(left) || anyNA(right)) {
          return(as.character(data$contrast[index]))
        }

        paste0(
          make_cell_label(effect_factors, unname(left[effect_factors])),
          " - ",
          make_cell_label(effect_factors, unname(right[effect_factors]))
        )
      },
      character(1)
    )

    data
  }

  quote_name <- function(x) {
    paste0("`", gsub("`", "", x, fixed = TRUE), "`")
  }

  find_column <- function(data, candidates, pattern = NULL) {
    exact <- intersect(candidates, names(data))
    if (length(exact) > 0) {
      return(exact[1])
    }

    if (!is.null(pattern)) {
      matched <- grep(pattern, names(data), value = TRUE, ignore.case = TRUE)
      if (length(matched) > 0) {
        return(matched[1])
      }
    }

    NA_character_
  }

  standardize_anova_table <- function(
    anova_table,
    correction_label,
    effect_size_name = "pes"
  ) {
    table_df <- as.data.frame(anova_table)
    table_df$Effect <- rownames(table_df)
    rownames(table_df) <- NULL

    effect_column <- find_column(table_df, c("Effect", "effect"))
    df1_column <- find_column(
      table_df,
      c("num Df", "num.Df", "Df", "df1"),
      "num.*df|^df1$"
    )
    df2_column <- find_column(
      table_df,
      c("den Df", "den.Df", "Df.res", "df2"),
      "den.*df|^df2$"
    )
    mse_column <- find_column(table_df, c("MSE", "mse"), "^mse$")
    f_column <- find_column(
      table_df,
      c("F", "F.value", "F value"),
      "^f($|[ ._])"
    )
    p_column <- find_column(
      table_df,
      c("Pr(>F)", "p.value", "p"),
      "pr\\(>f\\)|p.value|^p$"
    )
    es_column <- find_column(
      table_df,
      effect_size_name,
      paste0("^", effect_size_name, "$")
    )

    output <- data.frame(
      Effect = clean_effect_name(as.character(table_df[[effect_column]])),
      df1 = if (!is.na(df1_column)) table_df[[df1_column]] else NA_real_,
      df2 = if (!is.na(df2_column)) table_df[[df2_column]] else NA_real_,
      MSE = if (!is.na(mse_column)) table_df[[mse_column]] else NA_real_,
      F = if (!is.na(f_column)) table_df[[f_column]] else NA_real_,
      p.value = if (!is.na(p_column)) table_df[[p_column]] else NA_real_,
      stringsAsFactors = FALSE
    )

    effect_values <- if (!is.na(es_column)) {
      as.numeric(table_df[[es_column]])
    } else {
      rep(NA_real_, nrow(output))
    }

    # afexの出力形式が変わってpes列を取得できない場合に備え、
    # F値と自由度から偏イータ二乗を復元する。
    if (effect_size_name == "pes") {
      fallback <- with(
        output,
        (F * df1) / (F * df1 + df2)
      )
      replace_index <- is.na(effect_values) & is.finite(fallback)
      effect_values[replace_index] <- fallback[replace_index]
    }

    output[[effect_size_name]] <- effect_values
    output <- output[, c(
      "Effect", "df1", "df2", "MSE", "F",
      effect_size_name, "p.value"
    ), drop = FALSE]
    output$correction <- correction_label
    output
  }

  extract_sphericity_tables <- function(model) {
    if (is.null(model$Anova)) {
      return(list(mauchly = NULL, corrections = NULL, raw = NULL))
    }

    summary_object <- tryCatch(
      summary(model$Anova, multivariate = FALSE),
      error = function(e) NULL
    )

    if (is.null(summary_object)) {
      return(list(mauchly = NULL, corrections = NULL, raw = NULL))
    }

    mauchly <- NULL
    corrections <- NULL

    possible_mauchly_names <- c(
      "sphericity.tests",
      "sphericity.tests.table",
      "Mauchly's Test for Sphericity"
    )
    possible_correction_names <- c(
      "pval.adjustments",
      "sphericity.corrections",
      "Greenhouse-Geisser and Huynh-Feldt Corrections"
    )

    for (name in possible_mauchly_names) {
      if (!is.null(summary_object[[name]])) {
        mauchly <- as.data.frame(
          unclass(summary_object[[name]]),
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
        break
      }
    }

    for (name in possible_correction_names) {
      if (!is.null(summary_object[[name]])) {
        corrections <- as.data.frame(
          unclass(summary_object[[name]]),
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
        break
      }
    }

    if (!is.null(mauchly)) {
      mauchly$Effect <- clean_effect_name(rownames(mauchly))
      rownames(mauchly) <- NULL
      mauchly <- mauchly[, c("Effect", setdiff(names(mauchly), "Effect")), drop = FALSE]
    }

    if (!is.null(corrections)) {
      corrections$Effect <- clean_effect_name(rownames(corrections))
      rownames(corrections) <- NULL
      corrections <- corrections[, c("Effect", setdiff(names(corrections), "Effect")), drop = FALSE]
    }

    list(
      mauchly = mauchly,
      corrections = corrections,
      raw = summary_object
    )
  }

  run_levene_tests <- function(data, dv_name, within, between) {
    if (length(between) == 0) {
      return(NULL)
    }

    complete_rows <- stats::complete.cases(
      data[, unique(c(dv_name, within, between)), drop = FALSE]
    )
    test_data <- data[complete_rows, , drop = FALSE]

    if (nrow(test_data) == 0) {
      return(NULL)
    }

    run_one_test <- function(subset_data, label) {
      subset_group <- interaction(
        subset_data[, between, drop = FALSE],
        drop = TRUE,
        sep = ":"
      )

      if (nlevels(subset_group) < 2) {
        return(NULL)
      }

      result <- tryCatch(
        car::leveneTest(
          y = subset_data[[dv_name]],
          group = subset_group,
          center = stats::median
        ),
        error = function(e) NULL
      )

      if (is.null(result)) {
        return(NULL)
      }

      result_df <- as.data.frame(result)
      f_column <- find_column(result_df, c("F value", "F.value"), "^f")
      p_column <- find_column(result_df, c("Pr(>F)", "p.value"), "pr\\(>f\\)|p")
      df_column <- find_column(result_df, c("Df", "df"), "^df$")

      data.frame(
        Condition = label,
        df1 = if (!is.na(df_column)) result_df[[df_column]][1] else NA_real_,
        df2 = if (!is.na(df_column) && nrow(result_df) >= 2) result_df[[df_column]][2] else NA_real_,
        F = if (!is.na(f_column)) result_df[[f_column]][1] else NA_real_,
        p.value = if (!is.na(p_column)) result_df[[p_column]][1] else NA_real_,
        stringsAsFactors = FALSE
      )
    }

    if (length(within) == 0) {
      output <- run_one_test(test_data, "全体")
    } else {
      within_key <- interaction(
        test_data[, within, drop = FALSE],
        drop = TRUE,
        sep = ", "
      )
      split_data <- split(test_data, within_key)
      output_list <- lapply(names(split_data), function(label) {
        run_one_test(split_data[[label]], label)
      })
      output_list <- Filter(Negate(is.null), output_list)
      output <- if (length(output_list) > 0) do.call(rbind, output_list) else NULL
    }

    if (!is.null(output)) {
        rownames(output) <- NULL
    }

    output
  }

  summarize_emmeans_contrast <- function(
    contrast_object,
    effect_factors,
    compared_factor,
    fixed_factor = NULL
  ) {
    result_df <- as.data.frame(summary(contrast_object, infer = c(TRUE, TRUE), level = 1 - alpha))
    result_df <- add_comparison_label(
      data = result_df,
      effect_factors = effect_factors,
      compared_factor = compared_factor,
      fixed_factor = fixed_factor
    )
    result_df$contrast <- NULL


    result_df
  }

  run_main_effect_comparison <- function(model, factor_name) {
    emm <- emmeans::emmeans(
      model,
      specs = stats::as.formula(paste0("~ ", quote_name(factor_name))),
      at = if (.use_covariates) stats::setNames(as.list(rep(0, length(covariate))), covariate) else NULL
    )
    contrast <- emmeans::contrast(
      emm,
      method = "pairwise",
      adjust = "none"
    )

    list(
      emmeans = emm,
      contrast = contrast,
      table = summarize_emmeans_contrast(
        contrast,
        effect_factors = factor_name,
        compared_factor = factor_name
      )
    )
  }

  run_simple_comparison <- function(
    model,
    compared_factor,
    fixed_factor,
    effect_factors
  ) {
    formula_text <- paste0(
      "~ ", quote_name(compared_factor),
      " | ", quote_name(fixed_factor)
    )
    emm <- emmeans::emmeans(model, specs = stats::as.formula(formula_text),
      at = if (.use_covariates) stats::setNames(as.list(rep(0, length(covariate))), covariate) else NULL)
    contrast <- emmeans::contrast(
      emm,
      method = "pairwise",
      adjust = "none"
    )

    list(
      compared_factor = compared_factor,
      fixed_factor = fixed_factor,
      emmeans = emm,
      contrast = contrast,
      table = summarize_emmeans_contrast(
        contrast,
        effect_factors = effect_factors,
        compared_factor = compared_factor,
        fixed_factor = fixed_factor
      )
    )
  }

  run_follow_up <- function(model, effect, simple_factors = NULL) {
    factors <- split_effect_factors(effect)

    if (length(factors) == 1) {
      return(list(
        type = "main_effect",
        effect = effect,
        result = tryCatch(
          run_main_effect_comparison(model, factors[1]),
          error = function(e) list(error = conditionMessage(e))
        )
      ))
    }

    if (length(factors) == 2) {
      if (is.null(simple_factors)) {
        simple_factors <- factors
      }

      comparison_1 <- if (factors[1] %in% simple_factors) {
        tryCatch(
          run_simple_comparison(model, factors[1], factors[2], factors),
          error = function(e) list(error = conditionMessage(e))
        )
      } else {
        NULL
      }

      comparison_2 <- if (factors[2] %in% simple_factors) {
        tryCatch(
          run_simple_comparison(model, factors[2], factors[1], factors),
          error = function(e) list(error = conditionMessage(e))
        )
      } else {
        NULL
      }

      return(list(
        type = "two_way_interaction",
        effect = effect,
        factors = factors,
        simple_factors = simple_factors,
        first = comparison_1,
        second = comparison_2
      ))
    }

    list(
      type = "higher_order_interaction",
      effect = effect,
      factors = factors,
      message = paste0(
        length(factors),
        "要因以上の交互作用は自動展開せず、モデルを保存します。"
      )
    )
  }

  # ----------------------------------------------------------
  # パッケージと引数の確認
  # ----------------------------------------------------------

  required_packages <- c("afex", "emmeans", "car")
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_packages) > 0) {
    stop(
      "次のパッケージがインストールされていません: ",
      paste(missing_packages, collapse = ", "),
      "\ninstall.packages(c(",
      paste(sprintf('"%s"', missing_packages), collapse = ", "),
      ")) を実行してください。"
    )
  }

  required_elements <- c(
    "data", "id", "iv", "dv", "within", "between", "design"
  )
  missing_elements <- setdiff(required_elements, names(experiment))

  if (length(missing_elements) > 0) {
    stop(
      "experimentに次の要素がありません: ",
      paste(missing_elements, collapse = ", ")
    )
  }

  data <- experiment$data
  id <- experiment$id
  iv <- experiment$iv
  dv <- experiment$dv
  within <- experiment$within
  between <- experiment$between
  design <- experiment$design
  covariate <- if (!is.null(experiment$covariate)) experiment$covariate else character(0)

  if (!is.data.frame(data)) {
    stop("experiment$dataはデータフレームである必要があります。")
  }
  if (length(id) != 1) {
    stop("experiment$idには参加者ID列を1つ指定してください。")
  }
  if (length(iv) == 0) {
    stop("独立変数が指定されていません。")
  }
  if (length(dv) == 0) {
    stop("従属変数が指定されていません。")
  }
  if (!all(within %in% iv)) {
    stop("withinにはivに含まれる要因だけを指定してください。")
  }
  if (!all(between %in% iv)) {
    stop("betweenにはivに含まれる要因だけを指定してください。")
  }
  if (length(intersect(within, between)) > 0) {
    stop("同じ要因をwithinとbetweenの両方に指定することはできません。")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alphaには0より大きく1より小さい数値を指定してください。")
  }
  if (!is.numeric(followup_alpha) || length(followup_alpha) != 1 ||
      is.na(followup_alpha) || followup_alpha <= 0 || followup_alpha >= 1) {
    stop("followup_alphaには0より大きく1より小さい数値を指定してください。")
  }
  if (followup_alpha < alpha) {
    stop("followup_alphaにはalpha以上の値を指定してください。")
  }
  if (!is.logical(simple_effects) || length(simple_effects) != 1 || is.na(simple_effects)) {
    stop("simple_effectsにはTRUEまたはFALSEを指定してください。")
  }
  if (!correction %in% c("auto", "GG", "HF", "none")) {
    stop('correctionには"auto", "GG", "HF", "none"のいずれかを指定してください。')
  }
  if (!is.character(adjust) || length(adjust) != 1 || is.na(adjust)) {
    stop("adjustには多重比較補正法を文字列で指定してください。")
  }
  if (!is.numeric(digits) || length(digits) != 1 || is.na(digits) || digits < 0 || digits != floor(digits)) {
    stop("digitsには0以上の整数を指定してください。")
  }
  if (!is.numeric(p_digits) || length(p_digits) != 1 || is.na(p_digits) || p_digits < 1 || p_digits != floor(p_digits)) {
    stop("p_digitsには1以上の整数を指定してください。")
  }
  if (!effect_size %in% c("pes", "ges")) {
    stop('effect_sizeには"pes"または"ges"を指定してください。')
  }
  if (!type %in% c(2, 3)) {
    stop("typeには2または3を指定してください。")
  }

  if (!.use_covariates && length(covariate) > 0) {
    warning(
      "共変量「", paste(covariate, collapse = ", "),
      "」が指定されていますが、通常ANOVAでは調整しません。",
      "共変量を調整する場合はanalyze(method = \"ancova\")または\"permuco_ancova\"を使用してください。",
      call. = FALSE
    )
  }

  if (.use_covariates) {
    if (!is.character(covariate) || !length(covariate) || anyNA(covariate) ||
        any(!nzchar(covariate)) || anyDuplicated(c(id, iv, dv, covariate)) ||
        !setequal(iv, c(within, between))) {
      stop("ANCOVAには実験要因と重複しないexperiment$covariateの列名を指定してください。", call. = FALSE)
    }
    if (!is.logical(.covariate_interactions) || length(.covariate_interactions) != 1L ||
        is.na(.covariate_interactions)) {
      stop("covariate_interactionsにはTRUEまたはFALSEを指定してください。", call. = FALSE)
    }
    if (any(grepl("`", c(id, iv, dv, covariate), fixed = TRUE))) {
      stop("ANCOVAの列名にバッククォートを含めないでください。", call. = FALSE)
    }
  }
  required_columns <- unique(c(id, iv, dv, if (.use_covariates) covariate))
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0) {
    stop(
      "dataに次の列がありません: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  if (.use_covariates) {
    for (name in c(dv, covariate)) {
      if (!is.numeric(data[[name]]) || any(!is.finite(data[[name]]) & !is.na(data[[name]]))) {
        stop("ANCOVAの従属変数と共変量には有限の数値が必要です: ", name, call. = FALSE)
      }
    }
    for (name in c(between, covariate)) {
      counts <- vapply(split(data[[name]], data[[id]], drop = TRUE),
        function(x) length(unique(x[!is.na(x)])), integer(1))
      if (any(counts > 1L)) {
        stop("参加者ごとに一定の値を指定してください: ", name, call. = FALSE)
      }
    }
  }

  # afexでは要因をfactorとして扱う
  analysis_data <- data
  analysis_data[[id]] <- factor(analysis_data[[id]])
  for (factor_name in iv) {
    analysis_data[[factor_name]] <- factor(analysis_data[[factor_name]])
  }

  # Type III平方和に必要なsum-to-zero contrastはafex側で確認・設定される
  results <- list()

  # ----------------------------------------------------------
  # 従属変数ごとの分析
  # ----------------------------------------------------------

  for (target_name in dv) {
    current <- local({
      dv_name <- target_name
      complete_columns <- unique(c(id, iv, dv_name, if (.use_covariates) covariate))
      dv_data <- analysis_data[
        stats::complete.cases(analysis_data[, complete_columns, drop = FALSE]),
        ,
        drop = FALSE
      ]

      if (nrow(dv_data) == 0) {
        warning("従属変数「", dv_name, "」には分析可能な完全ケースがありません。")
        return(NULL)
      }

      centers <- numeric(0)
      excluded_ids <- character()
      if (.use_covariates) {
        # 不完全な反復測定の参加者は中心化前に除外する。
        if (length(within) && nrow(dv_data)) {
          cells <- interaction(dv_data[within], drop = FALSE)
          counts <- table(dv_data[[id]], cells)
          complete_ids <- rownames(counts)[rowSums(counts > 0L) == ncol(counts)]
          dv_data <- dv_data[as.character(dv_data[[id]]) %in% complete_ids, , drop = FALSE]
        }
        original_ids <- unique(as.character(analysis_data[[id]][!is.na(analysis_data[[id]])]))
        excluded_ids <- setdiff(original_ids, as.character(dv_data[[id]]))
        if (length(excluded_ids)) {
          warning("ANCOVA「", dv_name, "」の欠測・不完全測定により", length(excluded_ids),
            "名を除外しました。", call. = FALSE)
        }
        if (!nrow(dv_data)) stop("ANCOVAに使用できる完全な参加者がいません: ", dv_name, call. = FALSE)
        dv_data <- droplevels(dv_data)
        if (any(vapply(dv_data[iv], nlevels, integer(1)) < 2L)) {
          stop("ANCOVAの各実験要因には2水準以上が必要です。", call. = FALSE)
        }
        for (name in covariate) {
          values <- dv_data[[name]][!duplicated(dv_data[[id]])]
          if (length(unique(values)) < 2L) stop("共変量が一定値です: ", name, call. = FALSE)
          centers[name] <- mean(values)
          dv_data[[name]] <- dv_data[[name]] - centers[name]
        }
      }

      if (.use_covariates && .covariate_interactions && length(between)) {
        between_term <- paste(quote_name(between), collapse = " * ")
        fixed <- paste(c(between_term, paste0("(", between_term, ") * ", quote_name(covariate))), collapse = " + ")
        error <- if (length(within)) paste0(quote_name(id), " / (", paste(quote_name(within), collapse = " * "), ")") else quote_name(id)
        model_formula <- stats::as.formula(paste0(quote_name(dv_name), " ~ ", fixed, " + Error(", error, ")"))
        model <- afex::aov_car(model_formula, data = dv_data, type = type,
          observed = covariate, factorize = FALSE, return = "afex_aov")
      } else {
        model <- afex::aov_ez(
          id = id,
          dv = dv_name,
          data = dv_data,
          within = if (length(within) > 0) within else NULL,
          between = if (length(between) > 0) between else NULL,
          covariate = if (.use_covariates) covariate else NULL,
          observed = if (.use_covariates) covariate else NULL,
          type = type,
          factorize = FALSE,
          return = "afex_aov"
        )
      }

      requested_correction <- "none"

      # 偏イータ二乗（pes）は常に取得する。
      anova_numeric <- stats::anova(
        model,
        correction = requested_correction,
        es = "pes"
      )

      anova_table <- standardize_anova_table(
        anova_numeric,
        correction_label = requested_correction,
        effect_size_name = "pes"
      )

      # 必要に応じて一般化イータ二乗（ges）も追加する。
      anova_ges <- NULL
      if (effect_size == "ges") {
        anova_ges <- stats::anova(
          model,
          correction = requested_correction,
          es = "ges"
        )

        ges_table <- standardize_anova_table(
          anova_ges,
          correction_label = requested_correction,
          effect_size_name = "ges"
        )

        anova_table$ges <- ges_table$ges[
          match(anova_table$Effect, ges_table$Effect)
        ]

        anova_table <- anova_table[, c(
          "Effect", "df1", "df2", "MSE", "F", "pes", "ges",
          "p.value", "correction"
        ), drop = FALSE]
      }

      sphericity <- if (length(within) > 0) {
        extract_sphericity_tables(model)
      } else {
        list(mauchly = NULL, corrections = NULL, raw = NULL)
      }

      # 生の診断表を、共通の補正関数が受け取る効果別形式に整える。
      sphere_table <- anova_table[, "Effect", drop = FALSE]
      sphere_table$Within <- vapply(sphere_table$Effect, function(effect) {
        paste(intersect(split_effect_factors(effect), within), collapse = ":")
      }, character(1))
      sphere_table$GG <- sphere_table$HF <- sphere_table$W <- sphere_table$Mauchly.p <- NA_real_
      for (j in which(nzchar(sphere_table$Within))) {
        effect <- sphere_table$Effect[j]
        find_effect_row <- function(table) {
          if (is.null(table)) return(NA_integer_)
          candidates <- which(vapply(table$Effect, function(name) {
            setequal(setdiff(split_effect_factors(name), c("Intercept", "(Intercept)")), split_effect_factors(effect))
          }, logical(1)))
          if (length(candidates) == 1L) candidates else NA_integer_
        }
        k <- find_effect_row(sphericity$corrections)
        if (prod(vapply(dv_data[intersect(split_effect_factors(effect), within)], nlevels, integer(1)) - 1L) == 1) {
          sphere_table$GG[j] <- sphere_table$HF[j] <- 1
        } else if (!is.na(k)) {
          sphere_table$GG[j] <- sphericity$corrections[["GG eps"]][k]
          sphere_table$HF[j] <- min(1, sphericity$corrections[["HF eps"]][k])
        }
        k <- find_effect_row(sphericity$mauchly)
        if (!is.na(k)) {
          sphere_table$W[j] <- sphericity$mauchly[["Test statistic"]][k]
          sphere_table$Mauchly.p[j] <- sphericity$mauchly[["p-value"]][k]
        }
      }
      sphere_errors <- setNames(
          as.list(rep("球面性補正係数を取得できませんでした。",
            sum(nzchar(sphere_table$Within) & (!is.finite(sphere_table$GG) | !is.finite(sphere_table$HF))))),
          sphere_table$Effect[nzchar(sphere_table$Within) & (!is.finite(sphere_table$GG) | !is.finite(sphere_table$HF))])
      if (!is.null(.on_sphericity)) .on_sphericity(dv_name, sphere_table, sphere_errors)
      levene <- if (.use_covariates) NULL else run_levene_tests(
        data = dv_data,
        dv_name = dv_name,
        within = within,
        between = between
      )

      compute_followups <- function(effects, forced, main_effects) {
        contrasts <- list()
        comparisons <- list()
        for (effect in effects) {
          factors <- split_effect_factors(effect)
          if (.use_covariates && !all(factors %in% iv)) next
          compared <- if (effect %in% forced) intersect(factors, main_effects) else NULL
          follow_up <- run_follow_up(model, effect, compared)
          comparisons[[effect]] <- follow_up
          items <- switch(follow_up$type,
            main_effect = list(result = follow_up$result),
            two_way_interaction = Filter(Negate(is.null), list(first = follow_up$first, second = follow_up$second)),
            list(list(error = follow_up$message)))
          for (index in seq_along(items)) {
            item <- items[[index]]
            name <- if (length(items) == 1L) effect else paste0(effect, "/", item$compared_factor)
            if (!is.null(item$error)) {
              contrasts[[name]] <- list(error = item$error, source_effect = effect, source_slot = names(items)[index])
              next
            }
            contrasts[[name]] <- list(
              raw = item$contrast, comparison_grid = item$contrast,
              unadjusted = item$table, formatted = item$table,
              compared_factor = item$compared_factor, source_effect = effect,
              source_slot = names(items)[index], error = NULL)
          }
        }
        attr(contrasts, "comparisons") <- comparisons
        contrasts
      }

      result <- list(
        analysis = if (.use_covariates) "ancova" else "anova",
        dependent_variable = dv_name,
        design = design,
        within = within,
        between = between,
        covariate = covariate,
        covariate_adjusted = .use_covariates,
        covariate_centers = centers,
        covariate_interactions = .use_covariates && .covariate_interactions,
        within_covariate_interactions = .use_covariates && length(within) > 0L,
        followup_factors = iv,
        excluded_ids = excluded_ids,
        analysis_data = dv_data,
        model = model,
        anova = anova_numeric,
        anova_ges = anova_ges,
        anova_table = anova_table,
        sphericity_diagnostics = sphericity,
        sphericity_table = sphere_table,
        sphericity_errors = sphere_errors,
        levene = levene,
        compute_followups = compute_followups,
        contrasts = list(),
        anova_table_uncorrected = anova_table,
        anova_type = as.character(type),
        formula = stats::formula(model$lm),
        alpha = alpha,
        followup_alpha = followup_alpha,
        simple_effects = simple_effects,
        correction = correction,
        applied_within_correction = requested_correction,
        adjust = "none",
        effect_size = effect_size,
        primary_effect_size = "pes",
        effect_size_label = "partial_eta_squared",
        effect_size_basis = "original_scale",
        type = type,
        digits = digits,
        p_digits = p_digits
      )

      class(result) <- c(
        if (.use_covariates) "analysis_ancova_result" else "analysis_anova_result",
        "analysis_test_result",
        "list"
      )

      result
    })
    if (!is.null(current)) {
      if (!is.null(.on_result)) current <- .on_result(current)
      results[[target_name]] <- current
    }
  }

  attr(results, "stage") <- "test"
  attr(results, "type") <- if (.use_covariates) "ancova" else "anova"
  attr(results, "experiment") <- experiment
  class(results) <- c(if (.use_covariates) "analysis_ancova" else "analysis_anova", "analysis_test", "list")
  invisible(results)
}
