# ============================================================
# ART / ART-C
# 統計計算は ARTool、car、stats、emmeans に委譲する。
# 未補正の検定表・ε・下位検定用関数を返す。制御・補正・表示はhelperが担当。
# GG/HFは効果別ART順位データから求める拡張。ARTool組込みの補正ではない。
# 独自コードはデータの取得・整形、モデル式の生成、標準関数の呼び出し、
# 比較結果からの行抽出、および表示を担当する。
#
# 注意：ARTool の art.con() は Error(ID / (A * B)) のような
# 入れ子の反復測定誤差項と交互作用の組み合わせで失敗することがある。
# 反復測定の複数要因比較では、公式 ART-C の順位データを取得し、
# 元の要因表現で同じ固定効果空間・同じ Error() を aov に指定する。
# この互換性対処は実データでの実行検証が必要。検定量は自作しない。
# ============================================================

t_art <- function(
    experiment,
    alpha = 0.05,
    adjust = "bonferroni",
    digits = 3,
    p_digits = 3,
    followup_alpha = 0.10,
    simple_effects = FALSE,
    type = NULL,
    sphericity = c("none", "GG", "HF"),
    .on_result = NULL,
    .on_sphericity = NULL
) {
  if (!requireNamespace("ARTool", quietly = TRUE)) {
    stop('ARTool が必要です。install.packages("ARTool") を実行してください。', call. = FALSE)
  }
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    stop('emmeans が必要です。install.packages("emmeans") を実行してください。', call. = FALSE)
  }

  sphericity <- match.arg(sphericity)

  needed <- c("data", "id", "iv", "dv", "within", "between", "design")
  missing <- setdiff(needed, names(experiment))
  if (length(missing)) {
    stop("experiment に必要な要素がありません: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!is.data.frame(experiment$data) || length(experiment$id) != 1L ||
      !length(experiment$iv) || !length(experiment$dv)) {
    stop("experiment の data / id / iv / dv を確認してください。", call. = FALSE)
  }
  if (!setequal(experiment$iv, c(experiment$within, experiment$between)) ||
      length(intersect(experiment$within, experiment$between))) {
    stop("experiment$within と $between が $iv に一致しません。", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha <= 0 || alpha >= 1 ||
      !is.numeric(followup_alpha) || length(followup_alpha) != 1L ||
      is.na(followup_alpha) || followup_alpha < alpha || followup_alpha >= 1) {
    stop("alpha および followup_alpha は 0 < alpha <= followup_alpha < 1 にしてください。", call. = FALSE)
  }
  if (!is.logical(simple_effects) || length(simple_effects) != 1L ||
      is.na(simple_effects)) {
    stop("simple_effects には TRUE または FALSE を指定してください。", call. = FALSE)
  }
  if (!is.character(adjust) || length(adjust) != 1L || is.na(adjust)) {
    stop("adjust には emmeans の多重比較補正法を指定してください。", call. = FALSE)
  }
  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) ||
      digits < 0 || digits != floor(digits) ||
      !is.numeric(p_digits) || length(p_digits) != 1L || is.na(p_digits) ||
      p_digits < 1 || p_digits != floor(p_digits)) {
    stop("digits と p_digits は正しい桁数を指定してください。", call. = FALSE)
  }

  if (length(experiment$within) && !requireNamespace("car", quietly = TRUE)) {
    stop('球面性の推定には car が必要です。install.packages("car") を実行してください。',
         call. = FALSE)
  }

  within <- experiment$within
  between <- experiment$between
  iv <- experiment$iv
  id <- experiment$id
  dv <- experiment$dv
  covariate <- if (is.null(experiment$covariate)) character(0) else experiment$covariate

  if (length(covariate)) {
    warning("ART は共変量を調整しません: ", paste(covariate, collapse = ", "), call. = FALSE)
  }

  # ARTool の Error() モデルは Type I に制限される。
  anova_type <- if (length(within)) "I" else "III"
  if (!is.null(type)) {
    choices <- c("1" = "I", "2" = "II", "3" = "III", "I" = "I", "II" = "II", "III" = "III")
    requested <- unname(choices[as.character(type)])
    if (length(requested) != 1L || is.na(requested)) {
      stop("type は 1、2、3、I、II、III のいずれかです。", call. = FALSE)
    }
    if (length(within) && requested != "I") {
      stop("Error() を使う反復測定 ART では Type I のみ利用できます。", call. = FALSE)
    }
    anova_type <- requested
  }

  # ------------------------------------------------------------
  # 以下の補助関数は表示用の整形とエラー記録を担当する。
  # ------------------------------------------------------------
  # 遅延評価の式を同じ条件で実行し、失敗した処理の名前だけを記録する。
  contrast_stage <- NULL
  at_stage <- function(label, expr) {
    contrast_stage <<- label
    force(expr)
  }
  call_text <- function(x) {
    if (is.null(x)) return(NULL)
    paste(deparse(x, width.cutoff = 120L), collapse = " ")
  }
  runtime_info <- list(
    R = R.version.string,
    ARTool = as.character(utils::packageVersion("ARTool")),
    emmeans = as.character(utils::packageVersion("emmeans"))
  )

  split_effect <- function(effect) trimws(strsplit(gsub("`", "", effect, fixed = TRUE), ":", fixed = TRUE)[[1]])

  column <- function(df, options) {
    name <- intersect(options, names(df))
    if (length(name)) df[[name[1]]] else NULL
  }

  anova_display <- function(a) {
    df <- as.data.frame(a, check.names = FALSE, stringsAsFactors = FALSE)
    effects <- column(df, c("Term", "Effect", "effect"))
    if (is.null(effects)) effects <- rownames(df)
    error <- column(df, c("Error", "error"))
    df1 <- column(df, c("Df", "df1", "num Df"))
    df2 <- column(df, c("Df.res", "df2", "den Df"))
    f <- column(df, c("F value", "F", "F.value"))
    p <- column(df, c("Pr(>F)", "p.value", "p"))
    if (is.null(f) || is.null(p)) {
      stop("ARTool の ANOVA 表から F / p 値を取得できませんでした。", call. = FALSE)
    }

    # 偏 eta^2 は ARTool 公式の効果量説明に掲載された平方和の式を使用。
    # 平方和がなければ独自の近似で補完せず、NA とする。
    ss <- column(df, c("Sum Sq", "Sum.Sq"))
    ss_res <- column(df, c("Sum Sq.res", "Sum.Sq.res"))
    pes <- rep(NA_real_, nrow(df))
    if (!is.null(ss) && !is.null(ss_res)) {
      den <- ss + ss_res
      valid <- is.finite(den) & den > 0
      pes[valid] <- ss[valid] / den[valid]
    }

    data.frame(
      Effect = gsub("`", "", as.character(effects), fixed = TRUE),
      Error = if (is.null(error)) rep(NA_character_, nrow(df)) else as.character(error),
      df1 = if (is.null(df1)) rep(NA_real_, nrow(df)) else as.numeric(df1),
      df2 = if (is.null(df2)) rep(NA_real_, nrow(df)) else as.numeric(df2),
      F = as.numeric(f), pes = pes, p.value = as.numeric(p),
      stringsAsFactors = FALSE
    )
  }

  # 比較名の例："A1,B1 - A2,B1"。ARTool の出力ラベルを読むだけ。
  read_levels <- function(label, n_factors) {
    parts <- trimws(strsplit(label, ",", fixed = TRUE)[[1]])
    if (length(parts) != n_factors) return(rep(NA_character_, n_factors))
    parts
  }

  # 数値の比較・検定は実施しない。ARTool の比較表から
  # 「1要因だけ違う」行を選び、emmeans の emmGrid を部分抽出する。
  select_simple_pairs <- function(contrast_grid, factors, compared_factors) {
    reference <- at_stage("比較表の評価：summary(adjust = none)",
      as.data.frame(summary(contrast_grid, adjust = "none")))
    if (!"contrast" %in% names(reference)) {
      stop("ART-C の contrast 列を取得できません。", call. = FALSE)
    }
    contrast_stage <<- "単純比較ラベルの解析・抽出"
    parsed <- lapply(as.character(reference$contrast), function(label) {
      halves <- strsplit(label, " - ", fixed = TRUE)[[1]]
      if (length(halves) != 2L) return(NULL)
      left <- read_levels(halves[1], length(factors))
      right <- read_levels(halves[2], length(factors))
      if (anyNA(left) || anyNA(right)) return(NULL)
      changed <- which(left != right)
      if (length(changed) != 1L || !(factors[changed] %in% compared_factors)) return(NULL)
      fixed <- setdiff(seq_along(factors), changed)
      list(
        comparison = paste0(paste(left, collapse = ""), " - ", paste(right, collapse = "")),
        compared_factor = factors[changed],
        fixed_conditions = paste(paste0(factors[fixed], "=", left[fixed]), collapse = ", "),
        fixed_factor = paste(factors[fixed], collapse = ":"),
        fixed_level = paste(left[fixed], collapse = ",")
      )
    })
    rows <- which(!vapply(parsed, is.null, logical(1)))
    if (!length(rows)) {
      stop("ART-C から単純効果の比較行を抽出できませんでした。", call. = FALSE)
    }
    # ここで行だけを絞る。補正対象の比較群はこの行集合としてhelperへ渡す。
    metadata <- do.call(rbind, lapply(parsed[rows], function(x) {
      as.data.frame(x, stringsAsFactors = FALSE)
    }))
    list(grid = at_stage("単純比較の行抽出：emmGrid[rows]", contrast_grid[rows]),
         metadata = metadata)
  }

  make_contrast_table <- function(grid, factors, compared_factors) {
    if (length(factors) == 1L) {
      df <- at_stage("主効果比較の未補正集計：summary",
        as.data.frame(summary(grid, adjust = "none", infer = c(TRUE, TRUE), level = 1 - alpha)))
      df$comparison <- as.character(df$contrast)
      df$contrast <- NULL
    } else {
      selection <- select_simple_pairs(grid, factors, compared_factors)
      df <- at_stage("単純比較の未補正集計：summary",
        as.data.frame(summary(selection$grid, adjust = "none", infer = c(TRUE, TRUE), level = 1 - alpha)))
      df$comparison <- selection$metadata$comparison
      df$compared_factor <- selection$metadata$compared_factor
      df$fixed_conditions <- selection$metadata$fixed_conditions
      df$fixed_factor <- selection$metadata$fixed_factor
      df$fixed_level <- selection$metadata$fixed_level
      df$contrast <- NULL
    }
    contrast_stage <<- "比較結果の表示列の整形"
    rownames(df) <- NULL
    list(table = df, grid = if (length(factors) == 1L) grid else selection$grid)
  }

  # 効果ごとの ART 順位データを、参加者 x 被験者内条件の行列へ整形する。
  # εとMauchly検定は car の公開APIに委譲し、Type I のF値は置き換えない。
  prepare_sphericity_layout <- function(dat) {
    ids <- unique(as.character(dat[[id]]))
    id_index <- match(as.character(dat[[id]]), ids)
    level_sets <- lapply(dat[within], levels)
    sizes <- lengths(level_sets)
    if (any(sizes < 2L)) stop("被験者内要因は2水準以上必要です。", call. = FALSE)
    idata <- expand.grid(level_sets, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    for (name in within) idata[[name]] <- factor(idata[[name]], levels = level_sets[[name]])
    # expand.gridと同じ順序（先頭要因が最も速く変化）で条件番号を生成する。
    condition_index <- rep(1L, nrow(dat))
    stride <- 1L
    for (name in within) {
      condition_index <- condition_index + (as.integer(dat[[name]]) - 1L) * stride
      stride <- stride * nlevels(dat[[name]])
    }
    if (anyNA(id_index) || anyNA(condition_index)) {
      stop("IDまたは被験者内条件に欠損があります。", call. = FALSE)
    }
    cell_index <- id_index + (condition_index - 1L) * length(ids)
    cell_counts <- matrix(tabulate(cell_index, nbins = length(ids) * nrow(idata)), nrow = length(ids))
    if (any(cell_counts == 0L)) {
      stop("被験者内条件に欠測があります。完全な反復測定データが必要です。", call. = FALSE)
    }
    if (length(unique(as.integer(cell_counts))) != 1L) {
      stop("球面性評価のセル内反復は、全参加者×条件で回数が等しい必要があります。", call. = FALSE)
    }
    first <- match(ids, as.character(dat[[id]]))
    subject_data <- dat[first, between, drop = FALSE]
    for (name in between) {
      values <- as.character(subject_data[[name]])[id_index]
      if (!identical(as.character(dat[[name]]), values)) {
        stop("被験者間要因がID内で変化しています: ", name, call. = FALSE)
      }
    }
    rownames(subject_data) <- NULL
    list(ids = ids, id_index = id_index, condition_index = condition_index,
         idata = idata, subject_data = subject_data, cell_index = cell_index, cell_counts = cell_counts)
  }

  effect_epsilon <- function(model, effect, layout) {
    within_effect <- intersect(split_effect(effect), within)
    q <- prod(vapply(model$data[within_effect], nlevels, integer(1)) - 1L)
    if (q == 1L) {
      return(list(GG = 1, HF = 1, HF_raw = 1, W = NA_real_,
                  Mauchly.p = NA_real_, status = "1 df: epsilon = 1"))
    }
    ranks <- model$aligned.ranks[[effect]]
    if (is.null(ranks) || length(ranks) != nrow(model$data) || any(!is.finite(ranks))) {
      stop("効果に対応するART順位データを取得できません: ", effect, call. = FALSE)
    }
    responses <- matrix(NA_real_, nrow = length(layout$ids), ncol = nrow(layout$idata))
    # analyze()は同一ID×条件の元データを事前に平均してから渡す。
    rank_sums <- rowsum(matrix(ranks, ncol = 1L), layout$cell_index, reorder = FALSE)
    positions <- as.integer(rownames(rank_sums))
    responses[positions] <- rank_sums[, 1L] / layout$cell_counts[positions]
    if (anyNA(responses)) stop("順位データの条件配置に欠損があります。", call. = FALSE)
    colnames(responses) <- paste0("cell", seq_len(ncol(responses)))
    # 数式内の応答行列を公開APIのlmへ渡す。betweenの完全要因配置で残差を求める。
    rhs <- if (length(between)) paste(sprintf("`%s`", between), collapse = " * ") else "1"
    mlm_formula <- stats::as.formula(paste("responses ~", rhs), env = environment())
    codes <- if (length(between)) stats::setNames(rep(list("contr.sum"), length(between)), between) else NULL
    fit <- stats::lm(mlm_formula, data = layout$subject_data, contrasts = codes)
    if (fit$df.residual <= 0L) stop("球面性推定の残差自由度がありません。", call. = FALSE)
    idesign <- stats::as.formula(paste("~", paste(sprintf("`%s`", within), collapse = " * ")))
    # ここではcarのType IIIのF・pは採用せず、残差に基づくεだけを取得する。
    repeated <- car::Anova(fit, idata = layout$idata, idesign = idesign,
                           type = "III", icontrasts = c("contr.sum", "contr.poly"))
    details <- summary(repeated, multivariate = FALSE, univariate = TRUE)
    adjustments <- details$pval.adjustments
    candidates <- which(vapply(rownames(adjustments), function(term) {
      parts <- setdiff(split_effect(term), "(Intercept)")
      setequal(parts, within_effect)
    }, logical(1)))
    if (length(candidates) != 1L) {
      stop("対象の被験者内効果に対応するε行が見つかりません: ", effect, call. = FALSE)
    }
    j <- candidates[1L]
    gg <- as.numeric(adjustments[j, "GG eps"])
    hf_raw <- as.numeric(adjustments[j, "HF eps"])
    if (!is.finite(gg) || !is.finite(hf_raw) || gg < 1/q - 1e-8 ||
        gg > 1 + 1e-8 || hf_raw <= 0) {
      stop("有効なGG/HFのεを推定できません: ", effect, call. = FALSE)
    }
    gg <- min(1, max(1/q, gg)) # 許容した浮動小数点誤差だけを境界へ戻す。
    hf <- min(1, hf_raw)       # car標準処理と同様、HFが1を超えた場合は1を使用。
    tests <- details$sphericity.tests
    test_row <- match(rownames(adjustments)[j], rownames(tests))
    w <- mp <- NA_real_
    if (!is.na(test_row)) {
      w <- as.numeric(tests[test_row, "Test statistic"])
      mp <- as.numeric(tests[test_row, "p-value"])
    }
    list(GG = gg, HF = hf, HF_raw = hf_raw, W = w, Mauchly.p = mp, status = "ok")
  }

  sphericity_results <- function(model, original) {
    output <- original
    output$GG <- output$HF <- NA_real_
    output$W <- output$Mauchly.p <- output$HF_raw <- NA_real_
    output$Within <- vapply(original$Effect, function(e) {
      paste(intersect(split_effect(e), within), collapse = ":")
    }, character(1))
    output$status <- ifelse(nzchar(output$Within), "pending", "between: not applicable")
    errors <- list()
    if (!length(within)) return(list(table = output, errors = errors))
    layout <- tryCatch(prepare_sphericity_layout(model$data), error = function(e) e)
    for (j in which(nzchar(output$Within))) {
      effect <- original$Effect[j]
      epsilon <- tryCatch({
        if (inherits(layout, "error")) stop(layout)
        effect_epsilon(model, effect, layout)
      }, error = function(e) e)
      if (inherits(epsilon, "error")) {
        errors[[effect]] <- conditionMessage(epsilon)
        output$status[j] <- conditionMessage(epsilon)
        next
      }
      for (name in c("GG", "HF", "HF_raw", "W", "Mauchly.p", "status")) {
        output[[name]][j] <- epsilon[[name]]
      }

    }
    list(table = output, errors = errors)
  }

  # ART-C の整列・順位変換は artlm.con() に委譲する。
  # 結合要因を元の完全要因配置へ再パラメータ化し、Error() は維持する。
  repeated_artc <- function(model, effect, factors) {
    contrast_stage <<- "ART-C 用順位モデルの生成：artlm.con"
    artc_model <- ARTool::artlm.con(model, effect)

    contrast_stage <<- "ART-C 用順位データの取得"
    model_call <- attr(artc_model, "call")
    # aovlist は formula() で直接扱わず、aov が保存した terms を読む。
    model_terms <- attr(artc_model, "terms", exact = TRUE)
    model_env <- attr(model_terms, ".Environment", exact = TRUE)
    if (is.null(model_call$data) || !is.environment(model_env)) {
      stop("ART-C モデルから順位データの参照を取得できません。", call. = FALSE)
    }
    # ARTool がモデル式の環境に保持した、実際に推定へ渡したデータを読む。
    ranked_data <- eval(model_call$data, envir = model_env)
    if (!is.data.frame(ranked_data) || !".y" %in% names(ranked_data) ||
        !is.numeric(ranked_data$.y) || nrow(ranked_data) != nrow(model$data) ||
        !all(c(id, iv) %in% names(ranked_data))) {
      stop("ART-C の順位データの構造を確認できません。", call. = FALSE)
    }
    # 行の並びと実験条件が同一であることを確認する（値は変更しない）。
    for (name in c(id, iv)) {
      if (!identical(ranked_data[[name]], model$data[[name]])) {
        stop("ART-C の順位データの行・条件が元データと一致しません: ", name,
             call. = FALSE)
      }
    }

    contrast_stage <<- "元の要因表現・同じ誤差層で推定：aov"
    refit_formula <- stats::update.formula(model$formula, .y ~ .)
    # emmeans が aov の call を再評価できるよう、データと式を保持する環境を指定。
    environment(refit_formula) <- environment()
    contrast_codes <- stats::setNames(rep(list("contr.sum"), length(iv)), iv)
    refit <- stats::aov(refit_formula, data = ranked_data, contrasts = contrast_codes)

    contrast_stage <<- "元の要因表現で周辺平均を生成：emmeans"
    marginal <- emmeans::emmeans(refit, specs = factors)
    contrast_stage <<- "組合せ水準の標準ペア比較：contrast"
    # 比較名の区切りだけを従来の ART-C 出力へ揃え、終了時に復元する。
    old_emmeans <- getOption("emmeans")
    on.exit(options(emmeans = old_emmeans), add = TRUE)
    emmeans::emm_options(sep = ",")
    grid <- emmeans::contrast(marginal, method = "pairwise", adjust = "none")
    list(grid = grid, rank_model = artc_model, comparison_model = refit)
  }

  # ------------------------------------------------------------
  # experiment の整形（値は変更せず、型だけ揃える）
  # ------------------------------------------------------------
  data <- experiment$data
  required_cols <- unique(c(id, iv, dv, covariate))
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols)) {
    stop("次のデータ列がありません: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  data[[id]] <- factor(data[[id]])
  data[iv] <- lapply(data[iv], factor)
  fixed <- paste(sprintf("`%s`", iv), collapse = " * ")
  error_term <- if (length(within)) {
    paste0("Error(`", id, "` / (", paste(sprintf("`%s`", within), collapse = " * "), "))")
  } else NULL

  results <- list()
  for (target_name in dv) {
    results[[target_name]] <- local({
      target <- target_name
      formula_text <- paste0("`", target, "` ~ ", fixed,
                             if (is.null(error_term)) "" else paste0(" + ", error_term))
      model_formula <- stats::as.formula(formula_text)

      model <- ARTool::art(model_formula, data = data)
      art_anova <- stats::anova(model, type = anova_type)
      table_uncorrected <- anova_display(art_anova)
      sphere <- sphericity_results(model, table_uncorrected)
      if (!is.null(.on_sphericity)) .on_sphericity(target, sphere$table, sphere$errors)
      diagnostics <- summary(model)
      compute_followups <- function(followups, forced, main_effects) {
        contrasts <- list()
        for (effect in followups) {
          factors <- split_effect(effect)
          compared <- if (effect %in% forced) intersect(factors, main_effects) else factors
          grid <- NULL
          rank_model <- comparison_model <- NULL
          method_used <- if (length(within) && length(factors) > 1L) {
            "ART-C ranks + original-factor aov"
          } else "ARTool::art.con"
          error_details <- NULL
          warning_details <- list()
          contrast_stage <<- "ARTool::art.con：モデル生成・比較"
          result <- tryCatch(
            withCallingHandlers({
              if (length(within) && length(factors) > 1L) {
                rebuilt <- repeated_artc(model, effect, factors)
                grid <- rebuilt$grid
                rank_model <- rebuilt$rank_model
                comparison_model <- rebuilt$comparison_model
              } else {
                grid <- ARTool::art.con(model, effect, method = "pairwise", adjust = "none")
              }
              contrast_stage <<- "比較結果の整形"
              prepared <- make_contrast_table(grid, factors, compared)
              list(raw = grid, comparison_grid = prepared$grid,
                   unadjusted = prepared$table, formatted = prepared$table, error = NULL)
            }, warning = function(w) {
              warning_details[[length(warning_details) + 1L]] <<- list(
                stage = contrast_stage, message = conditionMessage(w),
                call = call_text(conditionCall(w))
              )
              # 標準の警告表示は抑制しない。
            }, error = function(e) {
              # tryCatch でスタックが巻き戻される前に取得する。
              error_details <<- list(
                stage = contrast_stage,
                message = conditionMessage(e),
                condition_class = class(e),
                call = call_text(conditionCall(e)),
                calls = vapply(sys.calls(), call_text, character(1)),
                effect = effect, formula = formula_text,
                runtime = runtime_info
              )
            }), error = function(e) {
              list(raw = grid, formatted = NULL, error = conditionMessage(e))
            }
          )
          result$method <- method_used
          result$rank_model <- rank_model
          result$comparison_model <- comparison_model
          result$error_details <- error_details
          result$warnings <- warning_details
          contrasts[[effect]] <- result
        }
        contrasts
      }

      result <- list(
        analysis = "art", dependent_variable = target,
        design = experiment$design, within = within, between = between,
        covariate = covariate, covariate_adjusted = FALSE,
        analysis_data = data[, unique(c(id, iv, covariate, target)), drop = FALSE],
        formula = model_formula, model = model, diagnostics = diagnostics,
        anova = art_anova, anova_table = table_uncorrected,
        anova_table_uncorrected = table_uncorrected,
        sphericity = sphericity, sphericity_table = sphere$table,
        sphericity_errors = sphere$errors,
        sphericity_basis = "effect_specific_aligned_rank_residuals",
        sphericity_inference = "GG/HF extension; not an ARTool built-in option",
        effect_size = "pes", effect_size_basis = "aligned_rank_transformed",
        anova_type = anova_type, model_structure = if (length(within)) "repeated_measures_aov" else "factorial_lm",
        error_term = error_term, runtime_info = runtime_info,
        compute_followups = compute_followups, contrasts = list(),
        alpha = alpha, followup_alpha = followup_alpha,
        simple_effects = simple_effects, adjust = adjust,
        digits = digits, p_digits = p_digits
      )
      class(result) <- c("analysis_art_result", "analysis_test_result", "list")
      result
    })
    if (!is.null(.on_result)) results[[target_name]] <- .on_result(results[[target_name]])
  }
  attr(results, "stage") <- "test"
  attr(results, "type") <- "art"
  attr(results, "experiment") <- experiment
  class(results) <- c("analysis_art", "analysis_test", "list")
  invisible(results)
}
