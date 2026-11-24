# ============================================================
# 分析前処理：記述統計
#
# 公開関数：
#   p_describe()
#
# p_prepare()が返すanalysis_experimentを受け取り、従属変数ごとに
# 全体および要因の全組合せ別の記述統計を作成します。
#
# 各従属変数の結果には次のデータを保存します。
#   overall   : 全体の記述統計
#   tables    : 全体・各要因・全交互作用の記述統計表
#   cells     : 全要因を組み合わせた条件セル別の記述統計
#   raw_data  : 元の観測値をlong形式に整えたデータ
# ============================================================


# ------------------------------------------------------------
# 内部関数：表示用の丸め
# ------------------------------------------------------------

.p_round_descriptive_table <- function(data, digits) {
  output <- data
  integer_columns <- c("n", "observations", "missing")

  for (name in names(output)) {
    if (is.numeric(output[[name]]) && !name %in% integer_columns) {
      output[[name]] <- round(output[[name]], digits)
    }
  }

  output
}


# ------------------------------------------------------------
# 内部関数：数値ベクトルの要約
# ------------------------------------------------------------

.p_summarize_values <- function(values, observations, missing, ci) {
  x <- values[!is.na(values)]
  n <- length(x)

  if (n == 0) {
    return(data.frame(
      n = 0L,
      observations = as.integer(observations),
      missing = as.integer(missing),
      mean = NA_real_,
      sd = NA_real_,
      se = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      median = NA_real_,
      q1 = NA_real_,
      q3 = NA_real_,
      iqr = NA_real_,
      min = NA_real_,
      max = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  mean_value <- mean(x)
  sd_value <- if (n >= 2) stats::sd(x) else NA_real_
  se_value <- if (n >= 2) sd_value / sqrt(n) else NA_real_

  critical_value <- if (n >= 2) {
    stats::qt(1 - (1 - ci) / 2, df = n - 1)
  } else {
    NA_real_
  }

  margin <- if (!is.na(se_value)) {
    critical_value * se_value
  } else {
    NA_real_
  }

  quartiles <- stats::quantile(
    x,
    probs = c(0.25, 0.75),
    names = FALSE,
    type = 7
  )

  data.frame(
    n = as.integer(n),
    observations = as.integer(observations),
    missing = as.integer(missing),
    mean = mean_value,
    sd = sd_value,
    se = se_value,
    ci_lower = mean_value - margin,
    ci_upper = mean_value + margin,
    median = stats::median(x),
    q1 = quartiles[1],
    q3 = quartiles[2],
    iqr = stats::IQR(x, type = 7),
    min = min(x),
    max = max(x),
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# 内部関数：複数列からグループキーを作成
# ------------------------------------------------------------

.p_interaction_key <- function(data, variables) {
  if (length(variables) == 0) {
    return(factor(rep("all", nrow(data))))
  }

  do.call(
    interaction,
    c(
      data[variables],
      list(
        drop = TRUE,
        lex.order = TRUE,
        sep = "\r"
      )
    )
  )
}


# ------------------------------------------------------------
# 内部関数：参加者 × 条件セルへ集約
# ------------------------------------------------------------

.p_participant_data <- function(data, id, dv, factors) {
  variables <- c(id, factors)
  selected <- data[, c(variables, dv), drop = FALSE]
  groups <- .p_interaction_key(selected, variables)
  rows <- split(seq_len(nrow(selected)), groups, drop = TRUE)

  output <- lapply(
    rows,
    function(index) {
      values <- selected[[dv]][index]
      identifiers <- selected[index[1], variables, drop = FALSE]

      cbind(
        identifiers,
        data.frame(
          .value = if (all(is.na(values))) {
            NA_real_
          } else {
            mean(values, na.rm = TRUE)
          },
          .observations = as.integer(sum(!is.na(values))),
          .missing = as.integer(sum(is.na(values))),
          stringsAsFactors = FALSE
        )
      )
    }
  )

  if (length(output) == 0) {
    empty <- selected[FALSE, variables, drop = FALSE]
    empty$.value <- numeric(0)
    empty$.observations <- integer(0)
    empty$.missing <- integer(0)
    return(empty)
  }

  output <- do.call(rbind, output)
  rownames(output) <- NULL
  output
}


# ------------------------------------------------------------
# 内部関数：指定要因別の記述統計
# ------------------------------------------------------------

.p_summarize_by <- function(data, id, dv, factors, ci) {
  person <- .p_participant_data(data, id, dv, factors)

  if (length(factors) == 0) {
    return(
      .p_summarize_values(
        values = person$.value,
        observations = sum(person$.observations),
        missing = sum(person$.missing),
        ci = ci
      )
    )
  }

  groups <- .p_interaction_key(person, factors)
  rows <- split(seq_len(nrow(person)), groups, drop = TRUE)

  output <- lapply(
    rows,
    function(index) {
      factor_values <- person[index[1], factors, drop = FALSE]

      cbind(
        factor_values,
        .p_summarize_values(
          values = person$.value[index],
          observations = sum(person$.observations[index]),
          missing = sum(person$.missing[index]),
          ci = ci
        )
      )
    }
  )

  output <- do.call(rbind, output)
  rownames(output) <- NULL
  output
}


# ------------------------------------------------------------
# 内部関数：要因の全組合せを作成
# ------------------------------------------------------------

.p_factor_sets <- function(iv) {
  output <- list()

  for (size in seq_along(iv)) {
    output <- c(
      output,
      utils::combn(iv, size, simplify = FALSE)
    )
  }

  names(output) <- vapply(
    output,
    paste,
    collapse = ":",
    FUN.VALUE = character(1)
  )

  output
}


# ------------------------------------------------------------
# 内部関数：入力オブジェクトと表示引数の検証
# ------------------------------------------------------------

.p_validate_describe_input <- function(experiment, ci, digits, show) {
  required <- c(
    "data", "id", "iv", "dv",
    "within", "between", "design"
  )

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
  covariate <- if (!is.null(experiment$covariate)) experiment$covariate else character(0)

  if (!is.data.frame(data)) {
    stop("experiment$dataはデータフレームである必要があります。", call. = FALSE)
  }

  if (length(id) != 1) {
    stop("experiment$idには参加者ID列を1つ指定してください。", call. = FALSE)
  }

  if (length(iv) == 0) {
    stop("独立変数が指定されていません。", call. = FALSE)
  }

  if (length(dv) == 0) {
    stop("従属変数が指定されていません。", call. = FALSE)
  }

  if (!is.numeric(ci) || length(ci) != 1 || is.na(ci) || ci <= 0 || ci >= 1) {
    stop("ciには0より大きく1より小さい数値を指定してください。", call. = FALSE)
  }

  if (!is.numeric(digits) || length(digits) != 1 || is.na(digits) ||
      digits < 0 || digits != floor(digits)) {
    stop("digitsには0以上の整数を指定してください。", call. = FALSE)
  }

  if (!is.logical(show) || length(show) != 1 || is.na(show)) {
    stop("showにはTRUEまたはFALSEを指定してください。", call. = FALSE)
  }

  missing_columns <- setdiff(unique(c(id, iv, dv, covariate)), names(data))
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

  invisible(TRUE)
}


# ------------------------------------------------------------
# 内部関数：閲覧用describeオブジェクトの作成
#
# p_describe()で計算済みの結果を再計算せず、View()しやすい
# 従属変数を第1階層、要因・交互作用を第2階層として整理する。
# ------------------------------------------------------------

.p_make_describe_object <- function(results, iv, digits) {
  sets <- .p_factor_sets(iv)

  describe <- lapply(
    names(results),
    function(dv_name) {
      result <- results[[dv_name]]

      # 各従属変数の中に、overall、各要因・交互作用、cellsを格納する。
      dv_describe <- list(
        overall = .p_round_descriptive_table(result$overall, digits)
      )

      # 各要因・交互作用は、A、B、A_Bのような単純名で格納する。
      for (set_name in names(sets)) {
        table <- result$tables[[set_name]]

        if (!is.null(table)) {
          view_name <- gsub(":", "_", set_name, fixed = TRUE)
          dv_describe[[view_name]] <- .p_round_descriptive_table(table, digits)
        }
      }

      dv_describe$cells <- .p_round_descriptive_table(result$cells, digits)
      dv_describe
    }
  )

  names(describe) <- names(results)
  describe
}


# ------------------------------------------------------------
# 公開関数：記述統計の作成
# ------------------------------------------------------------

p_describe <- function(
  experiment,
  ci = 0.95,
  digits = 3,
  show = TRUE
) {
  .p_validate_describe_input(
    experiment = experiment,
    ci = ci,
    digits = digits,
    show = show
  )

  data <- experiment$data
  id <- experiment$id
  iv <- experiment$iv
  dv <- experiment$dv
  covariate <- if (!is.null(experiment$covariate)) experiment$covariate else character(0)

  sets <- .p_factor_sets(iv)
  full_key <- paste(iv, collapse = ":")
  results <- list()

  for (dv_name in dv) {
    overall <- .p_summarize_by(
      data = data,
      id = id,
      dv = dv_name,
      factors = character(0),
      ci = ci
    )

    tables <- list(overall = overall)

    for (set_name in names(sets)) {
      tables[[set_name]] <- .p_summarize_by(
        data = data,
        id = id,
        dv = dv_name,
        factors = sets[[set_name]],
        ci = ci
      )
    }

    raw_data <- data[, unique(c(id, iv, covariate, dv_name)), drop = FALSE]
    names(raw_data)[names(raw_data) == dv_name] <- ".value"
    raw_data$dependent_variable <- dv_name

    result <- list(
      analysis = "descriptive",
      dependent_variable = dv_name,
      design = experiment$design,
      id = id,
      iv = iv,
      within = experiment$within,
      between = experiment$between,
      covariate = covariate,
      ci = ci,
      basis = "participant_means",
      overall = overall,
      tables = tables,
      cells = tables[[full_key]],
      raw_data = raw_data
    )

    # 既存の参照方法を維持するため、A、A:Bなども直接格納する。
    for (set_name in names(sets)) {
      result[[set_name]] <- tables[[set_name]]
    }

    class(result) <- c("analysis_descriptive_result", "list")
    results[[dv_name]] <- result

    if (show) {
      cat("\n========================================\n")
      cat("従属変数：", dv_name, "\n", sep = "")
      cat("実験計画：", experiment$design, "\n", sep = "")
      cat("集計単位：参加者平均\n")
      cat("信頼区間：", ci * 100, "%\n", sep = "")
      cat("========================================\n")

      cat("\n【全体】\n")
      print(
        .p_round_descriptive_table(overall, digits),
        row.names = FALSE
      )

      for (set_name in names(sets)) {
        cat(
          "\n【",
          gsub(":", " × ", set_name, fixed = TRUE),
          "】\n",
          sep = ""
        )
        print(
          .p_round_descriptive_table(tables[[set_name]], digits),
          row.names = FALSE
        )
      }
    }
  }

  attr(results, "stage") <- "preprocess"
  attr(results, "type") <- "descriptive"
  attr(results, "ci") <- ci
  attr(results, "experiment") <- experiment
  class(results) <- c("analysis_descriptive", "list")

  invisible(results)
}
