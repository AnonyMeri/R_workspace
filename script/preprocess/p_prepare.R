# ============================================================
# 分析前処理：データ整形
#
# 公開関数：
#   p_prepare()
#
# 入力データを検証・整形し、後続処理で共通して使用する
# analysis_experimentオブジェクトを作成します。
#
# 各独立変数について、同一参加者ID内で複数水準を取るかを確認し、
# 被験者内要因（within）/ 被験者間要因（between）を自動推定します。
# その推定結果から実験計画（between / within / mixed）を判定します。
#
# covariateはANCOVA等で使用する連続共変量として、ivとは別に
# numericのまま保持します。NULLまたは複数列を指定できます。
# ============================================================

.p_check_column_names <- function(x, name, allow_empty = FALSE) {
  if (is.null(x) && !allow_empty) {
    stop(name, "を指定するか、preprocess(data, select = \"gui\")で選択してください。", call. = FALSE)
  }
  if (!is.character(x) || anyNA(x) || any(!nzchar(x))) {
    stop(name, "には列名を文字列で指定してください。", call. = FALSE)
  }
  if (!allow_empty && length(x) == 0) {
    stop(name, "には少なくとも1つの列を指定してください。", call. = FALSE)
  }
  invisible(TRUE)
}

.p_validate_prepare_data <- function(data) {
  if (!is.data.frame(data)) {
    stop("dataはデータフレームである必要があります。", call. = FALSE)
  }
  if (!nrow(data) || !ncol(data)) {
    stop("dataには1行以上・1列以上のデータが必要です。", call. = FALSE)
  }
  if (anyNA(names(data)) || any(!nzchar(names(data))) || anyDuplicated(names(data))) {
    stop("データの列名には空欄や重複のない名前を指定してください。", call. = FALSE)
  }
  invisible(TRUE)
}

.p_validate_prepare_input <- function(data, id, iv, dv, covariate) {
  .p_validate_prepare_data(data)

  .p_check_column_names(id, "id")
  .p_check_column_names(iv, "iv")
  .p_check_column_names(dv, "dv")
  .p_check_column_names(covariate, "covariate", allow_empty = TRUE)

  if (length(id) != 1) {
    stop("idには参加者ID列を1つ指定してください。", call. = FALSE)
  }

  required_columns <- unique(c(id, iv, dv, covariate))
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0) {
    stop(
      "次の列がデータ内に存在しません: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  role_sets <- list(id = id, iv = iv, dv = dv, covariate = covariate)
  duplicated_roles <- character(0)
  role_names <- names(role_sets)
  if (length(role_names) >= 2) {
    for (i in seq_len(length(role_names) - 1L)) {
      for (j in (i + 1L):length(role_names)) {
        duplicated_roles <- c(
          duplicated_roles,
          intersect(role_sets[[i]], role_sets[[j]])
        )
      }
    }
  }
  duplicated_roles <- unique(duplicated_roles)

  if (length(duplicated_roles) > 0) {
    stop(
      "次の列が複数の役割に指定されています: ",
      paste(duplicated_roles, collapse = ", "),
      call. = FALSE
    )
  }

  if (anyNA(data[[id]])) {
    stop("参加者ID列に欠損値があります。", call. = FALSE)
  }

  factor_missing <- iv[vapply(data[iv], anyNA, logical(1))]
  if (length(factor_missing) > 0) {
    stop(
      "次の独立変数に欠損値があります: ",
      paste(factor_missing, collapse = ", "),
      call. = FALSE
    )
  }

  non_numeric_dv <- dv[!vapply(data[dv], is.numeric, logical(1))]
  if (length(non_numeric_dv) > 0) {
    stop(
      "次の従属変数が数値型ではありません: ",
      paste(non_numeric_dv, collapse = ", "),
      call. = FALSE
    )
  }

  non_numeric_covariate <- covariate[
    !vapply(data[covariate], is.numeric, logical(1))
  ]
  if (length(non_numeric_covariate) > 0) {
    stop(
      "次の共変量が数値型ではありません: ",
      paste(non_numeric_covariate, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.p_detect_factor_roles <- function(data, id, iv) {
  is_within <- vapply(
    iv,
    function(factor_name) {
      level_counts <- tapply(
        as.character(data[[factor_name]]),
        data[[id]],
        function(x) length(unique(x))
      )

      any(level_counts > 1, na.rm = TRUE)
    },
    logical(1)
  )

  list(
    within = iv[is_within],
    between = iv[!is_within]
  )
}

.p_detect_design <- function(within, between) {
  if (length(within) > 0 && length(between) > 0) return("mixed")
  if (length(within) > 0) return("within")
  "between"
}

p_prepare <- function(
  data,
  id,
  iv,
  dv,
  covariate = NULL
) {
  if (is.null(covariate)) covariate <- character(0)

  iv <- unique(iv)
  dv <- unique(dv)
  covariate <- unique(covariate)

  .p_validate_prepare_input(
    data = data,
    id = id,
    iv = iv,
    dv = dv,
    covariate = covariate
  )

  analysis_data <- data
  analysis_data[[id]] <- factor(analysis_data[[id]])
  analysis_data[iv] <- lapply(analysis_data[iv], factor)

  factor_roles <- .p_detect_factor_roles(
    data = analysis_data,
    id = id,
    iv = iv
  )
  within <- factor_roles$within
  between <- factor_roles$between

  order_index <- do.call(
    order,
    c(analysis_data[c(id, iv)], list(na.last = TRUE))
  )
  analysis_data <- analysis_data[order_index, , drop = FALSE]
  rownames(analysis_data) <- NULL

  factor_levels <- stats::setNames(lapply(analysis_data[iv], levels), iv)

  experiment <- list(
    data = analysis_data,
    id = id,
    iv = iv,
    dv = dv,
    within = within,
    between = between,
    covariate = covariate,
    design = .p_detect_design(within, between),
    factor_levels = factor_levels,
    n_participants = nlevels(analysis_data[[id]])
  )

  class(experiment) <- c("analysis_experiment", "list")
  experiment
}
