# ============================================================
# 回帰分析の計算処理。
# helperが方法・対象・引数を確定し、次の一方だけを呼び出す。
#   r_regression()             : 通常の線形回帰・既存モデルの傾き分析
#   r_permutation_regression() : 指定済みの設計と置換設定による傾き分析
# GUI・表示・有意な対象の抽出・設定の継承・多重比較補正はhelperが担当。
# 計算はlm / confint / emtrends / lmpermへ委譲し、未補正p値・CIを返す。
# ============================================================

# 共通の推定処理
.r_check_estimation <- function(ci, compare_slopes) {
  if (!is.numeric(ci) || length(ci) != 1L || !is.finite(ci) || ci <= 0 || ci >= 1) {
    stop("ciは0より大きく1より小さい数値にしてください。", call. = FALSE)
  }
  if (!is.logical(compare_slopes) || length(compare_slopes) != 1L || is.na(compare_slopes)) {
    stop("compare_slopesにはTRUEまたはFALSEを指定してください。", call. = FALSE)
  }
}

.r_estimates <- function(grid, ci, trend = NULL) {
  table <- as.data.frame(summary(grid, infer = c(TRUE, TRUE), adjust = "none", level = ci))
  if (!is.null(trend)) {
    column <- paste0(trend, ".trend")
    if (!column %in% names(table)) stop("傾きの推定結果を取得できませんでした。", call. = FALSE)
    names(table)[names(table) == column] <- "estimate"
  }
  table
}

.r_regression_result <- function(..., permutation = NULL) {
  result <- list(analysis = "regression", ...)
  classes <- c("analysis_regression_result", "list")
  if (!is.null(permutation)) {
    result$inference <- "permutation"
    result$permutation <- permutation
    result$confidence_interval_method <- "parametric_linear_model"
    classes <- c("analysis_permutation_regression_result", classes)
  }
  structure(result, class = classes)
}

.r_quote_name <- function(x) paste0("`", gsub("`", "", x, fixed = TRUE), "`")

.r_lm_tables <- function(model, ci) {
  analysis_data <- stats::model.frame(model)
  if (stats::df.residual(model) <= 0L) stop("残差自由度がありません。説明変数の数や分析対象数を確認してください。", call. = FALSE)
  response <- stats::model.response(analysis_data)
  if (!is.numeric(response) || is.matrix(response)) stop("単一の数値の目的変数を持つモデルを指定してください。", call. = FALSE)
  sm <- summary(model)
  coefficients <- stats::coef(sm)
  interval <- stats::confint(model, level = ci)
  terms <- names(stats::coef(model))
  coefficient_table <- data.frame(term = terms, estimate = unname(stats::coef(model)),
    SE = coefficients[match(terms, rownames(coefficients)), "Std. Error"],
    df = stats::df.residual(model),
    t.ratio = coefficients[match(terms, rownames(coefficients)), "t value"],
    p.value = coefficients[match(terms, rownames(coefficients)), "Pr(>|t|)"],
    lower.CL = interval[match(terms, rownames(interval)), 1L],
    upper.CL = interval[match(terms, rownames(interval)), 2L], row.names = NULL)
  f <- sm$fstatistic
  model_table <- data.frame(N = stats::nobs(model), R2 = sm$r.squared,
    adjusted.R2 = sm$adj.r.squared, F = if (length(f)) unname(f[1L]) else NA_real_,
    df1 = if (length(f)) unname(f[2L]) else NA_real_, df2 = stats::df.residual(model),
    p.value = if (length(f)) stats::pf(f[1L], f[2L], f[3L], lower.tail = FALSE) else NA_real_)
  list(coefficient_table = coefficient_table, model_table = model_table)
}

# 通常回帰：データからの推定、または既存lm / afex_aovの利用
r_regression <- function(data = NULL, formula = NULL, model = NULL,
                         trend = NULL, by = NULL, compare_slopes = TRUE, ci = .95) {
  .r_check_estimation(ci, compare_slopes)
  if (!is.null(model) && !is.null(formula)) {
    stop("formulaとmodelは、どちらか一方を指定してください。", call. = FALSE)
  }
  supplied_model <- !is.null(model)
  if (!supplied_model) {
    if (!is.data.frame(data) || !inherits(formula, "formula") || length(formula) != 3L) {
      stop("データフレームと、目的変数を含むモデル式を指定してください。", call. = FALSE)
    }
    response <- stats::model.response(stats::model.frame(formula, data, na.action = stats::na.pass))
    if (!is.numeric(response) || is.matrix(response) || any(is.infinite(response))) {
      stop("目的変数には有限の数値を指定してください（欠測は除外します）。", call. = FALSE)
    }
    model <- stats::lm(formula, data = data, na.action = stats::na.exclude, model = TRUE)
  }
  ordinary <- inherits(model, "lm") && !inherits(model, c("glm", "mlm", "aov"))
  repeated <- inherits(model, "afex_aov")
  if (!ordinary && !repeated) {
    stop("modelにはlmまたはafex_aovを指定してください。置換ANCOVAはrelate(分析結果)で指定してください。", call. = FALSE)
  }
  if (repeated && is.null(trend)) {
    stop("afex_aovモデルを使う場合は、傾きを調べる共変量をtrendで指定してください。", call. = FALSE)
  }
  if (!is.null(data) && !is.data.frame(data)) stop("dataにはデータフレームを指定してください。", call. = FALSE)
  analysis_data <- if (ordinary) stats::model.frame(model) else model$data$long
  reference_data <- if (!is.null(data)) data else analysis_data
  tables <- if (ordinary) .r_lm_tables(model, ci) else
    list(coefficient_table = NULL, model_table = NULL)
  if (length(by) && is.null(trend)) stop("条件別の傾きを調べる場合はtrendも指定してください。", call. = FALSE)
  slope_table <- slope_grid <- NULL
  comparisons <- list()
  if (!is.null(trend)) {
    if (!is.character(trend) || length(trend) != 1L || is.na(trend) ||
        !trend %in% names(reference_data) || !is.numeric(reference_data[[trend]])) {
      stop("trendには、モデルで使用した数値の説明変数を1つ指定してください。", call. = FALSE)
    }
    if (!is.null(by) && (!is.character(by) || anyNA(by) || anyDuplicated(by) ||
        any(!by %in% names(reference_data)))) stop("byには条件を表す列名を指定してください。", call. = FALSE)
    if (trend %in% by) stop("trendとbyには異なる変数を指定してください。", call. = FALSE)
    if (any(!vapply(reference_data[by], function(x) is.factor(x) || is.character(x), logical(1)))) {
      stop("byには因子または文字列の条件変数を指定してください。", call. = FALSE)
    }
    rhs <- if (ordinary) all.vars(stats::delete.response(stats::terms(model))) else
      unique(c(model$within, model$between, model$covariate))
    # afexの変数情報はモデルにより保持形式が異なるため、基準データも利用する。
    if (ordinary && any(!c(trend, by) %in% rhs)) stop("trendとbyにはモデルの説明変数を指定してください。", call. = FALSE)
    if (!requireNamespace("emmeans", quietly = TRUE)) stop("条件別の傾きの計算にはemmeansが必要です。", call. = FALSE)
    if (repeated && !requireNamespace("afex", quietly = TRUE)) stop("afex_aovモデルの利用にはafexが必要です。", call. = FALSE)
    slope_grid <- emmeans::emtrends(model, specs = if (length(by)) by else ~ 1,
      var = trend, data = reference_data)
    slope_table <- .r_estimates(slope_grid, ci, trend)
    if (compare_slopes && length(by)) for (factor in by) {
      grid <- emmeans::contrast(slope_grid, method = "pairwise", by = setdiff(by, factor), adjust = "none")
      comparisons[[factor]] <- list(table = .r_estimates(grid, ci), grid = grid,
        compared_factor = factor, fixed_factors = setdiff(by, factor))
    }
  }
  .r_regression_result(dependent_variable = if (ordinary)
      paste(deparse(stats::formula(model)[[2L]]), collapse = " ") else attr(model, "dv"),
    source = if (supplied_model) "existing_model" else "data", model = model,
    model_class = class(model), analysis_data = analysis_data, formula = if (ordinary) stats::formula(model) else model$call,
    coefficient_table = tables$coefficient_table, model_table = tables$model_table,
    trend = trend, by = by, slope_table = slope_table, slope_grid = slope_grid,
    slope_comparisons = comparisons, ci = ci,
    diagnostics = if (ordinary) list(aliased = names(stats::coef(model))[is.na(stats::coef(model))],
      residuals = stats::residuals(model), fitted = stats::fitted(model),
      excluded_rows = model$na.action) else list())
}

# 置換回帰用の応答：参加者ごとの条件平均・条件差を作るためのデータ
.r_permutation_design <- function(data, formula, design, trend, by) {
  target <- design$dv
  within <- design$within
  between <- design$between
  covariates <- design$covariate
  if (!is.data.frame(data) || length(target) != 1L || !target %in% names(data) ||
      length(trend) != 1L || !trend %in% covariates || any(!by %in% c(within, between))) {
    stop("置換ANCOVAの分析データ・共変量・条件の情報を確認してください。", call. = FALSE)
  }
  needed <- unique(c(target, within, between, covariates))
  if (any(!needed %in% names(data)) || anyNA(data[needed])) stop("置換回帰に必要な列がないか、欠測があります。", call. = FALSE)
  if (any(!vapply(data[c(target, covariates)], function(x) is.numeric(x) && all(is.finite(x)), logical(1)))) {
    stop("目的変数・共変量には有限の数値が必要です。", call. = FALSE)
  }
  data[c(within, between)] <- lapply(data[c(within, between)], factor)
  id <- design$id
  if (length(within) && (length(id) != 1L || !id %in% names(data))) {
    stop("反復測定には有効な参加者IDが必要です。", call. = FALSE)
  }
  if (is.null(id)) {
    id <- ".relation_subject"
    while (id %in% names(data)) id <- paste0(id, "_")
    data[[id]] <- seq_len(nrow(data))
  }
  if (anyNA(data[[id]])) stop("参加者IDに欠測があります。", call. = FALSE)
  ids <- unique(as.character(data[[id]]))
  first <- match(ids, as.character(data[[id]]))
  for (name in c(between, covariates)) {
    counts <- vapply(split(data[[name]], as.character(data[[id]])), function(x) length(unique(x)), integer(1))
    if (any(counts != 1L)) stop("参加者内で値が変わる群間要因・共変量には対応していません：", name, call. = FALSE)
  }
  subject <- data[first, c(between, covariates), drop = FALSE]
  rownames(subject) <- NULL
  response_name <- ".relation_response"
  while (response_name %in% names(subject)) response_name <- paste0(response_name, "_")
  if (length(within)) {
    cells <- interaction(data[within], drop = FALSE, lex.order = TRUE)
    cell_index <- as.integer(cells)
    grid <- data[match(seq_len(nlevels(cells)), cell_index), within, drop = FALSE]
    rownames(grid) <- NULL
    means <- stats::aggregate(data[[target]], list(subject = as.character(data[[id]]), cell = cell_index), mean)
    responses <- matrix(NA_real_, length(ids), nlevels(cells))
    responses[cbind(match(means$subject, ids), means$cell)] <- means$x
    if (anyNA(responses)) stop("参加者ごとに全被験者内条件の測定が必要です。", call. = FALSE)
  } else {
    if (anyDuplicated(as.character(data[[id]]))) stop("対応なしの置換回帰には参加者ごとに1行のデータが必要です。", call. = FALSE)
    responses <- matrix(data[[target]][first], ncol = 1L)
    grid <- data.frame(.overall = 1L)
  }
  # 元の固定効果を条件別応答へ射影し、群間要因と他の共変量を保持する。
  allowed <- c(between, covariates)
  original_terms <- attr(stats::terms(formula, specials = "Error"), "term.labels")
  original_terms <- original_terms[!grepl("^Error\\(", original_terms)]
  projected <- lapply(original_terms, function(term) {
    vars <- all.vars(parse(text = term))
    if (any(!vars %in% c(within, allowed))) stop("親モデルに未対応の説明変数があります。", call. = FALSE)
    allowed[allowed %in% vars]
  })
  rhs <- unique(vapply(projected[lengths(projected) > 0L], function(x) paste(.r_quote_name(x), collapse = ":"), character(1)))
  if (!length(rhs)) stop("置換回帰の説明変数を作成できません。", call. = FALSE)
  regression_formula <- stats::as.formula(paste(.r_quote_name(response_name), "~", paste(rhs, collapse = " + ")))
  for (name in between) contrasts(subject[[name]]) <- stats::contr.sum(nlevels(subject[[name]]))
  list(data = data, subject = subject, responses = responses, grid = grid,
    subject_ids = ids, response_name = response_name, formula = regression_formula)
}

# 置換行列は分析ごとに1つ確定し、傾き・傾き差の全検定で共用する。
.r_permutation_matrix <- function(parent, n, np, type, reuse) {
  choices <- if (is.matrix(parent)) list(parent) else if (is.list(parent)) parent else list()
  eligible <- choices[vapply(choices, function(x) is.matrix(x) && nrow(x) == n && ncol(x) == np, logical(1))]
  reused <- reuse && length(eligible) > 0L
  list(P = if (reused) eligible[[1L]] else permuco::Pmat(n = n, np = np, type = type), reused = reused)
}

# 1つの線形仮説Kβを1係数に変換し、同じモデル空間でlmpermに渡す。
.r_permutation_hypothesis <- function(base, row, K, weights, response_name, method, P) {
  X <- base$X
  pivot <- which.max(abs(K))
  if (!length(pivot) || !is.finite(K[pivot]) || abs(K[pivot]) < 1e-12) {
    row$p.value <- NA_real_
    return(list(row = row, test = NULL))
  }
  transform <- diag(ncol(X)); transform[pivot, ] <- K
  changed <- X %*% solve(transform)
  columns <- paste0(".relation_b", seq_len(ncol(changed)))
  frame <- as.data.frame(changed); names(frame) <- columns
  frame[[response_name]] <- stats::model.response(stats::model.frame(base$fit))
  f <- stats::as.formula(paste(.r_quote_name(response_name), "~ 0 +", paste(columns, collapse = " + ")))
  permutation_fit <- permuco::lmperm(f, data = frame, P = P, method = method)
  table <- permutation_fit$table
  name <- columns[pivot]
  if (!name %in% rownames(table) || !"resampled Pr(>|t|)" %in% names(table)) stop("置換回帰の検定結果を取得できません。", call. = FALSE)
  if (!isTRUE(all.equal(as.numeric(row$estimate), as.numeric(table[name, "Estimate"]), tolerance = 1e-7))) {
    stop("傾きの推定値と置換検定の仮説が一致しません。", call. = FALSE)
  }
  row$p.parametric <- row$p.value
  row$p.value <- as.numeric(table[name, "resampled Pr(>|t|)"])
  test <- list(weights = weights, linear_contrast = K,
    coefficient = name, observed = table[name, , drop = FALSE], distribution = permutation_fit$distribution[, name],
    formula = f, lm = base$fit)
  list(row = row, test = test)
}

# 置換回帰：helperが確定した設計・置換法・回数・seedをそのまま使う。
r_permutation_regression <- function(data, formula, model, design, permutation,
                                     trend, by, compare_slopes = TRUE, ci = .95) {
  .r_check_estimation(ci, compare_slopes)
  if (!requireNamespace("permuco", quietly = TRUE) || !requireNamespace("emmeans", quietly = TRUE)) {
    stop("置換回帰の傾き分析にはpermucoとemmeansが必要です。", call. = FALSE)
  }
  prepared <- .r_permutation_design(data, formula, design, trend, by)
  data <- prepared$data
  subject <- prepared$subject
  responses <- prepared$responses
  grid <- prepared$grid
  response_name <- prepared$response_name
  regression_formula <- prepared$formula
  within <- design$within
  np <- permutation$np
  method <- permutation$method
  seed <- permutation$seed
  type <- permutation$type
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
    if (had_seed) previous_seed <- get(".Random.seed", .GlobalEnv)
    on.exit({if (had_seed) assign(".Random.seed", previous_seed, .GlobalEnv) else if
      (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)}, add = TRUE)
    set.seed(seed)
  }
  fits <- new.env(parent = emptyenv())
  tests <- list()
  within_by <- intersect(by, within)
  between_by <- intersect(by, design$between)
  weights_for <- function(cell) {
    keep <- rep(TRUE, nrow(grid))
    for (name in within_by) keep <- keep & as.character(grid[[name]]) == as.character(cell[[name]])
    if (!any(keep)) stop("指定した被験者内条件がありません。", call. = FALSE)
    as.numeric(keep) / sum(keep)
  }
  fit_for <- function(weights) {
    key <- paste(format(weights, digits = 17), collapse = ",")
    if (exists(key, fits, inherits = FALSE)) return(get(key, fits))
    frame <- subject
    frame[[response_name]] <- as.vector(responses %*% weights)
    fit <- stats::lm(regression_formula, frame, na.action = stats::na.fail)
    X <- stats::model.matrix(fit)
    if (fit$rank != ncol(X) || stats::df.residual(fit) <= 0L) stop("条件別回帰のモデルを一意に推定できません。", call. = FALSE)
    emm <- emmeans::emtrends(fit, specs = if (length(between_by)) between_by else ~ 1,
      var = trend, data = frame)
    result <- list(fit = fit, frame = frame, grid = emm, X = X, table = .r_estimates(emm, ci, trend))
    assign(key, result, fits)
    result
  }
  test_one <- function(left, right = NULL, compared_factor = NULL) {
    weights <- weights_for(left)
    if (!is.null(right) && compared_factor %in% within) weights <- weights - weights_for(right)
    base <- fit_for(weights)
    reference <- base$table
    find_row <- function(cell) {
      selected <- rep(TRUE, nrow(reference))
      for (name in between_by) selected <- selected & as.character(reference[[name]]) == as.character(cell[[name]])
      index <- which(selected)
      if (length(index) != 1L) stop("傾きに対応する条件を特定できません。", call. = FALSE)
      index
    }
    index <- find_row(left)
    if (!is.null(right) && compared_factor %in% design$between) {
      contrast <- rep(0, nrow(reference)); contrast[index] <- 1; contrast[find_row(right)] <- -1
      query <- emmeans::contrast(base$grid, list(difference = contrast), adjust = "none")
      row <- .r_estimates(query, ci)
      K <- as.numeric(methods::slot(query, "linfct")[1L, ])
    } else {
      row <- reference[index, , drop = FALSE]
      K <- as.numeric(methods::slot(base$grid, "linfct")[index, ])
    }
    tested <- .r_permutation_hypothesis(base, row, K, weights, response_name, method, permutations$P)
    if (!is.null(tested$test)) tests[[length(tests) + 1L]] <<- tested$test
    tested$row
  }
  cells <- if (length(by)) expand.grid(lapply(data[by], levels), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) else data.frame(.overall = 1L)
  initial <- fit_for(weights_for(cells[1L, , drop = FALSE]))
  n <- nrow(initial$X)
  if (method == "huh_jhun") n <- n - ncol(initial$X) + 1L
  permutations <- .r_permutation_matrix(permutation$P, n, np, type, reuse = !length(within))
  slopes <- vector("list", nrow(cells))
  stats_columns <- c("estimate", "SE", "df", "t.ratio", "lower.CL", "upper.CL", "p.value")
  for (i in seq_len(nrow(cells))) {
    row <- test_one(cells[i, , drop = FALSE])
    slopes[[i]] <- cbind(cells[i, by, drop = FALSE], row[, intersect(stats_columns, names(row)), drop = FALSE])
  }
  comparisons <- list()
  if (compare_slopes) for (factor in by) {
    fixed <- setdiff(by, factor)
    groups <- if (length(fixed)) unique(cells[fixed]) else data.frame(.overall = 1L)
    compared <- list()
    for (g in seq_len(nrow(groups))) {
      keep <- rep(TRUE, nrow(cells))
      for (name in fixed) keep <- keep & cells[[name]] == groups[[name]][g]
      rows <- which(keep)
      if (length(rows) < 2L) next
      pairs <- utils::combn(rows, 2L)
      for (j in seq_len(ncol(pairs))) {
        left <- cells[pairs[1L, j], , drop = FALSE]; right <- cells[pairs[2L, j], , drop = FALSE]
        row <- test_one(left, right, factor)
        labels <- data.frame(contrast = paste(left[[factor]], right[[factor]], sep = " - "), stringsAsFactors = FALSE)
        compared[[length(compared) + 1L]] <- cbind(labels, left[, fixed, drop = FALSE], row[, intersect(stats_columns, names(row)), drop = FALSE])
      }
    }
    comparisons[[factor]] <- list(table = do.call(rbind, compared), grid = NULL,
      compared_factor = factor, fixed_factors = fixed)
  }
  .r_regression_result(dependent_variable = design$dv,
    source = "existing_model", model = model, model_class = class(model),
    analysis_data = data, formula = formula, coefficient_table = NULL, model_table = NULL,
    trend = trend, by = by, slope_table = do.call(rbind, slopes), slope_grid = NULL,
    slope_comparisons = comparisons, ci = ci,
    permutation = list(method = method, parent_method = permutation$parent_method, np = np, seed = seed,
      type = type, P = permutations$P, reused_parent_P = permutations$reused, unit = "participant",
      repeated_measures = length(within) > 0L, tests = tests),
    diagnostics = list(subject_ids = prepared$subject_ids, within_cell_means = responses,
      between_formula = regression_formula))
}
