# 置換ANOVA / ANCOVA。検定は permuco::aovperm() に委譲する。
# https://doi.org/10.18637/jss.v099.i15
# 共変量は参加者ごとの連続値。平均中心化して投入する。
# covariate_interactions = TRUE は各共変量と実験要因の交互作用を含む。
# 複数共変量同士の交互作用は含めない。
t_permuco <- function(
    experiment, alpha = 0.05, followup_alpha = 0.10,
    digits = 3, p_digits = 3, np = 5000, seed = NULL,
    permutation_method = NULL, use_covariates = FALSE,
    covariate_interactions = TRUE,
    .on_result = NULL, .followup_seed = seed
) {
  if (!requireNamespace("permuco", quietly = TRUE)) {
    stop('permuco が必要です。install.packages("permuco") を実行してください。', call. = FALSE)
  }
  needed <- c("data", "id", "iv", "dv", "within", "between", "design")
  if (length(setdiff(needed, names(experiment)))) {
    stop("experimentのdata / id / iv / dv / within / between / designを確認してください。", call. = FALSE)
  }
  scalar_number <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!scalar_number(alpha) || !scalar_number(followup_alpha) ||
      alpha <= 0 || followup_alpha < alpha || followup_alpha >= 1) {
    stop("0 < alpha <= followup_alpha < 1にしてください。", call. = FALSE)
  }
  if (!scalar_number(np) || np < 2 || np != floor(np)) {
    stop("npには2以上の整数を指定してください。", call. = FALSE)
  }
  if (!scalar_number(digits) || digits < 0 || digits != floor(digits) ||
      !scalar_number(p_digits) || p_digits < 1 || p_digits != floor(p_digits)) {
    stop("digitsは0以上、p_digitsは1以上の整数にしてください。", call. = FALSE)
  }
  for (flag in list(use_covariates, covariate_interactions)) {
    if (!is.logical(flag) || length(flag) != 1L || is.na(flag)) {
      stop("共変量の指定にはTRUEまたはFALSEを使用してください。", call. = FALSE)
    }
  }
  if (!is.null(seed)) {
    if (!scalar_number(seed) || seed < 0 || seed > .Machine$integer.max || seed != floor(seed)) {
      stop("seedには0以上の整数を指定してください。", call. = FALSE)
    }
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) previous_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) assign(".Random.seed", previous_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(list = ".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
  }
  id <- experiment$id
  iv <- experiment$iv
  within <- experiment$within
  between <- experiment$between
  dv <- experiment$dv
  covariate <- if (use_covariates) experiment$covariate else character(0)
  if (use_covariates && !length(covariate)) {
    stop("permuco_ancovaにはexperiment$covariateの指定が必要です。", call. = FALSE)
  }
  if (!is.data.frame(experiment$data) || length(id) != 1L ||
      !length(iv) || !length(dv) || anyDuplicated(c(id, iv, covariate)) ||
      anyDuplicated(iv) || !setequal(iv, c(within, between)) ||
      length(intersect(within, between))) {
    stop("experimentのデータと要因指定を確認してください。", call. = FALSE)
  }
  required <- unique(c(id, iv, dv, covariate))
  if (length(setdiff(required, names(experiment$data)))) {
    stop("データに必要な列がありません。", call. = FALSE)
  }
  if (any(grepl("`", required, fixed = TRUE))) {
    stop("列名にバッククォートを含めないでください。", call. = FALSE)
  }
  data <- experiment$data
  if (!nrow(data) || anyNA(data[, required, drop = FALSE])) {
    stop("分析に使用する列に欠測があります。分析対象を整えて再実行してください。", call. = FALSE)
  }
  for (name in c(dv, covariate)) {
    if (!is.numeric(data[[name]]) || any(!is.finite(data[[name]]))) {
      stop("従属変数と共変量には有限の数値が必要です: ", name, call. = FALSE)
    }
  }
  data[c(id, iv)] <- lapply(data[c(id, iv)], factor)
  if (any(vapply(data[iv], nlevels, integer(1)) < 2L)) {
    stop("各実験要因には2水準以上が必要です。", call. = FALSE)
  }
  for (name in c(between, covariate)) {
    counts <- vapply(split(data[[name]], data[[id]], drop = TRUE),
                     function(x) length(unique(x)), integer(1))
    if (any(counts != 1L)) {
      stop("参加者ごとに一定の値を指定してください: ", name, call. = FALSE)
    }
  }
  centers <- numeric(0)
  for (name in covariate) {
    centers[name] <- mean(data[[name]][!duplicated(data[[id]])])
    data[[name]] <- data[[name]] - centers[name]
    if (length(unique(data[[name]])) < 2L) {
      stop("共変量が一定値です: ", name, call. = FALSE)
    }
  }
  # 各条件が揃っていることを確認。analyze()はセル内反復を事前に平均する。
  if (length(within)) {
    cell <- interaction(data[within], drop = FALSE)
    counts <- table(data[[id]], cell)
    if (any(counts == 0L)) {
      stop("参加者ごとに全被験者内条件の測定が必要です。", call. = FALSE)
    }
  }
  selected_method <- if (length(within)) "Rd_kheradPajouh_renaud" else "freedman_lane"
  allowed <- if (length(within)) {
    c("Rd_kheradPajouh_renaud", "Rde_kheradPajouh_renaud")
  } else c("freedman_lane", "draper_stoneman", "manly", "terBraak", "kennedy", "huh_jhun", "dekker")
  if (!is.null(permutation_method)) {
    selected_method <- match.arg(permutation_method, allowed)
  }
  quote_name <- function(x) paste0("`", x, "`")
  fixed <- paste(quote_name(iv), collapse = " * ")
  if (length(covariate)) {
    cov_terms <- quote_name(covariate)
    if (covariate_interactions) cov_terms <- paste0("(", fixed, ") * ", cov_terms)
    fixed <- paste(c(fixed, cov_terms), collapse = " + ")
  }
  error <- if (length(within)) paste0(" + Error(", quote_name(id), " / (",
                   paste(quote_name(within), collapse = " * "), "))") else ""
  results <- list()
  # Identical within-cell ordering is required for whole-subject residual permutations.
  data <- data[do.call(order, data[c(id, within)]), , drop = FALSE]
  for (target in dv) {
    formula <- stats::as.formula(paste0(quote_name(target), " ~ ", fixed, error))
    model <- permuco::aovperm(formula, data = data, np = np,
                             method = selected_method, coding_sum = TRUE)
    table <- as.data.frame(model$table, check.names = FALSE)
    if (!nrow(table) || !"resampled P(>F)" %in% names(table)) {
      stop("permucoの分散分析表の構造を確認できません。", call. = FALSE)
    }
    table$Effect <- rownames(table)
    table$p.value <- table[["resampled P(>F)"]]
    results[[target]] <- list(
      analysis = if (use_covariates) "permuco_ancova" else "permuco",
      dependent_variable = target, design = experiment$design,
      within = within, between = between, covariate = covariate,
      covariate_adjusted = use_covariates, covariate_centers = centers,
      covariate_interactions = use_covariates && covariate_interactions,
      analysis_data = data[, unique(c(id, iv, covariate, target)), drop = FALSE],
      formula = formula, model = model, anova_table = table,
      anova_table_uncorrected = table,
      anova_type = "III", permutation_method = selected_method,
      np = model$np, requested_np = np, seed = seed, alpha = alpha, followup_alpha = followup_alpha,
      contrasts = list(), digits = digits, p_digits = p_digits,
      followup_factors = iv,
      compute_followups = local({
        followup_data <- data; followup_formula <- formula
        engine <- NULL
        function(effects, forced, main_effects) {
          if (is.null(engine)) engine <<- .t_permuco_followups(followup_data,
            followup_formula, id, iv, within, covariate, np, .followup_seed)
          engine(effects, forced, main_effects)
        }
      }),
      runtime_info = list(R = R.version.string,
                          permuco = as.character(utils::packageVersion("permuco")))
    )
    class(results[[target]]) <- c("analysis_permuco_result", "analysis_test_result", "list")
    if (!is.null(.on_result)) results[[target]] <- .on_result(results[[target]])
  }
  class(results) <- c("analysis_permuco", "analysis_test", "list")
  invisible(results)
}
