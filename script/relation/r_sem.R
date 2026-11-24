# ============================================================
# SEMの計算処理。モデル選択・GUI・表示・p値の補正はhelperが担当。
# 単一集団・個票データのパス解析／測定モデル／構造モデルをlavaan::semへ委譲する。
# 推定値・検定・CI・適合度・R²を標準の抽出関数から取得する。
# ============================================================

.r_sem_fit_table <- function(measures) {
  get <- function(name) if (name %in% names(measures)) unname(measures[[name]]) else NA_real_
  rows <- list()
  for (suffix in c("", ".scaled", ".robust")) {
    if (nzchar(suffix) && !any(paste0(c("chisq", "cfi", "tli", "rmsea"), suffix) %in% names(measures))) next
    row <- data.frame(variant = if (!nzchar(suffix)) "通常" else if (suffix == ".scaled") "尺度補正" else "頑健",
      chisq = get(paste0("chisq", suffix)), df = get(paste0("df", suffix)), p.value = get(paste0("pvalue", suffix)),
      CFI = get(paste0("cfi", suffix)), TLI = get(paste0("tli", suffix)), RMSEA = get(paste0("rmsea", suffix)),
      RMSEA.lower = get(paste0("rmsea.ci.lower", suffix)), RMSEA.upper = get(paste0("rmsea.ci.upper", suffix)),
      SRMR = if (!nzchar(suffix)) get("srmr") else NA_real_,
      AIC = if (!nzchar(suffix)) get("aic") else NA_real_, BIC = if (!nzchar(suffix)) get("bic") else NA_real_,
      stringsAsFactors = FALSE)
    rows[[length(rows) + 1L]] <- row
  }
  do.call(rbind, rows)
}

r_sem <- function(data, model, estimator = "ML", missing_data = "listwise", ordered = NULL,
                   std_lv = FALSE, fixed_x = TRUE, meanstructure = TRUE,
                   se = "default", bootstrap = 1000, seed = NULL, ci = .95) {
  if (!requireNamespace("lavaan", quietly = TRUE)) stop('SEMにはlavaanが必要です。install.packages("lavaan")で導入してください。', call. = FALSE)
  if (!is.data.frame(data) || !nrow(data) || anyDuplicated(names(data)) ||
      !is.character(model) || !length(model) || anyNA(model) || !nzchar(trimws(paste(model, collapse = "\n")))) {
    stop("個票のデータフレームと、lavaan形式のモデル式を指定してください。", call. = FALSE)
  }
  estimator <- match.arg(estimator, c("ML", "MLR", "WLSMV"))
  missing_data <- match.arg(missing_data, c("listwise", "fiml"))
  se <- match.arg(se, c("default", "bootstrap"))
  for (flag in list(std_lv, fixed_x, meanstructure)) if (!is.logical(flag) || length(flag) != 1L || is.na(flag)) stop("std_lv・fixed_x・meanstructureはTRUEまたはFALSEにしてください。", call. = FALSE)
  if (!is.numeric(ci) || length(ci) != 1L || !is.finite(ci) || ci <= 0 || ci >= 1) stop("ciは0より大きく1より小さい数値にしてください。", call. = FALSE)
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || seed < 0 || seed > .Machine$integer.max || seed != floor(seed)) stop("seedには0以上の整数を指定してください。", call. = FALSE)
    had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
    if (had_seed) previous_seed <- get(".Random.seed", .GlobalEnv)
    on.exit({if (had_seed) assign(".Random.seed", previous_seed, .GlobalEnv) else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)}, add = TRUE)
    set.seed(seed)
  }
  arguments <- list(model = paste(model, collapse = "\n"), data = data, estimator = estimator,
    missing = if (missing_data == "fiml") "ml" else "listwise", ordered = ordered,
    std.lv = std_lv, fixed.x = fixed_x, meanstructure = meanstructure)
  if (se == "bootstrap") {arguments$se <- "bootstrap"; arguments$bootstrap <- bootstrap}
  # sem()は呼び出し名を検査する版があるため、名前で解決して呼び出す。
  fit <- do.call("sem", arguments, envir = asNamespace("lavaan"))
  converged <- isTRUE(lavaan::lavInspect(fit, "converged"))
  errors <- list(); parameters <- measures <- fit_table <- r2 <- NULL
  admissible <- NA
  if (converged) {
    admissible <- isTRUE(lavaan::lavInspect(fit, "post.check"))
    parameters <- tryCatch(lavaan::parameterEstimates(fit, standardized = TRUE, ci = TRUE, level = ci, boot.ci.type = "perc"),
      error = function(e) {errors$parameters <<- conditionMessage(e); NULL})
    measures <- tryCatch(lavaan::fitMeasures(fit, fm.args = list(rmsea.ci.level = ci)),
      error = function(e) {errors$fit_measures <<- conditionMessage(e); NULL})
    if (length(measures)) fit_table <- .r_sem_fit_table(measures)
    squared <- tryCatch(lavaan::lavInspect(fit, "rsquare"), error = function(e) {errors$r2 <<- conditionMessage(e); NULL})
    if (length(squared)) r2 <- data.frame(variable = names(squared), R2 = unname(squared), stringsAsFactors = FALSE)
  }
  n <- sum(lavaan::lavInspect(fit, "nobs"))
  result <- list(analysis = "sem", model = fit, model_syntax = paste(model, collapse = "\n"), analysis_data = data,
    estimator = estimator, missing_data = missing_data, ordered = lavaan::lavNames(fit, "ov.ord"),
    std_lv = std_lv, fixed_x = fixed_x, meanstructure = meanstructure, se = se,
    bootstrap = if (se == "bootstrap") bootstrap else NULL, seed = seed, ci = ci,
    converged = converged, admissible = admissible, N_input = nrow(data), N = n,
    observed_variables = lavaan::lavNames(fit, "ov"), latent_variables = lavaan::lavNames(fit, "lv"),
    parameter_estimates = parameters, fit_measures = measures, fit_table = fit_table, r2_table = r2,
    options = lavaan::lavInspect(fit, "options"),
    diagnostics = list(errors = errors, rows_used = lavaan::lavInspect(fit, "case.idx")),
    runtime_info = list(R = R.version.string, lavaan = as.character(utils::packageVersion("lavaan"))))
  class(result) <- c("analysis_sem_result", "list")
  result
}
