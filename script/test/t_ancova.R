# 通常ANCOVA。t_anova.Rの共通計算部を利用し、表示・補正はtest_helperが担当。
# 必要パッケージ: afex, car, emmeans
# experiment$covariate: 参加者内で一定の数値列名（複数指定可）。
# 分析対象の参加者平均で共変量を中心化し、下位検定はその平均値で行う。
# 既定は群間の共通傾斜。TRUEは群間要因×各共変量を追加する。
# 反復測定はafexの標準仕様に従い、被験者内要因×共変量を含む。
# 共変量同士の交互作用・共変量の水準間比較は自動追加しない。
# 直接使用するときは、先に同じフォルダのt_anova.Rもsourceする。
# 通常はsource("test_helper.R")の後、analyze(..., method = "ancova")で実行。
# 例:
#   experiment$covariate <- c("age", "baseline")
#   result <- analyze(experiment, method = "ancova")
#   result <- analyze(experiment, method = c("anova", "ancova"))
# ANCOVAと置換ANCOVAを比較する場合は、共通の交互作用設定で別々に実行する。
#   a <- analyze(experiment, method = "ancova", covariate_interactions = TRUE)
#   b <- analyze(experiment, method = "permuco_ancova", covariate_interactions = TRUE, seed = 123)
# 既定ではpermuco_ancovaは全実験要因×共変量を含むため、モデルが異なる。

t_ancova <- function(
  experiment, alpha = 0.05, correction = "auto", adjust = "bonferroni",
  digits = 3, p_digits = 3, effect_size = "pes", type = 3,
  followup_alpha = 0.10, simple_effects = FALSE,
  covariate_interactions = FALSE, .on_result = NULL, .on_sphericity = NULL
) {
  if (!exists("t_anova", mode = "function")) {
    stop("t_ancova.Rと同じフォルダのt_anova.Rも読み込んでください。", call. = FALSE)
  }
  t_anova(
    experiment = experiment, alpha = alpha, correction = correction,
    adjust = adjust, digits = digits, p_digits = p_digits,
    effect_size = effect_size, type = type, followup_alpha = followup_alpha,
    simple_effects = simple_effects, .on_result = .on_result,
    .on_sphericity = .on_sphericity, .use_covariates = TRUE,
    .covariate_interactions = covariate_interactions
  )
}
