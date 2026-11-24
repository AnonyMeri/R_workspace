# R_workspaceから source("script/install_packages.R") で実行。
# 使用する方式の依存パッケージだけ導入する場合は、この一覧から選ぶ。
packages <- c("afex", "emmeans", "car", "ARTool", "permuco", "PMCMRplus",
  "lawstat", "rankFD", "npANCOVA", "nortest", "lavaan")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
still_missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_missing)) stop("導入できていないパッケージ：", paste(still_missing, collapse = ", "))
message("必要なパッケージを確認しました。GUIにはRのtcltk対応も必要です。")
