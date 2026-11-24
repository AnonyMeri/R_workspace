# ============================================================
# 分析スクリプトの一括読み込み
#
# 推奨配置：
#   script/
#   ├─ load_script.R
#   ├─ preprocess/
#   │  ├─ preprocess_helper.R
#   │  ├─ p_prepare.R
#   │  ├─ p_describe.R
#   │  └─ p_normality.R
#   ├─ test/
#   │  ├─ test_helper.R
#   │  ├─ t_anova.R
#   │  ├─ t_art.R
#   │  ├─ t_ancova.R
#   │  ├─ t_permuco.R
#   │  ├─ t_rank.R
#   │  ├─ t_brunner_munzel.R
#   │  ├─ t_ats.R
#   │  ├─ t_t.R
#   │  ├─ t_welch.R
#   │  └─ t_rank_ancova.R
#   └─ relation/
#      ├─ relation_helper.R
#      ├─ r_regression.R
#      ├─ r_correlation.R
#      └─ r_sem.R
#
# 使い方：
#   source("script/load_script.R", encoding = "UTF-8")
#
# 公開関数：
#   preprocess()
#   normality()
#   analyze()
#   relate()
#
# p_prepare(), p_describe(), p_normality(), t_anova(), t_art(),
# t_ancova(), t_permuco(), t_rank(), r_regression()および各補助関数はprivate environmentに保持し、
# Global Environmentには公開しません。
# ============================================================

local({

  # ----------------------------------------------------------
  # load_script.R自身の場所を取得
  # ----------------------------------------------------------

  source_file <- NULL
  frames <- sys.frames()

  for (index in rev(seq_along(frames))) {
    candidate <- frames[[index]]$ofile

    if (!is.null(candidate) && nzchar(candidate)) {
      source_file <- candidate
      break
    }
  }

  if (is.null(source_file)) {
    stop(
      "load_script.Rの場所を取得できませんでした。\n",
      'コンソールからsource("script/load_script.R", encoding = "UTF-8")を実行してください。',
      call. = FALSE
    )
  }

  script_directory <- dirname(
    normalizePath(
      source_file,
      winslash = "/",
      mustWork = TRUE
    )
  )


  # ----------------------------------------------------------
  # 工程別の統合helper
  # ----------------------------------------------------------

  helper_exports <- list(c("preprocess", "normality"), "analyze", "relate")
  names(helper_exports) <- c(
    file.path("preprocess", "preprocess_helper.R"),
    file.path("test", "test_helper.R"),
    file.path("relation", "relation_helper.R")
  )
  helper_relative_paths <- names(helper_exports)

  helper_paths <- file.path(
    script_directory,
    helper_relative_paths
  )

  missing_helpers <- helper_paths[
    !file.exists(helper_paths)
  ]

  if (length(missing_helpers) > 0) {
    stop(
      "次の統合helperが見つかりません：\n",
      paste0(
        "  - ",
        gsub("\\\\", "/", missing_helpers),
        collapse = "\n"
      ),
      "\n\nload_script.Rをscriptフォルダ直下に置き、",
      "preprocess、test、relationフォルダの配置を確認してください。",
      call. = FALSE
    )
  }


  # ----------------------------------------------------------
  # Global Environmentへ公開する関数
  # ----------------------------------------------------------

  public_functions <- unlist(helper_exports, use.names = FALSE)


  # ----------------------------------------------------------
  # helperごとに専用環境へ読み込み、担当する公開関数を取得する。
  # normality()の提供元はpreprocess_helper.Rとする。
  # relate()の提供元はrelation/relation_helper.Rとする。
  # relation内のr_*はrelation_helperが専用環境へ読み込む。
  # 全helperの検証が完了してから公開する。
  # ----------------------------------------------------------

  exports <- list()
  for (index in seq_along(helper_paths)) {
    module_environment <- new.env(parent = .GlobalEnv)
    source(
      helper_paths[index],
      local = module_environment,
      encoding = "UTF-8",
      chdir = FALSE
    )
    expected <- helper_exports[[index]]
    loaded <- vapply(expected, exists, logical(1), envir = module_environment,
      mode = "function", inherits = FALSE)
    if (any(!loaded)) {
      stop(
        helper_relative_paths[index], "から次の公開関数を読み込めませんでした：\n",
        paste0("  - ", expected[!loaded], "()", collapse = "\n"),
        "\n対応するhelperまたは子スクリプトを確認してください。",
        call. = FALSE
      )
    }
    for (function_name in expected) {
      exports[[function_name]] <- get(function_name, envir = module_environment, inherits = FALSE)
    }
  }
  for (function_name in public_functions) {
    assign(function_name, exports[[function_name]], envir = .GlobalEnv)
  }


  # ----------------------------------------------------------
  # 依存パッケージの確認は各helperへ委譲
  #
  # load_script.Rはhelperの読み込みと公開関数の登録だけを担当する。
  # 正規性検定のnortestは任意の追加検定を実行する際に確認する。
  # test系の依存パッケージはtest_helper.Rが読み込み時に確認し、
  # 各t_*()も実行時に最終確認する。
  # relation系の依存パッケージはrelation_helper.Rへ委譲する。
  # ----------------------------------------------------------


  # ----------------------------------------------------------
  # 読み込み結果を表示
  #
  # options(analysis.load.quiet = TRUE)で非表示にできます。
  # ----------------------------------------------------------

  quiet <- getOption(
    "analysis.load.quiet",
    FALSE
  )

  if (!isTRUE(quiet)) {
    cat(
      "\n========================================\n",
      "分析スクリプトを読み込みました。\n",
      "========================================\n",
      sep = ""
    )

    cat(
      "\n【scriptディレクトリ】\n",
      script_directory,
      "\n",
      sep = ""
    )

    cat(
      "\n【読み込んだ統合helper】\n",
      paste0(
        "  - ",
        gsub("\\\\", "/", helper_relative_paths),
        collapse = "\n"
      ),
      "\n",
      sep = ""
    )

    cat(
      "\n【公開関数】\n",
      paste0(
        "  - ",
        public_functions,
        "()",
        collapse = "\n"
      ),
      "\n",
      sep = ""
    )

  }

}, envir = new.env(parent = .GlobalEnv))
