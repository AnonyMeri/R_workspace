# ============================================================
# 前処理の統合インターフェース
#
# 読み込む子スクリプト：
#   p_prepare.R
#   p_describe.R
#   p_normality.R
#
# 利用者向け公開関数：
#   preprocess()
#   normality()（必要時に独立して実行）
#
# p_prepare(), p_describe(), p_normality() および .p_* 内部関数は
# private environmentへ読み込み、Global Environmentには公開しません。
# preprocess()実行時には、閲覧用のdescribeオブジェクトだけを
# Global Environmentへ自動作成します。
#
# 使い方：
#   source("script/preprocess/preprocess_helper.R", encoding = "UTF-8")
#
#   experiment <- preprocess(data)  # GUIで列を選択
#   experiment <- preprocess(data, select = "gui")
#
#   # 引数による直接指定も従来どおり利用できます。
#   experiment <- preprocess(
#     data = data,
#     id = "ID",
#     iv = c("A", "B"),
#     dv = c("Q1", "Q2"),
#     covariate = "PreScore"
#   )
#
#   normality_result <- normality(experiment)
#
#   View(describe$Q1$overall)
#   View(describe$Q1$A)
#   View(describe$Q1$A_B)
#   cat(experiment$selection$code)  # 同じ指定を再利用するコード
# ============================================================


# ------------------------------------------------------------
# helper自身の場所を取得
# ------------------------------------------------------------

.preprocess_source_file <- NULL
.preprocess_frames <- sys.frames()

for (.preprocess_index in rev(seq_along(.preprocess_frames))) {
  .preprocess_candidate <- .preprocess_frames[[.preprocess_index]]$ofile

  if (!is.null(.preprocess_candidate) && nzchar(.preprocess_candidate)) {
    .preprocess_source_file <- .preprocess_candidate
    break
  }
}

if (is.null(.preprocess_source_file)) {
  stop(
    "preprocess_helper.Rの場所を取得できませんでした。\n",
    'コンソールからsource("script/preprocess/preprocess_helper.R")を実行してください。',
    call. = FALSE
  )
}

.preprocess_directory <- dirname(
  normalizePath(
    .preprocess_source_file,
    winslash = "/",
    mustWork = TRUE
  )
)

.preprocess_child_paths <- file.path(
  .preprocess_directory,
  c("p_prepare.R", "p_describe.R", "p_normality.R")
)

.preprocess_missing_scripts <- .preprocess_child_paths[
  !file.exists(.preprocess_child_paths)
]

if (length(.preprocess_missing_scripts) > 0) {
  stop(
    "次の前処理スクリプトが見つかりません：\n",
    paste0(
      "  - ",
      .preprocess_missing_scripts,
      collapse = "\n"
    ),
    "\n\npreprocess_helper.R、p_prepare.R、p_describe.R、p_normality.Rを",
    "同じフォルダに配置してください。",
    call. = FALSE
  )
}


# ------------------------------------------------------------
# 子スクリプトはprivate environmentへ読み込む
#
# p_prepare.R / p_describe.R / p_normality.R内の補助関数も含めて、
# RStudioのGlobal Environment / Functionsには公開しない。
# ------------------------------------------------------------

.preprocess_public_environment <- environment()
.preprocess_internal_environment <- new.env(
  parent = .preprocess_public_environment
)

for (.preprocess_path in .preprocess_child_paths) {
  source(
    .preprocess_path,
    local = .preprocess_internal_environment,
    encoding = "UTF-8",
    chdir = FALSE
  )
}

.preprocess_required_functions <- c("p_prepare", "p_describe", "p_normality")
.preprocess_loaded <- vapply(
  .preprocess_required_functions,
  exists,
  logical(1),
  envir = .preprocess_internal_environment,
  mode = "function",
  inherits = FALSE
)

if (any(!.preprocess_loaded)) {
  stop(
    "次の前処理関数を読み込めませんでした：\n",
    paste0(
      "  - ",
      .preprocess_required_functions[!.preprocess_loaded],
      "()",
      collapse = "\n"
    ),
    call. = FALSE
  )
}


# ------------------------------------------------------------
# preprocess系の依存パッケージ確認
#
# 現在のp_prepare.R / p_describe.Rはbase R、stats、utilsのみを使用し、
# 基本処理とShapiro-Wilk検定に追加パッケージは不要です。
# Anderson-Darling / Lilliefors検定のみnortestが必要です。
# GUIのみR付属のtcltkを使用します。GUIを開けるデスクトップ環境が必要です。
#
# options(analysis.helper.quiet = TRUE)
# またはoptions(analysis.load.quiet = TRUE)で表示を抑制できる。
# ------------------------------------------------------------

.preprocess_helper_quiet <- isTRUE(getOption("analysis.helper.quiet", FALSE)) ||
  isTRUE(getOption("analysis.load.quiet", FALSE))

if (!.preprocess_helper_quiet) {
  cat(
    "\n========================================\n",
    "preprocess helper：依存パッケージ確認\n",
    "========================================\n",
    "基本処理・Shapiro-Wilk検定：追加パッケージ不要\n",
    "Anderson-Darling / Lilliefors検定：nortest（使用時のみ必要）\n",
    sep = ""
  )
}


# ------------------------------------------------------------
# 公開関数：データ整形と記述統計を一括実行
# ------------------------------------------------------------

preprocess <- local({
  internal <- .preprocess_internal_environment

  # GUIは指定の受付のみを担当。検証とexperiment作成はp_prepareに委譲する。
  choose_experiment <- function(data, initial) {
    internal$.p_validate_prepare_data(data)
    for (role in names(initial)) {
      value <- initial[[role]]
      if (is.null(value)) next
      internal$.p_check_column_names(value, role, allow_empty = TRUE)
      if (length(setdiff(value, names(data)))) {
        stop(role, "に指定した列がデータにありません。", call. = FALSE)
      }
      if (role == "id" && length(value) > 1L) {
        stop("idには参加者ID列を1つ指定してください。", call. = FALSE)
      }
    }
    if (!isTRUE(capabilities("tcltk")) ||
        !suppressWarnings(requireNamespace("tcltk", quietly = TRUE))) {
      stop("GUIにはRのtcltkが必要です。利用できない環境ではid・iv・dvを引数で指定してください。", call. = FALSE)
    }
    window <- tryCatch(tcltk::tktoplevel(), error = function(e) {
      stop("GUIを開けませんでした。デスクトップ版のR/RStudioで実行するか、id・iv・dvを引数で指定してください。", call. = FALSE)
    })
    on.exit({
      try(tcltk::tkgrab.release(window), silent = TRUE)
      try(tcltk::tkdestroy(window), silent = TRUE)
    }, add = TRUE)
    tcltk::tkwm.title(window, "分析する列の指定")
    tcltk::tkwm.minsize(window, 720, 380)
    tcltk::tkgrid.columnconfigure(window, 0, weight = 1)
    tcltk::tkgrid.rowconfigure(window, 1, weight = 1)
    note <- tcltk::tklabel(window,
      text = "チェックを入れて列を選択してください。参加者IDは1列、独立変数・従属変数・共変量は複数列を選択できます。",
      anchor = "w", justify = "left", wraplength = 760)
    tcltk::tkgrid(note, row = 0, column = 0, sticky = "ew", padx = 12, pady = 10)
    body <- tcltk::tkframe(window)
    tcltk::tkgrid(body, row = 1, column = 0, sticky = "nsew", padx = 8)
    tcltk::tkgrid.rowconfigure(body, 0, weight = 1)
    labels <- c(id = "参加者ID（1列）", iv = "独立変数（1列以上）",
      dv = "従属変数（数値・1列以上）", covariate = "共変量（数値・任意）")
    columns <- names(data)
    display_names <- vapply(columns, function(name) {
      paste0(name, "  [", paste(class(data[[name]]), collapse = "/"), "]")
    }, character(1))
    selections <- list()
    first_id_checkbox <- NULL
    for (i in seq_along(labels)) {
      role <- names(labels)[i]
      panel <- tcltk::tkframe(body)
      tcltk::tkgrid(panel, row = 0, column = i - 1L, sticky = "nsew", padx = 4)
      tcltk::tkgrid.columnconfigure(body, i - 1L, weight = 1)
      tcltk::tkgrid.rowconfigure(panel, 1, weight = 1)
      tcltk::tkgrid.columnconfigure(panel, 0, weight = 1)
      tcltk::tkgrid(tcltk::tklabel(panel, text = labels[i], anchor = "w"),
        row = 0, column = 0, columnspan = 2, sticky = "w", pady = 4)
      canvas <- tcltk::tkcanvas(panel, height = 300, width = 210,
        highlightthickness = 1, borderwidth = 0, yscrollincrement = 24)
      checks <- tcltk::tkframe(canvas)
      tcltk::tkcreate(canvas, "window", 0, 0, window = checks, anchor = "nw")
      variables <- setNames(lapply(columns, function(name) {
        tcltk::tclVar(if (name %in% initial[[role]]) "1" else "0")
      }), columns)
      selections[[role]] <- variables
      check_widgets <- list()
      for (j in seq_along(columns)) {
        checkbox <- tcltk::tkcheckbutton(checks, text = display_names[j],
          variable = variables[[j]], anchor = "w")
        if (role == "id") {
          # チェックボックスでもIDの選択は1列に限定する。
          command <- local({
            current_variables <- variables
            current_index <- j
            function() {
              if (tcltk::tclvalue(current_variables[[current_index]]) == "1") {
                for (other in setdiff(seq_along(current_variables), current_index)) {
                  tcltk::tclvalue(current_variables[[other]]) <- "0"
                }
              }
              invisible(NULL)
            }
          })
          tcltk::tkconfigure(checkbox, command = command)
          if (is.null(first_id_checkbox)) first_id_checkbox <- checkbox
        }
        tcltk::tkpack(checkbox, side = "top", anchor = "w", fill = "x")
        check_widgets[[j]] <- checkbox
      }
      vertical <- tcltk::tkscrollbar(panel, orient = "vertical", command = function(...) {})
      horizontal <- tcltk::tkscrollbar(panel, orient = "horizontal", command = function(...) {})
      # localで各欄のキャンバス・スクロールバー・チェック群を固定する。
      local({
        current_canvas <- canvas; current_checks <- checks
        current_v <- vertical; current_h <- horizontal
        tcltk::tkconfigure(current_v, command = function(...) tcltk::tkyview(current_canvas, ...))
        tcltk::tkconfigure(current_h, command = function(...) tcltk::tkxview(current_canvas, ...))
        tcltk::tkconfigure(current_canvas,
          yscrollcommand = function(...) tcltk::tkset(current_v, ...),
          xscrollcommand = function(...) tcltk::tkset(current_h, ...))
        tcltk::tkbind(current_checks, "<Configure>", function() {
          tcltk::tkconfigure(current_canvas, scrollregion = tcltk::tkbbox(current_canvas, "all"))
        })
        for (widget in c(list(current_canvas, current_checks), check_widgets)) {
          tcltk::tkbind(widget, "<Button-4>", function() tcltk::tkyview(current_canvas, "scroll", -1, "units"))
          tcltk::tkbind(widget, "<Button-5>", function() tcltk::tkyview(current_canvas, "scroll", 1, "units"))
          tcltk::tkbind(widget, "<MouseWheel>", function(D) {
            delta <- suppressWarnings(as.numeric(D))
            if (length(delta) == 1L && is.finite(delta) && delta != 0) {
              units <- -sign(delta) * max(1L, floor(abs(delta) / 120))
              tcltk::tkyview(current_canvas, "scroll", units, "units")
            }
            invisible(NULL)
          })
        }
      })
      tcltk::tkgrid(canvas, row = 1, column = 0, sticky = "nsew")
      tcltk::tkgrid(vertical, row = 1, column = 1, sticky = "ns")
      tcltk::tkgrid(horizontal, row = 2, column = 0, sticky = "ew")
    }
    error_text <- tcltk::tclVar("")
    tcltk::tkgrid(tcltk::tklabel(window, textvariable = error_text, foreground = "#a32020",
      anchor = "w", justify = "left", wraplength = 760),
      row = 2, column = 0, sticky = "ew", padx = 12, pady = 6)
    result <- NULL
    confirm <- function() {
      selected <- lapply(selections, function(variables) {
        columns[vapply(variables, function(variable) tcltk::tclvalue(variable) == "1", logical(1))]
      })
      prepared <- tryCatch(do.call(internal$p_prepare, c(list(data = data), selected)),
        error = function(e) {
          tcltk::tclvalue(error_text) <- conditionMessage(e)
          NULL
        })
      if (is.null(prepared)) return(invisible(NULL))
      result <<- prepared
      tcltk::tkdestroy(window)
      invisible(NULL)
    }
    cancel <- function() tcltk::tkdestroy(window)
    buttons <- tcltk::tkframe(window)
    tcltk::tkgrid(buttons, row = 3, column = 0, sticky = "e", padx = 12, pady = 10)
    tcltk::tkpack(tcltk::tkbutton(buttons, text = "キャンセル", command = cancel), side = "left", padx = 4)
    tcltk::tkpack(tcltk::tkbutton(buttons, text = "確定", command = confirm, default = "active"), side = "left", padx = 4)
    tcltk::tkwm.protocol(window, "WM_DELETE_WINDOW", cancel)
    tcltk::tkbind(window, "<Escape>", cancel)
    tcltk::tkbind(window, "<Return>", confirm)
    tcltk::tkfocus(first_id_checkbox)
    if (!as.logical(as.integer(tcltk::tkwinfo("viewable", window)))) {
      tcltk::tkwait.visibility(window)
    }
    tcltk::tkgrab.set(window)
    tcltk::tkwait.window(window)
    if (is.null(result)) stop("列の選択をキャンセルしました。experimentは作成していません。", call. = FALSE)
    result
  }

  function(
    data,
    id = NULL,
    iv = NULL,
    dv = NULL,
    covariate = NULL,
    ci = 0.95,
    digits = 3,
    show = TRUE,
    select = if (is.null(id) || is.null(iv) || is.null(dv)) "gui" else "arguments"
  ) {
    data_expression <- substitute(data)
    select <- match.arg(select, c("gui", "arguments"))
    # GUIを開く前に、指定方法に依存しない表示引数を確認する。
    if (!is.numeric(ci) || length(ci) != 1L || !is.finite(ci) || ci <= 0 || ci >= 1) {
      stop("ciには0より大きく1より小さい数値を指定してください。", call. = FALSE)
    }
    if (!is.numeric(digits) || length(digits) != 1L || !is.finite(digits) || digits < 0 || digits != floor(digits)) {
      stop("digitsには0以上の整数を指定してください。", call. = FALSE)
    }
    if (!is.logical(show) || length(show) != 1L || is.na(show)) {
      stop("showにはTRUEまたはFALSEを指定してください。", call. = FALSE)
    }
    roles <- list(id = id, iv = iv, dv = dv, covariate = covariate)
    experiment <- if (select == "gui") {
      choose_experiment(data, roles)
    } else {
      do.call(internal$p_prepare, c(list(data = data), roles))
    }
    roles <- experiment[c("id", "iv", "dv", "covariate")]
    if (!length(roles$covariate)) roles["covariate"] <- list(NULL)
    replay_call <- as.call(c(list(as.name("preprocess")), list(data = data_expression),
      roles, list(ci = ci, digits = digits, show = show, select = "arguments")))
    code <- deparse(replay_call, width.cutoff = 80L)
    code[1L] <- paste0("experiment <- ", code[1L])
    experiment$selection <- list(method = select, call = replay_call,
      code = paste(code, collapse = "\n"))

    descriptives <- internal$p_describe(
      experiment = experiment,
      ci = ci,
      digits = digits,
      show = show
    )

    experiment$descriptives <- descriptives

    # 閲覧専用の単純構造を、experimentとは別オブジェクトとして出力する。
    # Global Environmentにはdescribeという1オブジェクトだけを追加する。
    describe <- internal$.p_make_describe_object(
      results = descriptives,
      iv = experiment$iv,
      digits = digits
    )
    assign("describe", describe, envir = .GlobalEnv)

    attr(experiment, "stage") <- "preprocess"
    attr(experiment, "descriptive_ci") <- ci
    class(experiment) <- unique(c(
      "analysis_preprocessed",
      class(experiment)
    ))

    if (show && select == "gui") {
      cat("\n【再実行用コード】\n", experiment$selection$code, "\n", sep = "")
    }
    invisible(experiment)
  }
})


# ------------------------------------------------------------
# 公開関数：正規性検定（preprocessから自動実行しない）
# ------------------------------------------------------------

normality <- local({
  internal <- .preprocess_internal_environment

  function(
    experiment,
    methods = "shapiro",
    alpha = 0.05,
    digits = 3,
    p_digits = 3,
    show = TRUE
  ) {
    results <- internal$p_normality(
      experiment = experiment,
      methods = methods,
      alpha = alpha,
      digits = digits,
      p_digits = p_digits,
      show = show
    )
    attr(results, "interface") <- "normality"
    invisible(results)
  }
})


# ------------------------------------------------------------
# helperの初期化にだけ使用したオブジェクトを公開環境から削除
#
# private environmentはpreprocess()/normality()のクロージャ内に保持するため、
# 削除後も内部関数を呼び出せる。
# ------------------------------------------------------------

rm(
  list = c(
    ".preprocess_source_file",
    ".preprocess_frames",
    ".preprocess_index",
    ".preprocess_candidate",
    ".preprocess_directory",
    ".preprocess_child_paths",
    ".preprocess_missing_scripts",
    ".preprocess_public_environment",
    ".preprocess_internal_environment",
    ".preprocess_path",
    ".preprocess_required_functions",
    ".preprocess_loaded",
    ".preprocess_helper_quiet"
  ),
  envir = environment()
)
