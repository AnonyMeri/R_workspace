# R_workspace

実験データの前処理、正規性検定、分散分析・順位検定、回帰・相関・SEMを行うRスクリプトです。

`relate()`（回帰・相関・SEM）は現時点では仮実装で、仕様や処理内容について十分な推敲はまだ行っていません。今後の見直しを予定しています。

## フォルダと説明書

| 場所 | 内容・説明書 |
|---|---|
| `script/` | [スクリプトの構成、導入、必須ライブラリ](script/README.md) |
| `script/preprocess/` | [前処理、記述統計、正規性検定と参考資料](script/preprocess/README.md) |
| `script/test/` | [検定方法、下位検定、補正と関連論文](script/test/README.md) |
| `script/relation/` | [回帰・相関・SEMの使い方と関連資料](script/relation/README.md) |
| `test_data/` | [模擬CSVの一覧、列の意味、使用例](test_data/README.md) |
| `LICENSE` | [GNU General Public License v3.0の本文](LICENSE) |

現在のリポジトリには`examples/`・`tests/`は含まれていません。実行例は各フォルダのREADMEに記載しています。`script/test/`は統計検定の実装を置くフォルダです。

## 始め方

R/RStudioの作業ディレクトリを、このREADMEがある`R_workspace`に設定します。

```r
source("script/install_packages.R", encoding = "UTF-8") # 初回
source("script/load_script.R", encoding = "UTF-8")
```

GUIでは`experiment <- preprocess(data)`で列を選び、`result <- analyze(experiment)`で分析方法を選びます。引数から指定する例は各分析フォルダの説明を参照してください。

同梱の模擬データを使い、GUIなしで実行する最小例です。以下はR標準機能だけで計算できます。

```r
data <- read.csv("test_data/between_2.csv")
experiment <- preprocess(data, id = "ID", iv = "group", dv = "score",
  select = "arguments", show = FALSE)
normality_result <- normality(experiment, show = FALSE)
result <- analyze(experiment, method = "t", select = "arguments", show = TRUE)
```

GUIにはRの`tcltk`対応と画面環境が必要です。追加パッケージは分析方法ごとに異なります。[依存パッケージ一覧](script/README.md#必須ライブラリ)を参照してください。

## 実験CSVと模擬データ

`script/`と同じ階層に置いたCSVは、拡張子の大文字小文字によらず`.gitignore`で除外します。`test_data/`の模擬CSVは登録対象です。Rの作業状態、環境変数ファイル、`results/`の保存結果も除外します。

実験CSVを移動・削除する必要はありません。既に追跡されたCSVには`.gitignore`は効かないため、アップロード前に`git status --short`と`git ls-files`で対象を確認します。必要なら該当CSVを`git rm --cached -- ファイル名.csv`で追跡対象から外します。この操作はローカルファイルを削除しませんが、過去のコミットの内容を履歴から消すものでもありません。

## 開発・確認状況

作成中のスクリプト集です。実験計画や入力条件への対応範囲は各READMEを確認してください。確認内容と制約は[script/README.md](script/README.md#実行例と動作確認)にまとめています。別の環境で変更したスクリプトを取り込む際は、ファイルごとの差分を確認してください。

## ライセンス

[GNU General Public License v3.0](LICENSE)。
