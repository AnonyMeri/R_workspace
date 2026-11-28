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

置換ANOVA・ANCOVAでは、全体検定に加えて水準間の比較を行えます。比較方法や実行条件は[検定の説明](script/test/README.md)をご覧ください。

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

## データと分析結果

入力するCSVには、参加者ID、実験要因、測定値の列を用意します。共変量を使う場合は、その列も指定してください。列の例は[模擬データの説明](test_data/README.md)にあります。

分析結果は `saveRDS()` で保存できます。Gitを使う場合、ルートに置いたCSVや `results/` は新規の追跡対象から除外されます。すでに追跡済みのデータは自動では除外されないため、共有前に確認してください。

## ライセンス

[GNU General Public License v3.0](LICENSE)。
