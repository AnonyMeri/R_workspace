# R_workspace

実験データの前処理、正規性検定、分散分析・順位検定、回帰・相関・SEMを行うRスクリプトです。

## フォルダと説明書

| 場所 | 内容・説明書 |
|---|---|
| `script/` | [スクリプトの構成、導入、必須ライブラリ](script/README.md) |
| `script/preprocess/` | [前処理、記述統計、正規性検定と参考資料](script/preprocess/README.md) |
| `script/test/` | [検定方法、下位検定、補正と関連論文](script/test/README.md) |
| `script/relation/` | [回帰・相関・SEMの使い方と関連資料](script/relation/README.md) |
| `test_data/` | [模擬CSVの一覧、列の意味、再生成方法](test_data/README.md) |
| `examples/` | 分析例と模擬データの生成スクリプト |
| `tests/` | [一式の動作確認方法と確認記録](tests/README.md) |

## 始め方

R/RStudioの作業ディレクトリを、このREADMEがある`R_workspace`に設定します。

```r
source("script/install_packages.R", encoding = "UTF-8") # 初回
source("script/load_script.R", encoding = "UTF-8")
```

GUIでは`experiment <- preprocess(data)`で列を選び、`result <- analyze(experiment)`で分析方法を選びます。引数から指定する例は各分析フォルダの説明を参照してください。

## 実験CSVと模擬データ

`script/`と同じ階層に置いたCSVは、拡張子の大文字小文字によらず`.gitignore`で除外します。`test_data/`の模擬CSVは登録対象です。Rの作業状態、環境変数ファイル、`results/`の保存結果も除外します。

実験CSVを移動・削除する必要はありません。既に追跡されたCSVには`.gitignore`は効かないため、アップロード前に`git status --short`と`git ls-files`で対象を確認します。必要なら該当CSVを`git rm --cached -- ファイル名.csv`で追跡対象から外します。この操作はローカルファイルを削除しませんが、過去のコミットの内容を履歴から消すものでもありません。

## 既存環境への取り込み

共有済みのスクリプトから構成しています。デスクトップの現物との差分は未確認です。既存環境へ取り込む際は、並行して変更した`script/`との差分を確認してください。収録元と互換修正は[script/README.md](script/README.md#収録元と互換修正)に記載しています。

GitHubには`R_workspace`をルートとして登録します。ライセンスはまだ指定していません。
