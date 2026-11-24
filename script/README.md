# スクリプトの構成と実行環境

[全体説明](../README.md)に戻る。

## 導入・読み込み・再読み込み

作業ディレクトリは`script/`の親、`R_workspace`に設定します。

```r
source("script/install_packages.R", encoding = "UTF-8") # 初回のパッケージ導入
source("script/load_script.R", encoding = "UTF-8")       # 読み込み・再読み込み
```

`load_script.R`は3つのhelperを専用環境へ読み込み、必要な公開関数が揃ったことを確認してから登録します。パッケージを自動インストールする処理ではありません。

| 公開関数 | 提供元 | 詳細 |
|---|---|---|
| `preprocess()` | `preprocess/preprocess_helper.R` | [データ構造の作成と記述統計](preprocess/README.md) |
| `normality()` | 同上 | [独立して実行する正規性検定](preprocess/README.md#正規性検定) |
| `analyze()` | `test/test_helper.R` | [検定・分散分析・共変量の分析](test/README.md) |
| `relate()` | `relation/relation_helper.R` | [回帰・相関・SEM](relation/README.md) |

各helperがGUI、入力の検証、表示、結果の取りまとめを担当し、子スクリプトが計算を担当します。子の計算関数はグローバル環境へ公開しません。`preprocess()`は記述統計を閲覧するための`describe`オブジェクトも作成します。

## ファイル構成

| フォルダ | ファイル | 担当 |
|---|---|---|
| 直下 | `load_script.R` | 全helperの読み込みと4つの公開関数の登録 |
| 直下 | `install_packages.R` | 未導入の依存パッケージの導入 |
| `preprocess/` | `preprocess_helper.R` | 前処理のGUI、公開関数、表示の調整 |
| 同上 | `p_prepare.R` | 列の検証、要因の役割と実験計画の判定、experiment作成 |
| 同上 | `p_describe.R` | 記述統計と閲覧用オブジェクトの作成 |
| 同上 | `p_normality.R` | 条件別の正規性検定、集計、参考候補の作成 |
| `test/` | `test_helper.R` | 方法の選択、共通補正、段階表示、結果の統合 |
| 同上 | `t_t.R` | Student・対応ありt検定 |
| 同上 | `t_anova.R`、`t_ancova.R` | 通常ANOVA・ANCOVA |
| 同上 | `t_welch.R` | WelchとGames–Howell |
| 同上 | `t_art.R` | ART・ART-C |
| 同上 | `t_permuco.R` | 置換ANOVA・ANCOVA |
| 同上 | `t_rank.R` | 対応・水準数に応じた一般的な順位検定 |
| 同上 | `t_brunner_munzel.R` | Brunner–Munzel、ATSと共用する入力検証 |
| 同上 | `t_ats.R` | ATS・MCTP |
| 同上 | `t_rank_ancova.R` | Quadeの順位ANCOVA |
| `relation/` | `relation_helper.R` | 分析の選択、GUI、表示、補正、結果入力の制御 |
| 同上 | `r_regression.R` | 線形回帰、モデルからの傾きの推定、置換回帰 |
| 同上 | `r_correlation.R` | Pearson・Spearman・Kendallの相関 |
| 同上 | `r_sem.R` | lavaanによるSEMと結果の抽出 |

前処理の実フォルダ名は`preprocess/`です。フォルダや子スクリプトを改名する場合は、loaderとhelperの参照先も合わせて変更する必要があります。

## 必須ライブラリ

必要なパッケージは選択した分析によって異なります。下表は直接使用するパッケージです。依存先は通常`install.packages()`が合わせて導入します。

| 機能・方法 | 必須パッケージ | 補足 |
|---|---|---|
| 前処理、記述統計、Shapiro–Wilk、t検定、一般ランクの全体検定、基本回帰、相関 | R標準の`stats`・`utils`等 | 追加導入不要 |
| Anderson–Darling・Lilliefors | `nortest` | 追加の正規性検定を選ぶ場合 |
| ANOVA・ANCOVA | `afex`、`emmeans` | 球面性関連には依存先の`car`も使用 |
| ART・ART-C | `ARTool`、`emmeans`、`car` | `car`は被験者内の球面性評価で使用 |
| 置換ANOVA・ANCOVA | `permuco` | 下位検定は実装していない |
| Dunn、Games–Howell | `PMCMRplus` | 一般ランク・Welchの3群以上の下位検定 |
| Brunner–Munzel | `lawstat` | t近似 |
| ATS＋MCTP | `rankFD` | 全体検定と多重比較 |
| Quadeの順位ANCOVA | `npANCOVA` | 対応なし1要因 |
| モデルからの回帰・傾き比較 | `emmeans` | afexモデルでは`afex`、置換モデルでは`permuco`も使用 |
| SEM | `lavaan` | 推定・CI・適合度の抽出 |
| GUI | Rの`tcltk`対応 | Tcl/Tkと画面環境が必要 |

全方式を試す場合は`install_packages.R`を利用します。特定の方式だけ使う場合は、必要なパッケージだけを導入できます。読み込み時の確認に加え、各方式の実行時にも必要なパッケージを確認します。

```r
install.packages(c("afex", "emmeans", "car")) # ANOVA・ANCOVAの例
```

## 実行例と動作確認

[検定の実行例](../examples/analysis_examples.R)、[回帰・相関・SEMの実行例](../examples/relation_examples.R)を用意しています。動作確認の方法・結果・制約は[tests/README.md](../tests/README.md)にまとめています。

## 確認環境のパッケージ版

以下は動作確認に使用した版の記録です。最低バージョンや固定を要求する設定ではありません。

R version 4.2.2 Patched (2022-11-10 r83330)

- afex: 1.2.1
- emmeans: 1.8.4.1
- car: 3.1.1
- ARTool: 0.11.2
- permuco: 1.1.3
- PMCMRplus: 1.9.12
- lawstat: 3.6
- rankFD: 0.1.1
- npANCOVA: 0.2.0
- lavaan: 0.6.14

## 収録元と互換修正

共有ファイルの最新版を2026-11-24時点で組み合わせた構成です。デスクトップの現物は直接参照していません。

| 収録元の名前 | バージョン |
|---|---|
| `preprocess_helper(2).R` | 2 |
| `relation_helper.R` | 11 |
| `r_correlation.R` | 初版 |
| `r_regression.R` | 2 |
| `load_script(1).R` | 2 |
| `p_prepare.R` | 1 |
| `p_describe.R` | 初版 |
| `p_normality.R` | 2 |
| `r_sem.R` | 初版 |
| `test_script.zip` | 4 |

サーバー上の名前に付いた連番を除き、配置時のファイル名はhelperが参照する名前に統一しています。loaderはファイル一覧のコメントを更新しました。SEMは`lavaan::sem()`が呼び出し名を検査する版でも動作するよう、`do.call()`を名前指定に変更しました。それ以外の処理内容は収録元を維持しています（末尾の空白のみ整理）。

回帰・SEMを別途変更している場合は、デスクトップの最新版と収録ファイルを照合してから取り込んでください。
