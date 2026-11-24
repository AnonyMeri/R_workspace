# 前処理・記述統計・正規性検定

[全体説明](../../README.md)／[構成と必須ライブラリ](../README.md)。

## 前処理

`preprocess()`で列の役割を指定し、後続の分析で使う`analysis_experiment`を作成します。GUIでの列選択はチェックボタンです。列名の`iv_`・`dv_`等の接頭辞による自動選択は行いません。

```r
source("script/load_script.R", encoding = "UTF-8")
data <- read.csv("test_data/mixed_2x3.csv")
experiment <- preprocess(data, id = "ID", iv = c("group", "condition"),
  dv = "score", covariate = "covariate", show = FALSE, select = "arguments")

# GUIで選択する場合
# experiment <- preprocess(data)
```

| 引数 | 指定する列 |
|---|---|
| `id` | 参加者ID、1列 |
| `iv` | カテゴリ要因、1列以上 |
| `dv` | 数値の従属変数、1列以上 |
| `covariate` | 数値の共変量、任意・複数可 |

同じ参加者内で水準が変化する要因を`within`、一定の要因を`between`として判定し、`design`を`within`・`between`・`mixed`に分類します。判定はIDと要因の変化に基づくため、研究上の実験計画と一致するか確認してください。欠けている条件や不正なIDを自動的に補う処理ではありません。

返り値には`data`、`id`、`iv`、`dv`、`covariate`、`within`、`between`、`design`、`factor_levels`、`n_participants`等が入ります。要因はfactorに変換し、共変量は数値として保持します。共変量を実際にモデルへ入れるかは、後続で選択した分析方法が決めます。

## 記述統計

`preprocess()`は記述統計も計算し、`experiment$descriptives`へ保存します。閲覧用の`describe`オブジェクトも作成します。`show = FALSE`は表示を抑える指定です。`ci`は記述統計の信頼水準、`digits`は表示桁数です。

## 正規性検定

`normality()`は前処理と独立した公開関数です。`preprocess_helper.R`の読み込み時に利用可能になりますが、`preprocess()`から自動実行されるわけではありません。

```r
normality_result <- normality(experiment)
normality_all <- normality(experiment,
  methods = c("shapiro", "anderson", "lilliefors"), show = FALSE)
```

| `methods` | 方法 | 実行元 |
|---|---|---|
| `shapiro`（既定） | Shapiro–Wilk | `stats::shapiro.test()` |
| `anderson` | Anderson–Darling | `nortest::ad.test()` |
| `lilliefors` | Lilliefors | `nortest::lillie.test()` |
| `all` | 上記3方法 | 上記の関数 |

従属変数ごと、要因水準の組み合わせごとに検定します。同じID×条件に複数行がある場合は、参加者内の平均にまとめて評価します。標本数不足や一定値などで検定できないセルは、理由を記録して参考判定の分母から除外します。Shapiro–Wilkは有効観測数3〜5000、Anderson–Darlingは8以上、Lillieforsは5以上が対象です。

### 分析方法の参考候補

複数方法を指定した場合は`methods`の先頭を主要検定として使い、「棄却セル数／検定可能セル数」を計算します。

| 主要検定の棄却率 | 共変量なし | 共変量あり |
|---|---|---|
| 半数未満 | ANOVA | ANCOVA |
| 半数以上 | ART | 置換ANCOVA（permuco） |

検定可能なセルがない場合は判定を保留します。これは本スクリプトで決めた参考表示のルールであり、原著論文の推奨基準ではありません。棄却されないことは正規分布の証明ではなく、この表から分析を自動実行・自動変更することもありません。ANOVA・ANCOVAの判断では、モデル推定後の残差、分布、外れ値、標本数も確認します。

## 関連資料

- [R公式：shapiro.test](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/shapiro.test.html)：実装と標本数の条件。原著はShapiro, S. S., & Wilk, M. B. (1965). *An analysis of variance test for normality (complete samples).* Biometrika, 52(3/4), 591–611。
- [nortest公式リファレンス](https://www.stat.ethz.ch/CRAN/web/packages/nortest/refman/nortest.html)：`ad.test()`・`lillie.test()`の条件と参考文献。例えばStephens (1986), *Tests based on EDF statistics*, in *Goodness-of-Fit Techniques*、Dallal & Wilkinson (1986), *An analytic approximation to the distribution of Lilliefors' test for normality*, The American Statistician, 40, 294–296。

これらは各検定の根拠資料です。「半数以上で候補を切り替える」という独自の運用基準を裏づける引用としては扱いません。

## 模擬データと実行例

[test_data/README.md](../../test_data/README.md)に計画別CSVと列の説明があります。このページの前処理例と、[検定の実行例](../test/README.md#6-実行例再読み込み保存結果)を組み合わせて利用できます。
