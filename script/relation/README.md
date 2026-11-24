# 回帰・相関・SEM

[全体説明](../../README.md)／[構成と必須ライブラリ](../README.md)。公開関数は`relate()`です。

## 入力と分析の選択

| 分析 | 指定 | 入力 |
|---|---|---|
| 回帰 | `analysis = "regression"` | データと式、対応する推定済みモデル、または対応するANCOVA結果 |
| 相関 | `analysis = "correlation"` | データと2つ以上の数値変数 |
| SEM | `analysis = "sem"` | 個票データとlavaan形式のモデル式 |

式や変数の指定が不足する場合はGUIを使います。再現する実行例では`select = "arguments"`を明示します。`show = FALSE`で表示を抑え、返り値をオブジェクトに保存できます。

## 線形回帰と傾きの比較

```r
source("script/load_script.R", encoding = "UTF-8")
data <- read.csv("test_data/relation.csv")
regression_result <- relate(data, analysis = "regression",
  formula = outcome ~ x + mediator, show = FALSE, select = "arguments")
```

通常のデータ入力では線形モデルを推定します。`trend`・`by`・`compare_slopes`は傾きや条件別の傾きの比較に使う引数です。入力モデルの構造によって利用できる推定・比較が異なります。

`analyze()`のANCOVA・置換ANCOVAの結果を渡す経路では、元のモデル、対象効果、共変量、設定を確認して後続の回帰を組み立てます。両方式の結果が混在する場合は、一方の結果を指定します。

**順位ANCOVA（Quade）の結果入力は拒否します。** この方式では群×共変量の交互作用を検定しないため、同じように条件別回帰へ進める扱いにはしません。元のデータを渡して別途回帰を指定する経路は利用できます。

置換モデルからの傾き・傾き差は`permuco`を利用する処理です。通常の`lm`のp値と同じ推定を表示したものではありません。

## 相関分析

```r
correlation_result <- relate(data, analysis = "correlation",
  variables = c("x", "mediator", "outcome"), method = "pearson",
  adjust = "holm", show = FALSE, select = "arguments")
```

`method`は`pearson`・`spearman`・`kendall`から選びます。選んだ変数の全ペアを`stats::cor.test()`で検定します。独立した観測を行とする通常の二変量相関であり、反復測定相関や偏相関は実装していません。

`use = "pairwise.complete.obs"`はペアごと、`"complete.obs"`は選択した全変数に共通する完全ケースで欠測を除外します。CIは`cor.test()`が返す場合だけ保持します。`exact`は`NULL`・`TRUE`・`FALSE`、`alternative`は両側・正・負の片側を指定できます。p値の補正はhelperが担当します。補正したp値によってCIまで同時信頼区間になるわけではありません。

## SEM

```r
sem_result <- relate(data, analysis = "sem", sem_model = '
  mediator ~ a*x
  outcome ~ c*x + b*mediator
  indirect := a*b
  total := c + a*b
', estimator = "ML", show = FALSE, select = "arguments")
```

単一集団の個票データを使うパス解析・測定モデル・構造モデルを`lavaan::sem()`に委譲します。`estimator`は`auto`・`ML`・`MLR`・`WLSMV`を選べます。`auto`は順序変数があればWLSMV、それ以外はMLです。`ordered`に順序尺度の列名を指定でき、ordered factorの列も検出します。

欠測は既定の`missing_data = "listwise"`で行を除外し、`"fiml"`を選ぶとFIMLを利用します。WLSMVとFIMLは併用しません。この版のブートストラップ標準誤差はMLの場合に使います。`std_lv`・`fixed_x`・`meanstructure`等も指定できます。

返り値には`model`、`parameter_estimates`、`fit_measures`、`fit_table`、`r2_table`、`converged`、`admissible`、診断情報等を保存します。収束と不適解の確認は分けて扱います。適合度の抽出に失敗した場合も診断情報へ記録します。

## 補正と表示の設定

`alpha`・`followup_alpha`・`adjust`・`ci`・`digits`・`p_digits`等をhelperで検証します。`test/`とは別の公開インターフェースであり、`relate()`には`trend_alpha`を追加していません。検定側の引数名をそのまま流用しないでください。

## 関連論文・公式資料

- [R公式：lm](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/lm.html)：通常の線形モデルの実装。
- [R公式：cor.test](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cor.test.html)：Pearson・Spearman・Kendallの検定、CI、exact指定の条件。
- [emmeans公式リファレンス](https://search.r-project.org/CRAN/refmans/emmeans/html/emtrends.html)：モデルに基づく傾きの推定。モデルからの条件別回帰・傾き比較で使用。
- Frossard, J., & Renaud, O. (2021). *Permutation Tests for Regression, ANOVA, and Comparison of Signals: The permuco Package.* Journal of Statistical Software, 99(15), 1–32. [DOI: 10.18637/jss.v099.i15](https://doi.org/10.18637/jss.v099.i15)。置換モデルの一次資料。各種の置換法に条件があり、本スクリプト独自の後続処理すべてを論文が直接保証するという位置づけではありません。
- Rosseel, Y. (2012). *lavaan: An R Package for Structural Equation Modeling.* Journal of Statistical Software, 48(2), 1–36. [DOI: 10.18637/jss.v048.i02](https://doi.org/10.18637/jss.v048.i02)。SEMに使用する実装の論文。
- [lavaan公式チュートリアル](https://lavaan.ugent.be/tutorial/index.html)、[推定・欠測・標準誤差](https://lavaan.ugent.be/tutorial/est.html)、[カテゴリ変数](https://lavaan.ugent.be/tutorial/cat.html)：推定法とデータ型の対応。本スクリプトが公開する選択肢は上記に限定します。

## 模擬データ

[test_data/relation.csv](../../test_data/relation.csv)は300行の模擬データです。このページの回帰・相関・SEMのコード例で使用します。列の意味は[test_data/README.md](../../test_data/README.md)を参照してください。
