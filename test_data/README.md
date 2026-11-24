# 模擬データ

すべて`examples/generate_test_data.R`が乱数seed `20261124`から生成した人工データです。実験CSVの読み込み・加工や実際の参加者情報の使用はありません。

| CSV | 計画・行数 | 主な確認対象 |
|---|---|---|
| `between_2.csv` | 対応なし2群、48行 | Student/Welch t、Wilcoxon、Brunner–Munzel |
| `between_3.csv` | 対応なし3群、72行 | ANOVA、ART、permuco、Kruskal–Wallis/Dunn、ATS/MCTP、3方式のANCOVA |
| `unequal_variance_3.csv` | 対応なし3群、72行、群ごとに生成時の誤差SDが異なる | Welch/Games–Howell、ATS/MCTP |
| `within_2.csv` | 48人×2条件、96行 | 対応ありt、Wilcoxon符号付順位 |
| `within_3.csv` | 48人×3条件、144行 | 反復ANOVA、ART、permuco、Friedman/対応ありWilcoxon、球面性 |
| `mixed_2x3.csv` | 2群×3条件、48人、144行 | 混合ANOVA/ANCOVA、ART、permuco、交互作用 |
| `between_2x3.csv` | 対応なし2×3、108人 | 要因ANOVA/ANCOVA、ART、permuco |
| `relation.csv` | 対応なし300人、300行 | 回帰、相関、観測変数のSEM・間接効果 |

検定用CSVの`ID`は模擬参加者、`group`・`condition`・`A`・`B`はカテゴリ要因、`score`は連続従属変数、`covariate`は連続共変量です。同じIDの共変量は一定で、ID×条件につき1行です。列名の接頭辞による自動選択は行いません。

`relation.csv`の`x`、`mediator`、`outcome`、`auxiliary`は連続変数です。生成モデルは`mediator = .6*x + 誤差`、`outcome = .3*x + .7*mediator + 誤差`です。これは実際の因果関係の証拠や分析上の推奨モデルを意味しません。

これらは実行経路を確認するためのデータです。正規性の棄却、球面性補正の種類、有意差の有無を固定するテストや、統計手法の性能評価には使いません。
