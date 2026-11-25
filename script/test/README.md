# 検定・分析スクリプトの仕様と参考文献

更新日：2026-11-24

## 1. 対象と構成

この文書は、同梱の `test_helper.R` と10本の `t_*.R` の実装を対象とします。利用者の入口は `analyze()` です。各 `t_*.R` は計算を担当し、方法の選択、球面性補正の適用、多重比較補正、表示、結果の取りまとめはhelperが担当します。各計算法が持つ固有の補正は、そのパッケージに委譲します。

この文書は`script/test/`の分析方法と共通処理を説明します。[全体説明](../../README.md)、[実行例](#6-実行例再読み込み保存結果)、[模擬データ](../../test_data/README.md)も参照してください。

ファイル別の担当と必須ライブラリは[script/README.md](../README.md)にまとめています。

## 2. 方法の一覧と対応範囲

### GUIの並び順

| 1列目 | 2列目 | 3列目 |
|---|---|---|
| t検定 | ANOVA | — |
| Welch | ART（下位検定：ART-C） | permuco |
| ランク | Brunner–Munzel | ATS＋MCTP |

実行可能な方法をチェックボタンで複数選択できます。共変量の分析は、下段で「実行しない／ANCOVA／置換ANCOVA／順位ANCOVA」を縦1列に並べた単一選択です。共変量なしの方法とは併用できますが、共変量を扱う3方法同士の同時指定は、引数からでも停止します。

「主効果に基づく単純比較の追加」は「表示設定」にあります。この項目は表示だけでなく、計算する下位検定も増やす設定です。基本設定で下位検定の実行可否と計算範囲を、表示設定で下位検定の表示範囲を指定します。

### 実装済みの分析

| `method` | 対象 | 全体の検定 | 下位検定 | 共変量 |
|---|---|---|---|---|
| `t` | 1要因・2水準、対応あり／なし | 対応ありt、または等分散を仮定したStudentのt | なし | 使用しない |
| `anova` | 被験者間・被験者内・混合の要因計画 | afexによるANOVA | emmeansによる水準間比較・2要因の単純比較 | 使用しない |
| `rank` | 1要因、対応あり／なし、2水準以上 | 下表で自動選択 | Dunnまたは対応ありWilcoxon | 使用しない |
| `welch` | 対応なし・1要因・2群以上 | 2群：Welchのt、3群以上：Welchの一元配置ANOVA | 3群以上：Games–Howell | 使用しない |
| `art` | 被験者間・被験者内・混合の要因計画 | ART後のANOVA | ART-Cを用いた比較 | 使用しない |
| `permuco` | 被験者間・被験者内・混合の要因計画 | 置換ANOVA | 実装しない | 使用しない |
| `brunner_munzel` | 対応なし・1要因・2群 | Brunner–Munzel、t近似、両側 | なし | 使用しない |
| `ats` | 対応なし・1要因・3群以上 | ATS、F近似 | 全ペアMCTP、Fisher近似 | 使用しない |
| `ancova` | 共変量のある要因計画 | 通常ANCOVA | 共変量の中心値でのemmeans比較 | 使用する |
| `permuco_ancova` | 共変量のある要因計画 | 置換ANCOVA | 実装しない | 使用する |
| `rank_ancova` | 対応なし・1要因・2群以上 | Quadeの順位ANCOVA | 3群以上：各2群で再推定したQuade法 | 使用する |

`nonparam_anova` は `c("art", "permuco")` の別名です。通常ANOVA/ANCOVAの検定タイプは既定III、選択可能なのはII/IIIです。反復測定ARTは `Error()` を使うためType I、対応なしARTはType IIIです。同じ効果名でも、検定タイプ・順位変換・モデルが異なるので、F値やp値の一致を前提にしません。

ANOVA/ARTの全体モデルには複数要因を指定できますが、このスクリプトの下位検定は主効果と2要因交互作用を対象とします。3要因以上の交互作用の下位検定は自動展開しません。

### 「ランク」の内部選択

| 対応 | 水準数 | 検定 | 下位検定 |
|---|---|---|---|
| なし | 2 | Wilcoxon順位和（Mann–Whitney） | なし |
| あり | 2 | Wilcoxon符号付順位 | なし |
| なし | 3以上 | Kruskal–Wallis | Dunnの全ペア |
| あり | 3以上 | Friedman | IDで対応させたWilcoxon符号付順位の全ペア |

構造による選択です。正規性や等分散性の検定結果によって、別系統の方法に自動で切り替える処理ではありません。Wilcoxonは `exact = NULL, correct = TRUE` とし、同順位・ゼロ差などに関するRの警告も結果に保存します。

### 方法を選ぶ際の区別

- t・ANOVA・Welch・通常ANCOVAは、元の数値尺度の平均に関するモデルです。対応なしの不等分散と、反復測定における球面性は異なる仮定です。GG/HFは対応なし群間の不等分散の代わりにはなりません。
- ARTは、効果ごとの整列・順位変換を用いて要因効果を検定します。元尺度の平均差と同じ量を検定するわけではありません。単純な順位変換をそのまま多要因ANOVAに入れる方法とも異なります。
- permucoは、指定したモデル・置換法に基づく帰無分布からp値を求めます。「データが非正規ならあらゆる置換が許される」という意味ではありません。被験者内の対応や共変量を考慮する置換法の条件が必要です。
- 一般的な順位検定も無条件に仮定がないわけではありません。例えば、順位和を位置差として解釈する場合には分布形状の条件が必要です。Brunner–Munzelは2群の相対効果を、ATS＋MCTPは指定した非パラメトリック相対効果を対象にし、不等分散・異なる分布形状を許す系統です。
- Quade法はpermucoの別名ではありません。順位で共変量を調整する1要因の方法です。共変量の周辺分布が群間で共通という条件などがあり、群×共変量の交互作用は扱いません。

## 3. 共通の球面性補正

### 設定は1つ

GUIの「球面性補正」、引数の `correction` を、ANOVA・ANCOVA・ARTで共用します。数値計算によるεやMauchly検定の取得は各エンジン、適用判断・自由度とp値の補正・表示はhelperの `apply_sphericity()` が担当します。

| 指定 | 対象となる効果の処理 |
|---|---|
| `auto`（既定） | Mauchlyの `p >= alpha` なら未補正。`p < alpha` ならGGの推定εを見て、`.75` 以下はGG、`.75` を超えればHF |
| `GG` | Mauchlyの判定にかかわらずGGを適用 |
| `HF` | Mauchlyの判定にかかわらずHFを適用 |
| `none` | 未補正 |

対象は、被験者内要因を含み、その被験者内部分の対比自由度が1を超える効果です。被験者内の対比が1つの効果や純粋な被験者間効果には適用しません。3群×2時点の交互作用は分子自由度が2でも、被験者内の対比は1つなので対象外です。「2水準要因を含む」だけで交互作用全体を対象外にはしません。

`auto` でMauchlyのp値を取得できない場合は、推定できるGGのεで同じGG/HF選択を行い、その旨を通知・保存します。選択基準となるGGのεがない場合、または選択されたGG/HFのεがない場合は、補正後の自由度・p値をNAにし、「算出不可」と表示します。例えばGGのε > .75でもHFのεを取得できない場合は、別の補正に黙って置き換えません。未補正pを補正済みpとして流用しません。HFのεが1を超える場合は上限1を使います。

補正は `df1 × ε`、`df2 × ε` とし、F値は保ったままF分布の上側確率を再計算します。元の自由度とp値も保存します。Mauchly非棄却は球面性が成立する証明ではありません。εの`.75`を目安にGG/HFを選ぶ方法は、[MedCalcの公式資料](https://www.medcalc.org/en/manual/repeated-measures-anova.php)でもGirden (1992) に基づく目安として説明されています。境界`.75`は本実装ではGG側に含めます。Mauchly判定とこの目安を組み合わせる自動分岐は、本スクリプトで選んだ運用方針です。常にGGを使う分析計画なら明示的に `correction = "GG"` としてください。

### ARTについての境界

ARTのMauchly検定とεは、各効果の整列順位データから推定します。ANOVA/ANCOVAとARTは同じ適用ルールを使いますが、評価データが異なるため、同じ効果でも自動判定やεが一致するとは限りません。

**ARTへのGG/HF適用は、このスクリプト独自の拡張です。ART原論文やARTool標準実装が、この拡張の妥当性まで保証したという位置づけにはしません。** GUIと結果にも拡張であることを表示します。

permuco/置換ANCOVAの置換pにはGG/HFを重ねません。t・Welch・一般ランク・Brunner–Munzel・ATS・Quade法にも、この球面性補正は適用しません。GUIでは対象方法と被験者内要因がある場合にだけ設定欄を有効にします。

### 旧引数との互換性

ART専用だった `sphericity` は既定を `NULL` とし、通常は使用しません。引数で明示した場合のみ、ARTについて `correction` より優先します。GUIでは旧値を共通欄の初期値に移し、隠れたART専用指定を残しません。

ANOVAにGG、ARTに未補正を明示する例は次のとおりです。

```r
result <- analyze(experiment, method = c("anova", "art"),
                  correction = "GG", sphericity = "none")
```

## 4. 多重比較と信頼区間

**多重比較には「自動」を追加していません。** 既定は `adjust = "bonferroni"` のままです。`auto` は球面性補正の適用判断にだけ使います。

| 方法 | 補正の計算元 | 同時に補正する比較の範囲（family） |
|---|---|---|
| ANOVA/ANCOVA・主効果 | helper→emmeans | その主効果の全ペア |
| ANOVA/ANCOVA・2要因の単純比較 | helper→emmeans | 比較方向ごと、固定した他方の水準ごとに別family |
| ART-C・主効果 | helper→emmeans | その効果の全ペア |
| ART-C・交互作用 | helper→emmeans | その効果で抽出した単純比較をまとめて1family。両方向を抽出した場合は両方向を含む |
| ランク・3水準以上 | helper→`stats::p.adjust()` | その要因の全ペア。検定不能なペアも総比較数に含める |
| Quade・3群以上 | helper→`stats::p.adjust()` | 各2群で再推定した比較を全ペアまとめて補正。検定不能なペアも総比較数に含める |
| Welch・3群以上 | PMCMRplus／Games–Howell | 全ペアを考慮した固有の調整。追加補正しない |
| ATS後のMCTP | rankFD | 全ペアTukey型対比、Fisher近似の固有の調整。追加補正しない |
| t・2群Welch・Brunner–Munzel | 各検定 | 1比較なので追加補正なし |
| permuco／置換ANCOVA | — | この実装に下位検定はない |

**同じ `adjust` 名でも、ANOVAとARTでfamilyは完全には同じではありません。** 比較群の定義は維持しています。計算範囲と表示範囲の指定では、各効果内のfamilyを変更しません。familyは結果の `family` / `family_size` に保存します。異なる方法・従属変数・効果をまたぐ実験全体の多重性を一括で制御する設定ではありません。

共通補正の対応表をhelperの `adjustment_policy()` に集約し、GUIの選択肢と引数指定時の検証で同じ表を使います。既定はBonferroniのままです。

| 共通補正を使う方法 | 本スクリプトが受け付ける `adjust` |
|---|---|
| ANOVA・ANCOVA・ART-C | `bonferroni`, `holm`, `sidak`, `tukey`, `hochberg`, `hommel`, `BH`, `BY`, `none` |
| 一般ランク・順位ANCOVA | `bonferroni`, `holm`, `hochberg`, `hommel`, `BH`, `BY`, `none` |

`fdr` は `BH` の別名として受け付け、結果には `BH` として記録します。複数の方法を選んだ場合は、共通補正を使う方法の対応集合の共通部分だけをGUIに出します。引数指定でも同じ対応表で検証し、非対応の指定はモデル計算前に停止します。一般ランク／Quadeとの併用時にTukeyやSidakを指定することはできません。

これはこのラッパーの公開選択肢であり、emmeansが提供するすべての高度な補正を公開するという仕様ではありません。以前の、任意の補正文字列をGUIに追加して後段へ渡す処理は廃止しました。

Holm・BonferroniはFWERを制御します。Hochberg・Hommelには検定間の独立性／適切な正の依存の条件があり、BH・BYはFWERではなくFDRを制御する方法です。利用可能な名前が増えても、全方式を同じ目的・同じ条件で選べるという意味にはしません（[R公式p.adjust](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html)）。

Tukeyが比較の構成に適さない場合は、emmeansがSidakなどに置き換えることがあります。要求した `requested_adjust` と、実際の `p_adjustment` / `ci_adjustment` を分けて保存します。Games–Howell/MCTPと併用した場合、共通の `adjust` はそれらの固有補正を置き換えません。固有補正は `additional_adjustment = "none"` として保存し、追加補正しません。

### 下位検定の実行と表示

数値基準は次の2項目にまとめます。独立した「下位検定基準」の入力欄は設けません。

| 引数 | GUIの項目 | 既定 | 用途 |
|---|---|---|---|
| `alpha` | 有意水準 | `.05` | 全体検定・補正後の比較を有意と判定する上限 |
| `trend_alpha` | 有意傾向の上限 | `.10` | `alpha <= p < trend_alpha` を有意傾向として一覧に示す。「有意＋有意傾向」の計算対象を選ぶ基準も共用 |

条件は `0 < alpha <= trend_alpha < 1` です。例えば `.05` / `.10` なら、「有意のみ」は全体検定の `p < .05`、「有意＋有意傾向」は `p < .10`、「すべて」は有意性にかかわらず計算します。

旧引数 `followup_alpha` は、既存コードとの互換用に残します。単独で指定した場合は、同じ値を有意傾向と下位検定の共通基準として扱います。`trend_alpha` と両方を指定する場合は同値が必要で、異なる値なら計算前に停止します。独立した基準として使うことはできません。GUIの再実行コードには `trend_alpha` だけを記載します。

基本設定の「下位検定を実行する」は `posthoc` に対応し、既定は `TRUE` です。`FALSE` なら下位検定の遅延計算を呼び出さず、`contrasts` は空になります。全体検定・球面性評価・共変量の検定は実行します。GUIでは補正法・計算範囲・表示範囲・単純比較の追加を無効にします。有意傾向の判定は全体検定でも使うため、`trend_alpha` は有効のままです。

#### 計算範囲：基本設定

| `followup_scope` | GUIの表示 | 計算する下位検定 |
|---|---|---|
| `significant` | 有意な効果のみ | 全体検定で `p < alpha` の効果 |
| `followup`（既定） | 有意＋有意傾向の効果 | 全体検定で `p < trend_alpha` の効果 |
| `all` | すべて | 全体検定の有意性にかかわらず、その方法が下位検定に対応する全効果 |

計算対象は主効果・交互作用という効果単位です。ANCOVAの共変量や共変量との交互作用を、群の水準間比較へ自動展開しません。方法自体に下位検定がない場合や、2群の検定で下位検定を重ねない場合にも、新しい比較を追加しません。全体検定が算出不能な効果も自動的に別の検定へ切り替えません。

`simple_effects = TRUE` では、2要因交互作用自体が選択した計算基準外でも、構成する主効果がその計算基準内なら、該当方向の単純比較を追加します。`all` では対象効果をすでに計算するため、GUIの追加設定を無効にします。`posthoc = FALSE` なら、この追加設定を指定していても下位検定は行いません。

#### 表示範囲：表示設定

| `display_scope` | GUIの表示 | 表示する下位検定の結果 |
|---|---|---|
| `significant`（既定） | 有意な比較のみ | 各比較の補正後 `p < alpha` |
| `effect` | 有意な効果の全比較 | 全体検定で有意だった効果について、有意・非有意を含む全比較 |
| `all` | 計算済みの全結果 | 全体検定・個別比較の有意性にかかわらず、計算した全比較。算出不能なペアも残す |

表示側が計算範囲を広げることはありません。例えば `followup_scope = "followup"` と `display_scope = "all"` の組合せは、基準内の効果について計算した結果だけをすべて表示します。非有意な効果にも下位検定を行って表示したい場合は、計算側も `all` にします。

複数方法では表示条件を方法ごとに判定します。一方の方法で有意でも、別の方法の非有意な比較を「有意な比較のみ」の表に載せることはありません。表示行の順は共通にできる範囲で揃えますが、表示条件に該当しない行を空欄で埋めて増やすことはしません。

補正は、表示を絞る前にfamily全体へ掛けます。表示範囲を変更しても計算範囲・補正値・保存する比較は変わりません。保存結果には、表示条件から外れた計算済みの比較も残します。全体検定の表は、この下位検定用の表示フィルタでは絞りません。

現在の既定表示は「有意な比較のみ」です。表示範囲を導入する前の版の「補正後 `p < followup_alpha` を表示」とは異なります。有意な効果内の傾向の比較も見たい場合は「有意な効果の全比較」、傾向にとどまる効果の計算結果も見たい場合は「計算済みの全結果」を選びます。

### 信頼水準と尺度

信頼区間を計算する方法では、信頼水準を `1 - alpha` に揃えます。`trend_alpha` を `.10` にしても、自動で90%区間にはしません。ART-Cの補正後CIも保存・表示します。一般ランクとQuade法には、この実装では元尺度の平均差やCIを追加しません。

p値とCIの補正名は別に記録・表示します。冒頭に指定／固有の方式、実際のp値の補正、CIの補正をまとめます。後の効果で実際の方式が異なる場合だけ、その比較の注記に示します。1比較だけで調整不要な場合は実際の方式を `none` として記録することがあります。例えば、emmeansでHolmを指定した場合、p値はHolmでもCIはBonferroniとなります。補正pとCIが常に同じ補正方式という扱いはしません（emmeans公式 `summary.emmGrid` の仕様）。

ANOVA/ANCOVA・Welch・tの推定差は元尺度です。ART-Cの推定差は整列順位尺度です。Brunner–Munzelの推定値は `P(X < Y) + 0.5 P(X = Y)`、ATSの相対効果は等重みの参照分布に対する効果であり、MCTPではその差を比較します。これらを平均差や互いに同一の効果量として扱いません。

## 5. データと共変量

`preprocess()` で作成する `experiment` を受け取ります。最低限 `data`、`id`、`iv`、`dv`、`within`、`between`、`design` を持つ長い形式のリストを使い、共変量がある場合は `covariate` に列名を登録します。`iv_` / `dv_` の接頭辞による自動推定は、このtest側では行いません。

分析間の対象を揃えるには、前処理で参加者×条件の分析単位と欠損の扱いを揃えてください。

| 項目 | 現在の扱い |
|---|---|
| 同じID×条件に複数行 | 通常ANOVA/ANCOVAはafexが集約。一般ランク・t/Welch・Brunner–Munzel/ATS・Quadeは参加者×水準で平均集約 |
| ART/permucoの同じID×条件の複数行 | 分析前に従属変数をID×全要因の組合せごとに平均。反復回数が不揃いな場合も平均する |
| 通常ANOVA/ANCOVAの欠損 | DVごとに必要列の欠損を処理し、被験者内計画は完全な参加者を使用 |
| 一般ランク・t/Welch等の欠損 | 方法ごとに欠損を除外。対応ありは全水準の揃ったIDを全体・下位検定で使用 |
| ARTの欠損 | ARToolのエラーとして停止する場合がある。前処理で整える |
| permucoの欠損／欠けた被験者内条件 | 明示的に停止。自動の補完・不完全デザインへの切替はしない |
| IDと要因 | 対応なしで同じIDが複数群に入る構造は不可。被験者間要因はID内で一定 |
| DV／共変量 | 数値。無限値は不可。共変量は各参加者内で一定の値を想定 |

通常ANOVAが保存する `analysis_data` は、afexに渡す段階のデータです。afex内部で平均集約したデータと常に同一ではありません。比較方法間で、保存行数だけから標本数が一致すると判断しないでください。

ART・permucoは同じID×条件の観測が複数ある場合、従属変数ごとに平均してからモデルに渡します。集約した場合は元の行数、集約後のセル数、最小・最大反復回数を結果の `cell_replication` に保存します。対象の従属変数に欠損・無限値がある場合は停止します。置換ANCOVAの共変量は欠損がなく参加者内で一定である必要があります。条件が欠けた反復測定データは各方法の入力条件に従います。試行数が少ないセルの平均には不確かさがあるため、反復回数も確認してください。

### 共変量を使う3方法

- 通常ANCOVAは共変量を平均中心化します。既定では被験者間要因×共変量を含めません。`covariate_interactions = TRUE` で追加できます。被験者内要因×共変量はafexの反復測定モデルに従い常に含みます。下位検定は中心化後0、すなわち共変量の平均値で評価します。
- 置換ANCOVAも平均中心化します。既定では実験要因×各共変量の交互作用を含み、`covariate_interactions = FALSE` で除けます。共変量同士の交互作用は追加しません。対応なしでは `freedman_lane`、被験者内を含む場合は `Rd_kheradPajouh_renaud` を既定にしています。
- 順位ANCOVAは `npANCOVA::Quade()` を使います。対応なし1要因のみです。共変量との交互作用を要求した場合は停止します。共変量の群間分布に関する仮定を自動で検証・保証する機能はありません。

共変量列が登録されていても、これら以外の方法は調整に使用しません。Quade法の `covariate_table` は、群を含まない「順位DVを順位共変量で説明する補助回帰」の結果です。通常ANCOVAで群を調整した共変量の部分効果検定とは区別してください。Quadeの群効果と補助回帰は別に保存します。

## 6. 実行例・再読み込み・保存結果

```r
# 一括loaderで読み込み・再読み込み
source("script/load_script.R", encoding = "UTF-8")

# test部分だけ再読み込みする場合（helperと子スクリプトは同じフォルダ）
source("script/test/test_helper.R", encoding = "UTF-8")

# 対応なし1要因の模擬データからexperimentを作成
# このデータは下記の順位ANCOVA・一般ランクにも使用可能
data <- read.csv("test_data/between_3.csv")
experiment <- preprocess(data, id = "ID", iv = "group", dv = "score",
  covariate = "covariate", select = "arguments", show = FALSE)

# 方法を省略するとGUI
result <- analyze(experiment)

# 複数方法。同じ球面性の適用ルールを使用
result <- analyze(experiment, method = c("anova", "art", "permuco"),
                  correction = "auto", adjust = "holm", np = 5000, seed = 42)

# 共変量分析は3種類のうち1つを選択
result <- analyze(experiment, method = "ancova", correction = "auto")
result <- analyze(experiment, method = "permuco_ancova", np = 5000, seed = 42)
result <- analyze(experiment, method = "rank_ancova", adjust = "holm")

# GUIなし・表示なし
result <- analyze(experiment, method = "rank", select = "arguments", show = FALSE)

# 計算と表示を別々に指定
result <- analyze(experiment, method = "anova",
                  alpha = .05, trend_alpha = .10,
                  posthoc = TRUE, followup_scope = "followup", display_scope = "effect")

# 全効果の下位検定を計算し、計算済みの比較をすべて表示
result <- analyze(experiment, method = "art", followup_scope = "all", display_scope = "all")

# 全体検定だけ行う
result <- analyze(experiment, method = "anova", posthoc = FALSE)
dir.create("results", showWarnings = FALSE)
saveRDS(result, "results/analysis_result.rds")
```

パッケージが不足している場合は、該当方法の実行時に停止し、インストール方法を示します。例：`install.packages("permuco")`。GUIにはtcltkを利用できるデスクトップR/RStudioが必要です。

| 方法 | 主な外部パッケージ |
|---|---|
| t、2群Welch | R標準のstats |
| ANOVA/ANCOVA | afex、emmeans、car |
| ART/ART-C | ARTool、emmeans、car |
| permuco | permuco |
| 一般ランクのDunn下位検定 | PMCMRplus |
| Welchの3群以上 | stats、下位検定にPMCMRplus |
| Brunner–Munzel | lawstat |
| ATS＋MCTP | rankFD |
| 順位ANCOVA | npANCOVA |

### 段階表示

DVごとに、該当する段階だけを表示します。球面性→診断／t・分散分析→共変量込みの分析／順位検定→下位検定の順です。比較できる複数方法は同じ見出しにまとめ、表は方法ごとの縦ブロックにします。全結果の計算終了後にまとめて出す方式ではなく、段階ごとに表示します。ただし球面性評価にもモデル推定が必要です。

下位検定は効果単位で計算・表示します。補正・信頼水準・有意な効果の概要は冒頭に集約し、各比較群で有意／有意傾向のペアを示します。p値と偏η²は先頭の0を省略し、p値の左に効果量、右端に有意記号を置きます。Mauchlyのpの列見出しは球面性の表内では単に `p` とします。

分散分析・共分散分析の見出しでは、`auto`という設定名だけを出すのではなく、例えば「球面性補正：自動選択 → GG / HF / 未補正（効果別は補正列）」と、実際に使用した方式を列挙します。各効果の「補正」列にもGG／HF／未補正／算出不可を示します。被験者内要因がない場合は「球面性補正：対象なし」と表示します。表示だけの変更で、補正値やp値は変えません。

有意記号は慣例的な固定境界（`*** < .001`、`** < .01`、`* < .05`、`. < .10`）です。`alpha` を変更した場合の有意／傾向の一覧は、指定した基準に従います。p値は既定3桁、`.001` 未満は `<.001` で表示し、保存値は丸めません。`p_digits` は1以上の整数です。

### 戻り値

| 指定 | 参照例 |
|---|---|
| 単一方法・DV=`Y` | `result$Y` |
| 複数方法・DV=`Y` | `result$anova$Y`、`result$art$Y` |
| GUIから実行 | `attr(result, "selection")$code` に再実行コード |

各結果は数値表、モデル／元の検定オブジェクト、使用データ、警告・通知を保持します。内部の遅延計算関数 `compute_followups` は、返却時には除去します。

| 主な項目 | 内容 |
|---|---|
| `test_table` / `anova_table` | 全体検定。補正後pと、必要に応じ未補正p・自由度 |
| `anova_table_uncorrected` | ANOVA/ART等の未補正の全体表 |
| `sphericity_table` | W、Mauchlyのp、GG/HFなどの推定値と被験者内の対比自由度 `Within.df` |
| `sphericity_decisions` | 効果ごとの要求設定、Mauchlyのp、GG/HFの推定ε、適用補正、使用ε、選択規則、`.75`の閾値、判定理由 |
| `sphericity` | ANOVA/ANCOVA/ARTに要求した適用方針 |
| `contrasts[[効果名]]$unadjusted` | 未補正の比較表（固有補正法を除く） |
| `contrasts[[効果名]]$formatted` | 補正後の全比較。CI、family等も方法に応じ保存 |
| `confidence_level` / `p_adjustment` / `ci_adjustment` | 各比較の信頼水準・実際の補正。順位検定／QuadeでCIを算出しない場合は `ci_adjustment = "not_computed"` |
| `posthoc` / `followup_scope` / `display_scope` | 下位検定の実行可否、計算範囲、表示範囲 |
| `alpha` / `trend_alpha` | 有意水準、有意傾向と下位検定対象選択の共通上限 |
| `followup_alpha` | 旧結果との互換名。`trend_alpha` と同じ共通基準を保存 |
| `follow_up_effects` | その方法で選択した下位検定の計算対象効果。オフなら空 |
| `requested_adjust` / `adjustment_engine` | 要求した共通補正と計算元。固定補正に要求した共通設定を重ねたという意味にはしない |
| `family` / `family_size` | 比較ごとの補正範囲。比較表と比較オブジェクトに保持 |
| `covariate` / `covariate_adjusted` | 共変量の指定と、実際に調整したか |
| `covariate_table` / `covariate_regression` | Quade法の補助回帰。通常ANCOVAの共変量行とは意味が異なる |
| `warnings` / `messages` / `error_details` | 計算上の警告、説明、対応方法が保存する詳細エラー |

通常ANCOVA・置換ANCOVAの共変量と交互作用の検定行は、全体の `test_table` / `anova_table` に保存します。表示されない・有意でない共変量も削除しません。置換法の表では、残差行を保存したまま、表示時だけ省略して分母自由度を補います。

球面性の判定理由には `sphericity_rejected`、`sphericity_not_rejected`、`mauchly_unavailable_epsilon_selection`、`epsilon_unavailable`、`explicit_correction`、`not_requested`、`not_applicable` を使います。`selection_rule` にGGのεが`.75`以下／超であったことも保存します。検定不能な値はNAとして保持し、別検定へ黙って切り替えません。ただし全体モデルのエラーは方法によって呼び出し自体を停止するため、必ず部分結果が返るという仕様ではありません。

## 7. 論文・公式資料と実装の対応

以下は方法の原著・一次資料です。一般的なt・ANOVA・ランク・ANCOVA以外の手法について、方法の根拠とこの実装の範囲を分けて記します。論文の存在を、ラッパー独自の全処理の検証済みという意味には使いません。

### Welch、Games–Howell

1. Welch, B. L. (1947). *The Generalization of ‘Student’s’ Problem When Several Different Population Variances Are Involved.* Biometrika, 34(1–2), 28–35. [DOI: 10.1093/biomet/34.1-2.28](https://doi.org/10.1093/biomet/34.1-2.28)。不等分散下の平均比較の基礎。実装では2群に `stats::t.test(var.equal = FALSE)` を使います。
2. Welch, B. L. (1951). *On the Comparison of Several Mean Values: An Alternative Approach.* Biometrika, 38(3–4), 330–336. [DOI: 10.1093/biomet/38.3-4.330](https://doi.org/10.1093/biomet/38.3-4.330)。多群の不等分散平均比較。3群以上に `stats::oneway.test(var.equal = FALSE)` を使います。
3. Games, P. A., & Howell, J. F. (1976). *Pairwise Multiple Comparison Procedures with Unequal N’s and/or Variances: A Monte Carlo Study.* Journal of Educational Statistics, 1(2), 113–125. [DOI: 10.3102/10769986001002113](https://doi.org/10.3102/10769986001002113)。不等分散・不均衡標本数での多重比較。p値は `PMCMRplus::gamesHowellTest()`、区間は各ペアの標準誤差・Welch自由度とスチューデント化範囲分布で算出します。有限標本で厳密な被覆率を保証する区間という扱いはしません。

### ART、ART-C

4. Wobbrock, J. O., Findlater, L., Gergle, D., & Higgins, J. J. (2011). *The Aligned Rank Transform for Nonparametric Factorial Analyses Using Only ANOVA Procedures.* CHI 2011, 143–146. [DOI: 10.1145/1978942.1978963](https://doi.org/10.1145/1978942.1978963)、[著者公開PDF](https://faculty.washington.edu/wobbrock/pubs/chi-11.06.pdf)。要因ごとに整列して順位変換するARTの根拠です。
5. Elkin, L. A., Kay, M., Higgins, J. J., & Wobbrock, J. O. (2021). *An Aligned Rank Transform Procedure for Multifactor Contrast Tests.* UIST 2021, 754–768. [DOI: 10.1145/3472749.3474784](https://doi.org/10.1145/3472749.3474784)、[著者公開PDF](https://faculty.washington.edu/wobbrock/pubs/uist-21.pdf)。多要因の比較にはART-Cを使う根拠です。

[ARTool公式資料](https://depts.washington.edu/acelab/proj/art/index.html)に沿って `art()`、`art.con()`、`artlm.con()` を使用します。ただし、被験者内を含む多要因ART-Cでは、このラッパーがART-C順位データを元の要因表現と誤差層で `aov()` に再推定して比較します。この再推定と前述のGG/HFは独自の処理です。単にARToolの出力を無加工で表示しているという説明にはしません。

### permucoによる置換ANOVA／ANCOVA

6. Frossard, J., & Renaud, O. (2021). *Permutation Tests for Regression, ANOVA, and Comparison of Signals: The permuco Package.* Journal of Statistical Software, 99(15), 1–32. [DOI: 10.18637/jss.v099.i15](https://doi.org/10.18637/jss.v099.i15)。モデル・共変量・反復測定に応じた置換法とR実装の一次資料です。

計算は `permuco::aovperm()` に委譲します。実装するのは全体の要因・共変量の検定で、独自の置換下位検定は追加していません。[公式aovpermリファレンス](https://search.r-project.org/CRAN/refmans/permuco/html/aovperm.html)にある各置換法の条件に従います。

### Brunner–Munzel

7. Brunner, E., & Munzel, U. (2000). *The Nonparametric Behrens–Fisher Problem: Asymptotic Theory and a Small-Sample Approximation.* Biometrical Journal, 42(1), 17–25. [原著DOI](https://doi.org/10.1002/%28SICI%291521-4036%28200001%2942%3A1%3C17%3A%3AAID-BIMJ17%3E3.0.CO%3B2-U)。異なる分布・分散を許した2群の相対効果検定とt近似の根拠です。
8. Neubert, K., & Brunner, E. (2007). *A studentized permutation test for the non-parametric Behrens–Fisher problem.* Computational Statistics & Data Analysis, 51(10), 5192–5204. [DOI: 10.1016/j.csda.2006.05.024](https://doi.org/10.1016/j.csda.2006.05.024)。studentized置換版の関連原著で、lawstatの参考文献にも挙げられています。**この実装は置換版ではなく、lawstatによるt近似です。**

[lawstat公式リファレンス](https://search.r-project.org/CRAN/refmans/lawstat/html/brunner.munzel.test.html)を使用します。小標本や分散推定不能の状況で、近似を無条件に信頼したり、自動で別の検定に切り替えたりする実装ではありません。

### ATS＋MCTP

9. Brunner, E., Konietschke, F., Pauly, M., & Puri, M. L. (2017). *Rank-Based Procedures in Factorial Designs: Hypotheses About Non-Parametric Treatment Effects.* Journal of the Royal Statistical Society: Series B, 79(5), 1463–1485. [DOI: 10.1111/rssb.12222](https://doi.org/10.1111/rssb.12222)。非パラメトリック相対効果、pseudo-rankに基づく要因検定の理論です。
10. Konietschke, F., Hothorn, L. A., & Brunner, E. (2012). *Rank-based multiple test procedures and simultaneous confidence intervals.* Electronic Journal of Statistics, 6, 738–759. [DOI: 10.1214/12-EJS691](https://doi.org/10.1214/12-EJS691)。順位に基づく多重検定と同時信頼区間の根拠です。
11. Konietschke, F., Pauly, M., Bathke, A. C., Friedrich, S., & Brunner, E. (2023). *rankFD: An R Software Package for Nonparametric Analysis of General Factorial Designs.* The R Journal. [DOI: 10.32614/RJ-2023-029](https://doi.org/10.32614/RJ-2023-029)、[原文PDF](https://journal.r-project.org/articles/RJ-2023-029/RJ-2023-029.pdf)。ATSとMCTPを提供する公式実装と使用例です。

実装設定は `effect = "unweighted"`、`hypothesis = "H0p"`、ATSのF近似、全ペアのMCTPで `sci.method = "fisher"` です。rankFD自体は一般の要因計画を扱えますが、このスクリプトの入口は対応なし1要因・3群以上に限定します。MCTPの多変量tの数値計算に乱数が関係するため、再現時は `seed` を指定してください。「すべてのデータにATS＋MCTPを使うことを著者らが一律に推奨した」という意味ではなく、理論と公式実装に沿った構成です。

### Quadeの順位ANCOVA

12. Quade, D. (1967). *Rank Analysis of Covariance.* Journal of the American Statistical Association, 62(320), 1187–1200. [DOI: 10.1080/01621459.1967.10500925](https://doi.org/10.1080/01621459.1967.10500925)。順位ANCOVAの原著です。

[npANCOVA公式リファレンス](https://cran.r-project.org/web/packages/npANCOVA/refman/npANCOVA.html)の `Quade()` を使います。下位検定は同パッケージの `pairwise_npANCOVA(method = "Quade")` と同様に2群ずつ再順位化・再推定し、helperで補正します。これは反復測定の `stats::quade.test()` とは別の方法です。全体と下位検定で同じ順位尺度の差を推定しているとは解釈しません。

### 球面性・共通補正の資料

13. Greenhouse, S. W., & Geisser, S. (1959). *On Methods in the Analysis of Profile Data.* Psychometrika, 24, 95–112. [DOI: 10.1007/BF02289823](https://doi.org/10.1007/BF02289823)。GG補正の原著です。
14. Huynh, H., & Feldt, L. S. (1976). *Estimation of the Box Correction for Degrees of Freedom from Sample Data in Randomized Block and Split-Plot Designs.* Journal of Educational Statistics, 1(1), 69–82. [DOI: 10.3102/10769986001001069](https://doi.org/10.3102/10769986001001069)。HF補正の原著です。

多重比較は[R公式p.adjust](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html)、CIとp値の調整の違いは[emmeans公式summary.emmGrid](https://search.r-project.org/CRAN/refmans/emmeans/html/summary.emmGrid.html)を参照してください。基本手法の実装はR公式の `t.test`、`oneway.test`、`wilcox.test`、`kruskal.test`、`friedman.test` とafexの仕様に従います。

## 8. 検証状況と制約

現在のリポジトリには自動テスト一式は含まれていません。

下位検定の計算・表示範囲、補正、欠損処理は上記の実装仕様です。方法間で欠損処理や多重比較familyが異なる点、ARTへのGG/HF適用と被験者内ART-Cの再推定に独自処理を含む点に留意してください。独自拡張を含む全条件での検定サイズ・検出力を検証したという意味ではありません。
