# surviq-stats

Survey and thesis statistics in TypeScript, with no dependencies. These are the same formulas behind the free calculators at [SurvIQ](https://surviqapp.com/en/tools), checked against SPSS and SciPy.

```bash
npm install surviq-stats
```

```ts
import { cronbachAlpha, aikenV, kendallW, mannWhitney, sampleSize } from 'surviq-stats'

sampleSize(null, '95', 5)                 // 385 responses for ±5% at 95%
cronbachAlpha([[1, 2, 3], [2, 3, 4], [3, 3, 5], [4, 5, 5]]).alpha // 0.949
aikenV([4, 5, 4, 5, 5], 1, 5)             // { v: 0.9, lower: 0.699, upper: 0.972, judges: 5 }
kendallW([[5, 5, 4], [4, 4, 4], [3, 2, 3], [2, 3, 2], [1, 1, 2]]).w // 0.920
mannWhitney([12, 15, 14, 10, 13], [18, 17, 20, 16, 19]).exact       // 0.0079
```

Tables are arrays of rows. Unless a function says otherwise, rows are people (respondents) and columns are items. For expert judgment (`aikenV`, `contentValidity`, `kendallW`, `binomialJudges`), rows are items and columns are judges.

## What's inside

### Reliability and validity

| Function | What it does | Online calculator |
| --- | --- | --- |
| `cronbachAlpha(rows)` | Cronbach's alpha, with alpha if item deleted | [Cronbach's alpha](https://surviqapp.com/en/tools/cronbach-alpha-calculator) |
| `mcdonaldOmega(rows)` | McDonald's omega (one factor, minres) | [McDonald's omega](https://surviqapp.com/en/tools/mcdonald-omega-calculator) |
| `kr20(rows)` | KR-20 and KR-21 for right/wrong items, with difficulty and discrimination | [KR-20](https://surviqapp.com/en/tools/kr20-calculator) |
| `aikenV(ratings, low, high)` | Aiken's V with Penfield and Giacobbi score interval | [Aiken's V](https://surviqapp.com/en/tools/aiken-v-calculator) |
| `contentValidity(rows, vmax)` | Hernández-Nieto content validity coefficient (CVC) | [CVC](https://surviqapp.com/en/tools/content-validity-coefficient-calculator) |
| `kendallW(rows)` | Kendall's W with tie correction, chi-square and p-value | [Kendall's W](https://surviqapp.com/en/tools/kendalls-w-calculator) |
| `binomialJudges(rows)` | Exact binomial test for yes/no expert judgment | [Binomial test](https://surviqapp.com/en/tools/binomial-test-expert-judgment) |
| `cohenKappa(a, b)`, `fleissKappa(rows)` | Inter-rater agreement | [Kappa](https://surviqapp.com/en/tools/cohen-fleiss-kappa-calculator) |

### Group comparisons and normality

| Function | What it does | Online calculator |
| --- | --- | --- |
| `shapiroWilk(values)`, `lilliefors(values)` | Normality tests (Shapiro-Wilk up to 50, K-S Lilliefors above) | [Normality test](https://surviqapp.com/en/tools/shapiro-wilk-normality-test) |
| `independentTTest(a, b)`, `pairedTTest(a, b)`, `leveneTest(a, b)` | Student and Welch t-tests, Levene, Cohen's d | [t-test](https://surviqapp.com/en/tools/t-test-calculator) |
| `mannWhitney(a, b)` | Mann-Whitney U with exact p for small samples | [Mann-Whitney](https://surviqapp.com/en/tools/mann-whitney-u-test-calculator) |
| `wilcoxonSignedRank(before, after)` | Wilcoxon signed-rank, exact p without ties | [Wilcoxon](https://surviqapp.com/en/tools/wilcoxon-signed-rank-test-calculator) |
| `oneWayAnova(groups)` | One-way ANOVA, Welch, Bonferroni pairs, eta squared | [ANOVA](https://surviqapp.com/en/tools/one-way-anova-calculator) |
| `kruskalWallis(groups)` | Kruskal-Wallis with Dunn-Bonferroni pairs | [Kruskal-Wallis](https://surviqapp.com/en/tools/kruskal-wallis-test-calculator) |
| `chiSquareTest(table)` | Chi-square of independence with Cramér's V | [Chi-square](https://surviqapp.com/en/tools/chi-square-calculator) |
| `correlation(x, y, 'pearson' \| 'spearman')` | Pearson or Spearman correlation with p-value | [Correlation](https://surviqapp.com/en/tools/correlation-calculator) |
| `twoProportionTest(x1, n1, x2, n2, '95')` | Compare two percentages | [Compare percentages](https://surviqapp.com/en/tools/compare-two-percentages) |

### Survey metrics and sample size

| Function | What it does | Online calculator |
| --- | --- | --- |
| `sampleSize(population, confidence, marginPct)` | Sample size with finite population correction | [Sample size](https://surviqapp.com/en/tools/sample-size-calculator) |
| `marginOfError(n, population, confidence)` | Margin of error for a given sample | [Margin of error](https://surviqapp.com/en/tools/margin-of-error-calculator) |
| `npsScore(promoters, passives, detractors)` | Net Promoter Score | [NPS](https://surviqapp.com/en/tools/nps-calculator) |
| `npsCompare(a, b)` | Is the difference between two NPS significant? | [Compare two NPS](https://surviqapp.com/en/tools/nps-comparison-calculator) |
| `csatScore(counts)`, `cesScore(counts)` | CSAT and CES | [CSAT](https://surviqapp.com/en/tools/csat-calculator) |
| `likertLevels(items, low, high, levels)` | Score ranges (norms) for Likert scales | [Likert ranges](https://surviqapp.com/en/tools/likert-scale-score-ranges) |
| `vanWestendorp(rows)` | Van Westendorp price sensitivity | [Van Westendorp](https://surviqapp.com/en/tools/van-westendorp-calculator) |

Distributions are also exported: `normalCdf`, `normalQuantile`, `tTestP`, `fTestP`, `chiSquareP`, `tCritical`.

## En español

`surviq-stats` reúne las fórmulas de las [calculadoras gratis para tesis de SurvIQ](https://surviqapp.com/es/herramientas): alfa de Cronbach, omega de McDonald, KR-20, V de Aiken, CVC de Hernández-Nieto, W de Kendall, prueba binomial para juicio de expertos, kappa, Shapiro-Wilk, t de Student, U de Mann-Whitney, Wilcoxon, ANOVA, Kruskal-Wallis, chi-cuadrado, correlación, tamaño de muestra y NPS. Si no programas, usa las calculadoras: hacen lo mismo en el navegador y te dan la tabla APA en Word.

## Accuracy

Every test in `test/` compares results with SciPy or with published worked examples. If you find a difference with SPSS, R or SciPy, please [open an issue](https://github.com/enzocaminomontes-netizen/surviq-stats/issues) with the data.

## License

MIT. Made by [SurvIQ](https://surviqapp.com), an online survey tool with built-in analysis.
