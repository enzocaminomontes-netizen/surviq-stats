import { chiSquareTest, correlation, normalCdf } from './stats.js'

// Cálculos del piloto y de la aplicación del Modo Tesis: confiabilidad (alfa de
// Cronbach y omega de McDonald), puntajes por dimensión y variable, niveles
// (baremos), normalidad, correlaciones y chi cuadrado. Todo es puro (sin base de
// datos) para poder probarlo y usarlo igual en pantallas, exportaciones y la IA.

// ── Básicos ─────────────────────────────────────────────────

export const mean = (x: number[]) => x.reduce((a, b) => a + b, 0) / x.length

export function sd(x: number[]): number {
  if (x.length < 2) return 0
  const m = mean(x)
  return Math.sqrt(x.reduce((s, v) => s + (v - m) ** 2, 0) / (x.length - 1))
}

export function median(x: number[]): number {
  const s = [...x].sort((a, b) => a - b)
  const mid = Math.floor(s.length / 2)
  return s.length % 2 ? s[mid] : (s[mid - 1] + s[mid]) / 2
}

function pearson(x: number[], y: number[]): number | null {
  const mx = mean(x)
  const my = mean(y)
  let sxy = 0, sxx = 0, syy = 0
  for (let i = 0; i < x.length; i++) {
    sxy += (x[i] - mx) * (y[i] - my)
    sxx += (x[i] - mx) ** 2
    syy += (y[i] - my) ** 2
  }
  if (sxx === 0 || syy === 0) return null
  return sxy / Math.sqrt(sxx * syy)
}

// Inversa de la normal estándar (Acklam). Error relativo < 1.2e-9.
export function normalQuantile(p: number): number {
  const a = [-39.69683028665376, 220.9460984245205, -275.9285104469687, 138.357751867269, -30.66479806614716, 2.506628277459239]
  const b = [-54.47609879822406, 161.5858368580409, -155.6989798598866, 66.80131188771972, -13.28068155288572]
  const c = [-0.007784894002430293, -0.3223964580411365, -2.400758277161838, -2.549732539343734, 4.374664141464968, 2.938163982698783]
  const d = [0.007784695709041462, 0.3224671290700398, 2.445134137142996, 3.754408661907416]
  const low = 0.02425
  if (p <= 0) return -Infinity
  if (p >= 1) return Infinity
  if (p < low) {
    const q = Math.sqrt(-2 * Math.log(p))
    return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
  }
  if (p > 1 - low) {
    const q = Math.sqrt(-2 * Math.log(1 - p))
    return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
  }
  const q = p - 0.5
  const r = q * q
  return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1)
}

// La normal de stats.ts tiene un error de ~1e-7, suficiente para valores p.
const pnorm = normalCdf

// ── Confiabilidad ───────────────────────────────────────────

function alphaOf(cols: number[][]): number | null {
  const k = cols.length
  if (k < 2) return null
  const n = cols[0].length
  const itemVar = cols.reduce((s, c) => s + sd(c) ** 2, 0)
  const totals = Array.from({ length: n }, (_, i) => cols.reduce((s, c) => s + c[i], 0))
  const totalVar = sd(totals) ** 2
  if (totalVar === 0) return null
  return (k / (k - 1)) * (1 - itemVar / totalVar)
}

export type ItemReliability = {
  mean: number
  sd: number
  // Correlación ítem-total corregida (el total sin el ítem).
  citc: number | null
  alphaIfDeleted: number | null
}

export type Reliability = {
  alpha: number | null
  omega: number | null
  k: number
  n: number
  items: ItemReliability[]
}

// rows: una fila por persona, una columna por ítem (ya recodificados los inversos).
export function reliability(rows: number[][], withOmega = true): Reliability {
  const k = rows[0]?.length ?? 0
  const n = rows.length
  const cols = Array.from({ length: k }, (_, j) => rows.map(r => r[j]))
  const enough = n >= 2 && k >= 2
  const alpha = enough ? alphaOf(cols) : null
  const items = cols.map((c, j) => {
    const rest = rows.map(r => r.reduce((s, v, i) => (i === j ? s : s + v), 0))
    return {
      mean: n ? mean(c) : 0,
      sd: sd(c),
      citc: enough ? pearson(c, rest) : null,
      alphaIfDeleted: enough && k > 2 ? alphaOf(cols.filter((_, i) => i !== j)) : null,
    }
  })
  return { alpha, omega: withOmega && n >= 3 && k >= 3 ? omegaTotal(cols) : null, k, n, items }
}

// Omega de McDonald con un solo factor, por ejes principales (como hace SPSS con
// "ejes principales" y un factor; jamovi usa mínimos residuos y puede diferir en el
// tercer decimal).
export function omegaTotal(cols: number[][]): number | null {
  const k = cols.length
  const R: number[][] = Array.from({ length: k }, () => new Array<number>(k).fill(0))
  for (let i = 0; i < k; i++) {
    R[i][i] = 1
    for (let j = i + 1; j < k; j++) {
      const r = pearson(cols[i], cols[j])
      if (r === null) return null
      R[i][j] = R[j][i] = r
    }
  }
  // Comunalidad inicial: la mayor correlación de cada ítem.
  let h2 = R.map((row, i) => Math.max(...row.filter((_, j) => j !== i).map(Math.abs)))
  let loadings = new Array<number>(k).fill(0)
  for (let iter = 0; iter < 500; iter++) {
    const M = R.map((row, i) => row.map((v, j) => (i === j ? h2[i] : v)))
    // Primer autovector por el método de la potencia.
    let v = new Array<number>(k).fill(1 / Math.sqrt(k))
    let lambda = 0
    for (let p = 0; p < 300; p++) {
      const w = M.map(row => row.reduce((s, x, j) => s + x * v[j], 0))
      const norm = Math.sqrt(w.reduce((s, x) => s + x * x, 0))
      if (norm === 0) return null
      const next = w.map(x => x / norm)
      const delta = next.reduce((s, x, j) => s + Math.abs(x - v[j]), 0)
      v = next
      lambda = norm
      if (delta < 1e-10) break
    }
    if (!(lambda > 0)) return null
    loadings = v.map(x => x * Math.sqrt(lambda))
    const next = loadings.map(l => Math.min(l * l, 0.999))
    const change = next.reduce((s, x, i) => s + Math.abs(x - h2[i]), 0)
    h2 = next
    if (change < 1e-8) break
  }
  const sum = loadings.reduce((a, b) => a + b, 0)
  const signed = sum < 0 ? loadings.map(l => -l) : loadings
  const s = signed.reduce((a, b) => a + b, 0)
  const unique = signed.reduce((acc, l) => acc + (1 - l * l), 0)
  const omega = (s * s) / (s * s + unique)
  return Number.isFinite(omega) ? omega : null
}

// ── Normalidad ──────────────────────────────────────────────

export type Normality = { test: 'shapiro' | 'ks'; statistic: number; pValue: number; n: number }

// Shapiro-Wilk con la aproximación de Royston (algoritmo AS R94), de 3 a 5000 casos.
export function shapiroWilk(values: number[]): { w: number; pValue: number } | null {
  const n = values.length
  if (n < 3 || n > 5000) return null
  const x = [...values].sort((a, b) => a - b)
  const range = x[n - 1] - x[0]
  if (range === 0) return null
  const m = Array.from({ length: n }, (_, i) => normalQuantile((i + 1 - 0.375) / (n + 0.25)))
  const mm = m.reduce((s, v) => s + v * v, 0)
  const a = new Array<number>(n).fill(0)
  if (n === 3) {
    a[0] = -Math.SQRT1_2
    a[2] = Math.SQRT1_2
  } else {
    const u = 1 / Math.sqrt(n)
    const c = m.map(v => v / Math.sqrt(mm))
    const an = c[n - 1] + 0.221157 * u - 0.147981 * u ** 2 - 2.07119 * u ** 3 + 4.434685 * u ** 4 - 2.706056 * u ** 5
    if (n > 5) {
      const an1 = c[n - 2] + 0.042981 * u - 0.293762 * u ** 2 - 1.752461 * u ** 3 + 5.682633 * u ** 4 - 3.582633 * u ** 5
      const phi = (mm - 2 * m[n - 1] ** 2 - 2 * m[n - 2] ** 2) / (1 - 2 * an ** 2 - 2 * an1 ** 2)
      for (let i = 2; i < n - 2; i++) a[i] = m[i] / Math.sqrt(phi)
      a[0] = -an; a[1] = -an1; a[n - 2] = an1; a[n - 1] = an
    } else {
      const phi = (mm - 2 * m[n - 1] ** 2) / (1 - 2 * an ** 2)
      for (let i = 1; i < n - 1; i++) a[i] = m[i] / Math.sqrt(phi)
      a[0] = -an; a[n - 1] = an
    }
  }
  const mx = mean(x)
  const ssq = x.reduce((s, v) => s + (v - mx) ** 2, 0)
  const w = Math.min(1, (a.reduce((s, ai, i) => s + ai * x[i], 0) ** 2) / ssq)

  let pValue: number
  if (n === 3) {
    pValue = Math.max(0, (6 / Math.PI) * (Math.asin(Math.sqrt(w)) - Math.asin(Math.sqrt(0.75))))
  } else if (n <= 11) {
    const gamma = 0.459 * n - 2.273
    const mu = 0.544 - 0.39978 * n + 0.025054 * n ** 2 - 0.0006714 * n ** 3
    const sigma = Math.exp(1.3822 - 0.77857 * n + 0.062767 * n ** 2 - 0.0020322 * n ** 3)
    const inner = gamma - Math.log(1 - w)
    if (!(inner > 0)) pValue = 0
    else pValue = 1 - pnorm((-Math.log(inner) - mu) / sigma)
  } else {
    const ln = Math.log(n)
    const mu = 0.0038915 * ln ** 3 - 0.083751 * ln ** 2 - 0.31082 * ln - 1.5861
    const sigma = Math.exp(0.0030302 * ln ** 2 - 0.082676 * ln - 0.4803)
    pValue = w >= 1 ? 1 : 1 - pnorm((Math.log(1 - w) - mu) / sigma)
  }
  return { w, pValue: Math.min(1, Math.max(0, pValue)) }
}

// Kolmogorov-Smirnov con la corrección de Lilliefors (media y desviación de la
// muestra), con el valor p de Dallal y Wilkinson como en SPSS y R (nortest).
export function lilliefors(values: number[]): { d: number; pValue: number } | null {
  const n = values.length
  if (n < 5) return null
  const x = [...values].sort((a, b) => a - b)
  const m = mean(x)
  const s = sd(x)
  if (s === 0) return null
  let d = 0
  x.forEach((v, i) => {
    const p = pnorm((v - m) / s)
    d = Math.max(d, (i + 1) / n - p, p - i / n)
  })
  const [kd, nd] = n <= 100 ? [d, n] : [d * (n / 100) ** 0.49, 100]
  let pValue = Math.exp(-7.01256 * kd * kd * (nd + 2.78019) + 2.99587 * kd * Math.sqrt(nd + 2.78019) - 0.122119 + 0.974598 / Math.sqrt(nd) + 1.67997 / nd)
  if (pValue > 0.1) {
    const kk = (Math.sqrt(n) - 0.01 + 0.85 / Math.sqrt(n)) * d
    if (kk <= 0.302) pValue = 1
    else if (kk <= 0.5) pValue = 2.76773 - 19.828315 * kk + 80.709644 * kk ** 2 - 138.55152 * kk ** 3 + 81.218052 * kk ** 4
    else if (kk <= 0.9) pValue = -4.901232 + 40.662806 * kk - 97.490286 * kk ** 2 + 94.029866 * kk ** 3 - 32.355711 * kk ** 4
    else if (kk <= 1.31) pValue = 6.198765 - 19.558097 * kk + 23.186922 * kk ** 2 - 12.234627 * kk ** 3 + 2.423045 * kk ** 4
    else pValue = 0
  }
  return { d, pValue: Math.min(1, Math.max(0, pValue)) }
}

// Regla habitual en tesis: Shapiro-Wilk hasta 50 casos, Kolmogorov-Smirnov
// (Lilliefors) con más de 50.
export const SHAPIRO_MAX_N = 50

export function normality(values: number[]): Normality | null {
  if (values.length <= SHAPIRO_MAX_N) {
    const r = shapiroWilk(values)
    return r ? { test: 'shapiro', statistic: r.w, pValue: r.pValue, n: values.length } : null
  }
  const r = lilliefors(values)
  return r ? { test: 'ks', statistic: r.d, pValue: r.pValue, n: values.length } : null
}

// ── Puntajes, niveles y análisis del instrumento ────────────

export type ScoringMode = 'sum' | 'mean'
export type LevelMethod = 'equal' | 'percentile'

export type AnalysisItem = {
  id: string
  dimensionId: string | null
  reverse: boolean
  // Rango de respuesta del ítem (para invertir y para los baremos).
  min: number
  max: number
}

export type AnalysisScope = {
  id: string
  kind: 'total' | 'variable' | 'dimension'
  name: string
  // Variable a la que pertenece la dimensión.
  parentId: string | null
  itemIds: string[]
}

export type Respondent = {
  id: string
  // Respuesta cruda por ítem (sin invertir).
  values: Record<string, number | undefined>
  // Variables de grupo: id de la pregunta -> etiqueta elegida.
  groups: Record<string, string | undefined>
}

export type GroupVariable = { id: string; name: string; options: string[] }

export type Level = { name: string; from: number; to: number; count: number; pct: number }

export type ScopeResult = {
  scope: AnalysisScope
  // Personas con todos los ítems de esta escala respondidos.
  n: number
  k: number
  theoreticalMin: number
  theoreticalMax: number
  mean: number | null
  sd: number | null
  median: number | null
  observedMin: number | null
  observedMax: number | null
  reliability: Reliability
  levels: Level[]
  normality: Normality | null
  // Puntaje de cada persona (por id), para correlacionar y exportar.
  scores: Record<string, number>
  levelOf: Record<string, number>
}

export type CorrelationResult = {
  a: string
  b: string
  method: 'pearson' | 'spearman'
  r: number
  pValue: number
  n: number
}

export type CrosstabResult = {
  groupId: string
  groupName: string
  scopeId: string
  // Filas: opciones del grupo; columnas: niveles.
  table: number[][]
  rows: string[]
  chi2: number
  df: number
  pValue: number
  cramerV: number
  lowExpectedPct: number
  n: number
}

export type AnalysisOptions = {
  scoring: ScoringMode
  levelMethod: LevelMethod
  levelNames: string[]
  // Nivel de confianza de las pruebas (95 => alfa 0.05).
  confidence: '90' | '95' | '99'
  correlationMethod: 'auto' | 'pearson' | 'spearman'
  omega: boolean
}

export type Analysis = {
  respondents: number
  scopes: ScopeResult[]
  correlations: CorrelationResult[]
  crosstabs: CrosstabResult[]
  alphaLevel: number
}

export function recode(item: AnalysisItem, value: number): number {
  return item.reverse ? item.min + item.max - value : value
}

const round6 = (v: number) => Math.round(v * 1e6) / 1e6

// Cortes de los niveles. Rangos iguales: sobre el puntaje teórico (baremo clásico).
// Percentiles: sobre los puntajes de la muestra.
export function levelCuts(method: LevelMethod, levels: number, theoreticalMin: number, theoreticalMax: number, scores: number[]): { from: number; to: number }[] {
  if (method === 'percentile' && scores.length >= levels) {
    const s = [...scores].sort((a, b) => a - b)
    const q = (p: number) => {
      const pos = (s.length - 1) * p
      const lo = Math.floor(pos)
      return s[lo] + (s[Math.min(lo + 1, s.length - 1)] - s[lo]) * (pos - lo)
    }
    const cuts = Array.from({ length: levels - 1 }, (_, i) => round6(q((i + 1) / levels)))
    return Array.from({ length: levels }, (_, i) => ({
      from: i === 0 ? theoreticalMin : cuts[i - 1],
      to: i === levels - 1 ? theoreticalMax : cuts[i],
    }))
  }
  const width = (theoreticalMax - theoreticalMin) / levels
  return Array.from({ length: levels }, (_, i) => ({
    from: round6(theoreticalMin + i * width),
    to: round6(i === levels - 1 ? theoreticalMax : theoreticalMin + (i + 1) * width),
  }))
}

// Con puntajes enteros (sumas), el baremo se muestra en enteros sin solaparse:
// 10–30, 31–50, 51–70. El límite cuenta para el nivel de abajo, igual que levelIndex.
export function displayCuts(cuts: { from: number; to: number }[], integerScores: boolean): { from: number; to: number }[] {
  if (!integerScores) return cuts
  const out = cuts.map((c, i) => ({ from: i === 0 ? Math.ceil(c.from - 1e-9) : Math.floor(c.from + 1e-9) + 1, to: Math.floor(c.to + 1e-9) }))
  return out.every(c => c.from <= c.to) ? out : cuts
}

// Nivel de un puntaje: el primer corte cuyo "hasta" lo alcanza (los límites
// cuentan para el nivel de abajo, como en los baremos de las tesis).
export function levelIndex(cuts: { from: number; to: number }[], score: number): number {
  const i = cuts.findIndex(c => score <= c.to + 1e-9)
  return i === -1 ? cuts.length - 1 : i
}

export function analyze(input: {
  items: AnalysisItem[]
  scopes: AnalysisScope[]
  respondents: Respondent[]
  groups: GroupVariable[]
  options: AnalysisOptions
}): Analysis {
  const byId = new Map(input.items.map(i => [i.id, i]))
  const { options } = input
  const levelCount = Math.max(2, Math.min(5, options.levelNames.length))
  const alpha = 1 - Number(options.confidence) / 100

  const scopes: ScopeResult[] = input.scopes.map(scope => {
    const items = scope.itemIds.map(id => byId.get(id)).filter((i): i is AnalysisItem => Boolean(i))
    const k = items.length
    const complete = input.respondents.filter(r => items.every(it => typeof r.values[it.id] === 'number'))
    const rows = complete.map(r => items.map(it => recode(it, r.values[it.id] as number)))
    const sums = rows.map(row => row.reduce((a, b) => a + b, 0))
    const scoreValues = options.scoring === 'mean' ? sums.map(s => s / k) : sums
    const minSum = items.reduce((s, it) => s + it.min, 0)
    const maxSum = items.reduce((s, it) => s + it.max, 0)
    const theoreticalMin = options.scoring === 'mean' && k ? minSum / k : minSum
    const theoreticalMax = options.scoring === 'mean' && k ? maxSum / k : maxSum
    const cuts = levelCuts(options.levelMethod, levelCount, theoreticalMin, theoreticalMax, scoreValues)
    const integerScores = options.scoring === 'sum' && items.every(it => Number.isInteger(it.min) && Number.isInteger(it.max)) && scoreValues.every(v => Number.isInteger(v))
    const counts = new Array<number>(levelCount).fill(0)
    const scores: Record<string, number> = {}
    const levelOf: Record<string, number> = {}
    complete.forEach((r, i) => {
      scores[r.id] = scoreValues[i]
      const li = levelIndex(cuts, scoreValues[i])
      levelOf[r.id] = li
      counts[li]++
    })
    const n = complete.length
    return {
      scope,
      n,
      k,
      theoreticalMin,
      theoreticalMax,
      mean: n ? mean(scoreValues) : null,
      sd: n > 1 ? sd(scoreValues) : null,
      median: n ? median(scoreValues) : null,
      observedMin: n ? Math.min(...scoreValues) : null,
      observedMax: n ? Math.max(...scoreValues) : null,
      reliability: k ? reliability(rows, options.omega) : { alpha: null, omega: null, k: 0, n, items: [] },
      levels: displayCuts(cuts, integerScores).map((c, i) => ({ name: options.levelNames[i] ?? String(i + 1), from: c.from, to: c.to, count: counts[i], pct: n ? (counts[i] / n) * 100 : 0 })),
      normality: n >= 3 ? normality(scoreValues) : null,
      scores,
      levelOf,
    }
  })

  // Correlaciones: variable con variable, y cada dimensión con las otras variables.
  // Con una sola variable, dimensión con dimensión.
  const variables = scopes.filter(s => s.scope.kind === 'variable')
  const dims = scopes.filter(s => s.scope.kind === 'dimension')
  const pairs: [ScopeResult, ScopeResult][] = []
  for (let i = 0; i < variables.length; i++) for (let j = i + 1; j < variables.length; j++) pairs.push([variables[i], variables[j]])
  if (variables.length >= 2) {
    for (const d of dims) for (const v of variables) if (v.scope.id !== d.scope.parentId) pairs.push([d, v])
  } else {
    for (let i = 0; i < dims.length; i++) for (let j = i + 1; j < dims.length; j++) pairs.push([dims[i], dims[j]])
  }
  const isNormal = (s: ScopeResult) => s.normality !== null && s.normality.pValue > alpha
  const correlations: CorrelationResult[] = []
  for (const [a, b] of pairs) {
    const ids = Object.keys(a.scores).filter(id => id in b.scores)
    if (ids.length < 3) continue
    const method = options.correlationMethod === 'auto' ? (isNormal(a) && isNormal(b) ? 'pearson' : 'spearman') : options.correlationMethod
    const res = correlation(ids.map(id => a.scores[id]), ids.map(id => b.scores[id]), method)
    if (res) correlations.push({ a: a.scope.id, b: b.scope.id, method, r: res.r, pValue: res.pValue, n: res.n })
  }

  // Chi cuadrado: cada variable de grupo con los niveles de cada variable (o del total).
  const levelScopes = variables.length ? variables : scopes.filter(s => s.scope.kind === 'total')
  const crosstabs: CrosstabResult[] = []
  for (const g of input.groups) {
    for (const s of levelScopes) {
      const table = g.options.map(() => new Array<number>(levelCount).fill(0))
      let n = 0
      for (const r of input.respondents) {
        const opt = r.groups[g.id]
        const li = s.levelOf[r.id]
        const gi = opt === undefined ? -1 : g.options.indexOf(opt)
        if (gi === -1 || li === undefined) continue
        table[gi][li]++
        n++
      }
      // Sin filas o columnas vacías (la prueba no admite totales en cero).
      const keepRows = table.map((row, i) => [row, i] as const).filter(([row]) => row.some(v => v > 0))
      const keepCols = Array.from({ length: levelCount }, (_, j) => j).filter(j => keepRows.some(([row]) => row[j] > 0))
      const compact = keepRows.map(([row]) => keepCols.map(j => row[j]))
      const test = chiSquareTest(compact, options.confidence)
      if (!test) continue
      crosstabs.push({
        groupId: g.id, groupName: g.name, scopeId: s.scope.id, table, rows: g.options,
        chi2: test.chi2, df: test.df, pValue: test.pValue, cramerV: test.cramerV, lowExpectedPct: test.lowExpectedPct, n,
      })
    }
  }

  return { respondents: input.respondents.length, scopes, correlations, crosstabs, alphaLevel: alpha }
}

// Lo que viaja a la pantalla: sin los puntajes de cada persona.
export type ScopeView = Omit<ScopeResult, 'scores' | 'levelOf'>
export type AnalysisView = Omit<Analysis, 'scopes'> & { scopes: ScopeView[] }

export function analysisView(a: Analysis): AnalysisView {
  return { ...a, scopes: a.scopes.map(({ scores: _s, levelOf: _l, ...rest }) => rest) }
}
