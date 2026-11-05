import { chiSquareP, cronbachAlpha, fTestP, normalCdf, tCritical, tTestP } from './stats.js'

// Más calculadoras para tesis: omega de McDonald, CVC de Hernández-Nieto, kappa de
// Cohen y de Fleiss, t de Student y U de Mann-Whitney. Todo corre en el navegador.

const mean = (v: number[]) => v.reduce((a, b) => a + b, 0) / v.length
const sampleVar = (v: number[]) => {
  const m = mean(v)
  return v.reduce((s, x) => s + (x - m) ** 2, 0) / (v.length - 1)
}
const median = (v: number[]) => {
  const s = [...v].sort((a, b) => a - b)
  const mid = Math.floor(s.length / 2)
  return s.length % 2 ? s[mid] : (s[mid - 1] + s[mid]) / 2
}

// ---------------------------------------------------------------------------
// Omega de McDonald
// ---------------------------------------------------------------------------

function correlationMatrix(rows: number[][]): number[][] | null {
  const k = rows[0].length
  const cols = Array.from({ length: k }, (_, j) => rows.map(r => r[j]))
  const means = cols.map(mean)
  const sds = cols.map((c, j) => Math.sqrt(c.reduce((s, x) => s + (x - means[j]) ** 2, 0)))
  if (sds.some(s => s === 0)) return null
  return cols.map((a, i) => cols.map((b, j) => {
    if (i === j) return 1
    let s = 0
    for (let n = 0; n < a.length; n++) s += (a[n] - means[i]) * (b[n] - means[j])
    return s / (sds[i] * sds[j])
  }))
}

// Inversa por Gauss-Jordan; null si la matriz es singular.
function invert(m: number[][]): number[][] | null {
  const n = m.length
  const a = m.map((row, i) => [...row, ...row.map((_, j) => (i === j ? 1 : 0))])
  for (let c = 0; c < n; c++) {
    let piv = c
    for (let r = c + 1; r < n; r++) if (Math.abs(a[r][c]) > Math.abs(a[piv][c])) piv = r
    if (Math.abs(a[piv][c]) < 1e-12) return null
    ;[a[c], a[piv]] = [a[piv], a[c]]
    const d = a[c][c]
    for (let j = 0; j < 2 * n; j++) a[c][j] /= d
    for (let r = 0; r < n; r++) {
      if (r === c) continue
      const f = a[r][c]
      if (f) for (let j = 0; j < 2 * n; j++) a[r][j] -= f * a[c][j]
    }
  }
  return a.map(row => row.slice(n))
}

// Autovalor y autovector mayores de una matriz simétrica (método de Jacobi).
function topEigen(m: number[][]): { value: number; vector: number[] } {
  const n = m.length
  const a = m.map(r => [...r])
  const v: number[][] = a.map((_, i) => a.map((__, j) => (i === j ? 1 : 0)))
  for (let sweep = 0; sweep < 100; sweep++) {
    let off = 0
    for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) off += a[p][q] ** 2
    if (off < 1e-22) break
    for (let p = 0; p < n; p++) {
      for (let q = p + 1; q < n; q++) {
        if (Math.abs(a[p][q]) < 1e-15) continue
        const theta = (a[q][q] - a[p][p]) / (2 * a[p][q])
        const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1))
        const c = 1 / Math.sqrt(t * t + 1)
        const s = t * c
        for (let k = 0; k < n; k++) {
          const akp = a[k][p], akq = a[k][q]
          a[k][p] = c * akp - s * akq
          a[k][q] = s * akp + c * akq
        }
        for (let k = 0; k < n; k++) {
          const apk = a[p][k], aqk = a[q][k]
          a[p][k] = c * apk - s * aqk
          a[q][k] = s * apk + c * aqk
        }
        for (let k = 0; k < n; k++) {
          const vkp = v[k][p], vkq = v[k][q]
          v[k][p] = c * vkp - s * vkq
          v[k][q] = s * vkp + c * vkq
        }
      }
    }
  }
  let best = 0
  for (let i = 1; i < n; i++) if (a[i][i] > a[best][best]) best = i
  return { value: a[best][best], vector: v.map(r => r[best]) }
}

// Modelo de un factor por mínimos residuales (ejes principales iterados, que
// converge a la misma solución que "minres"). Devuelve las cargas estandarizadas.
export function oneFactorLoadings(r: number[][]): number[] {
  const k = r.length
  const inv = invert(r)
  let h = inv ? inv.map((row, i) => Math.min(0.995, Math.max(0.005, 1 - 1 / row[i]))) : r.map((row, i) => Math.max(...row.filter((_, j) => j !== i).map(Math.abs)))
  let loadings = new Array<number>(k).fill(0)
  for (let iter = 0; iter < 1000; iter++) {
    const reduced = r.map((row, i) => row.map((x, j) => (i === j ? h[i] : x)))
    const { value, vector } = topEigen(reduced)
    const scale = Math.sqrt(Math.max(value, 0))
    loadings = vector.map(x => x * scale)
    const next = loadings.map(l => Math.min(0.995, l * l))
    const diff = Math.max(...next.map((x, i) => Math.abs(x - h[i])))
    h = next
    if (diff < 1e-9) break
  }
  // El signo de un factor es arbitrario: que la mayoría de las cargas sean positivas.
  if (loadings.reduce((s, l) => s + l, 0) < 0) loadings = loadings.map(l => -l)
  return loadings
}

function omegaOf(loadings: number[]): number {
  const sum = loadings.reduce((s, l) => s + l, 0)
  const uniq = loadings.reduce((s, l) => s + (1 - l * l), 0)
  return (sum * sum) / (sum * sum + uniq)
}

// Filas = personas, columnas = ítems.
export function mcdonaldOmega(rows: number[][]) {
  if (rows.length < 3 || (rows[0]?.length ?? 0) < 3) return { error: 'size' as const }
  const r = correlationMatrix(rows)
  if (!r) return { error: 'variance' as const }
  const loadings = oneFactorLoadings(r)
  const omega = omegaOf(loadings)
  const alpha = cronbachAlpha(rows)
  const k = loadings.length
  const items = loadings.map((loading, j) => {
    let ifDeleted: number | null = null
    if (k > 3) {
      const keep = r.map((_, i) => i).filter(i => i !== j)
      ifDeleted = omegaOf(oneFactorLoadings(keep.map(a => keep.map(b => r[a][b]))))
    }
    return { loading, ifDeleted }
  })
  return { omega, alpha: alpha?.alpha ?? null, items, k, respondents: rows.length, negative: loadings.some(l => l < 0) }
}

// ---------------------------------------------------------------------------
// Coeficiente de validez de contenido (Hernández-Nieto, 2002)
// ---------------------------------------------------------------------------

// Filas = ítems, columnas = jueces, notas de 1 a vmax.
export function contentValidity(rows: number[][], vmax: number) {
  if (rows.length < 1 || (rows[0]?.length ?? 0) < 2) return { error: 'size' as const }
  if (rows.some(r => r.some(v => v < 1 || v > vmax || !Number.isInteger(v)))) return { error: 'range' as const }
  const judges = rows[0].length
  const pe = (1 / judges) ** judges
  const items = rows.map(r => {
    const avg = mean(r)
    const cvc = avg / vmax
    return { mean: avg, cvc, corrected: cvc - pe }
  })
  const total = mean(items.map(i => i.corrected))
  return { items, total, pe, judges }
}

// Hernández-Nieto: menos de 0,60 inaceptable, 0,60 a 0,70 deficiente, 0,71 a 0,80
// aceptable, 0,81 a 0,90 buena y más de 0,90 excelente.
export function cvcLevel(v: number): 'unacceptable' | 'deficient' | 'acceptable' | 'good' | 'excellent' {
  if (v > 0.9) return 'excellent'
  if (v > 0.8) return 'good'
  if (v > 0.7) return 'acceptable'
  if (v >= 0.6) return 'deficient'
  return 'unacceptable'
}

// ---------------------------------------------------------------------------
// Kappa de Cohen (2 evaluadores) y de Fleiss (3 o más)
// ---------------------------------------------------------------------------

// Tabla de categorías (texto o números): filas = casos, columnas = evaluadores.
export function parseCategoryTable(text: string): { rows: string[][]; skippedHeader: boolean } | { error: 'empty' | 'columns' } {
  const lines = text.split(/\r?\n/).filter(l => l.trim())
  if (!lines.length) return { error: 'empty' }
  const sep = lines[0].includes('\t') ? '\t' : lines[0].includes(';') ? ';' : lines[0].includes(',') ? ',' : /\s+/
  const rows = lines.map(l => l.trim().split(sep).map(c => c.trim()))
  let skippedHeader = false
  // Encabezado: primera fila con nombres que no se repiten en el resto de la tabla.
  if (rows.length > 2) {
    const body = new Set(rows.slice(1).flat().map(c => c.toLowerCase()))
    if (rows[0].every(c => c && !body.has(c.toLowerCase())) && rows[0].some(c => Number.isNaN(Number(c)))) {
      rows.shift()
      skippedHeader = true
    }
  }
  const width = rows[0].length
  if (rows.some(r => r.length !== width)) return { error: 'columns' }
  return { rows, skippedHeader }
}

export function kappaLevel(k: number): 'poor' | 'slight' | 'fair' | 'moderate' | 'substantial' | 'almostPerfect' {
  if (k < 0) return 'poor'
  if (k <= 0.2) return 'slight'
  if (k <= 0.4) return 'fair'
  if (k <= 0.6) return 'moderate'
  if (k <= 0.8) return 'substantial'
  return 'almostPerfect'
}

const categoriesOf = (rows: string[][]) => [...new Set(rows.flat())].sort((a, b) => {
  const na = Number(a), nb = Number(b)
  return Number.isFinite(na) && Number.isFinite(nb) ? na - nb : a.localeCompare(b)
})

export function cohenKappa(a: string[], b: string[]) {
  const n = a.length
  const cats = categoriesOf([a, b])
  const idx = new Map(cats.map((c, i) => [c, i]))
  const c = cats.length
  const p = Array.from({ length: c }, () => new Array<number>(c).fill(0))
  for (let i = 0; i < n; i++) p[idx.get(a[i])!][idx.get(b[i])!] += 1 / n
  const row = p.map(r => r.reduce((s, x) => s + x, 0))
  const col = cats.map((_, j) => p.reduce((s, r) => s + r[j], 0))
  const po = p.reduce((s, r, i) => s + r[i], 0)
  const pe = row.reduce((s, x, i) => s + x * col[i], 0)
  if (pe >= 1) return null
  const kappa = (po - pe) / (1 - pe)
  // Errores estándar de Fleiss, Cohen y Everitt (1969): bajo la hipótesis nula
  // (para el valor p) y sin suponerla (para el intervalo de confianza).
  const se0 = Math.sqrt((pe + pe * pe - row.reduce((s, x, i) => s + x * col[i] * (x + col[i]), 0)) / (n * (1 - pe) ** 2))
  let s1 = 0
  for (let i = 0; i < c; i++) {
    s1 += p[i][i] * (1 - (row[i] + col[i]) * (1 - kappa)) ** 2
    for (let j = 0; j < c; j++) if (i !== j) s1 += (1 - kappa) ** 2 * p[i][j] * (col[i] + row[j]) ** 2
  }
  s1 -= (kappa - pe * (1 - kappa)) ** 2
  const se = Math.sqrt(Math.max(s1, 0) / (n * (1 - pe) ** 2))
  const z = se0 > 0 ? kappa / se0 : 0
  const pValue = se0 > 0 ? 2 * (1 - normalCdf(Math.abs(z))) : 1
  return { kappa, po, pe, se, z, pValue, low: kappa - 1.96 * se, high: kappa + 1.96 * se, n, categories: cats }
}

export function fleissKappa(rows: string[][]) {
  const N = rows.length
  const m = rows[0].length
  const cats = categoriesOf(rows)
  const counts = rows.map(r => cats.map(c => r.filter(x => x === c).length))
  const pj = cats.map((_, j) => counts.reduce((s, r) => s + r[j], 0) / (N * m))
  const Pi = counts.map(r => (r.reduce((s, x) => s + x * x, 0) - m) / (m * (m - 1)))
  const P = mean(Pi)
  const Pe = pj.reduce((s, x) => s + x * x, 0)
  if (Pe >= 1) return null
  const kappa = (P - Pe) / (1 - Pe)
  const pq = pj.reduce((s, x) => s + x * (1 - x), 0)
  const se0 = (Math.sqrt(2) / (pq * Math.sqrt(N * m * (m - 1)))) * Math.sqrt(Math.max(pq * pq - pj.reduce((s, x) => s + x * (1 - x) * (1 - 2 * x), 0), 0))
  const z = se0 > 0 ? kappa / se0 : 0
  const pValue = se0 > 0 ? 2 * (1 - normalCdf(Math.abs(z))) : 1
  // Kappa de cada categoría (Fleiss, 1971).
  const perCategory = cats.map((cat, j) => {
    const q = 1 - pj[j]
    if (pj[j] === 0 || q === 0) return { category: cat, kappa: null as number | null }
    const k = 1 - counts.reduce((s, r) => s + r[j] * (m - r[j]), 0) / (N * m * (m - 1) * pj[j] * q)
    return { category: cat, kappa: k }
  })
  return { kappa, observed: P, expected: Pe, se0, z, pValue, n: N, raters: m, categories: cats, perCategory }
}

// Elige Cohen con 2 columnas y Fleiss con 3 o más. Las filas con celdas vacías se saltan.
export function interRaterAgreement(rows: string[][]) {
  const complete = rows.filter(r => r.every(c => c !== ''))
  const skipped = rows.length - complete.length
  if (complete.length < 2 || (complete[0]?.length ?? 0) < 2) return { error: 'size' as const }
  if (complete[0].length === 2) {
    const r = cohenKappa(complete.map(x => x[0]), complete.map(x => x[1]))
    if (!r) return { error: 'same' as const }
    return { method: 'cohen' as const, skipped, ...r }
  }
  const r = fleissKappa(complete)
  if (!r) return { error: 'same' as const }
  return { method: 'fleiss' as const, skipped, ...r }
}

// ---------------------------------------------------------------------------
// t de Student
// ---------------------------------------------------------------------------

const ci = (diff: number, se: number, df: number) => {
  const tc = tCritical(df)
  return { low: diff - tc * se, high: diff + tc * se }
}

// Levene centrado en la media, como SPSS.
export function leveneTest(a: number[], b: number[]) {
  const ma = mean(a), mb = mean(b)
  const za = a.map(x => Math.abs(x - ma))
  const zb = b.map(x => Math.abs(x - mb))
  const all = [...za, ...zb]
  const g = mean(all), ga = mean(za), gb = mean(zb)
  const N = all.length
  const between = a.length * (ga - g) ** 2 + b.length * (gb - g) ** 2
  const within = za.reduce((s, z) => s + (z - ga) ** 2, 0) + zb.reduce((s, z) => s + (z - gb) ** 2, 0)
  if (within === 0) return null
  const f = (between / 1) / (within / (N - 2))
  return { f, df1: 1, df2: N - 2, pValue: fTestP(f, 1, N - 2) }
}

export function independentTTest(a: number[], b: number[]) {
  if (a.length < 2 || b.length < 2) return null
  const n1 = a.length, n2 = b.length
  const m1 = mean(a), m2 = mean(b), v1 = sampleVar(a), v2 = sampleVar(b)
  if (v1 === 0 && v2 === 0) return null
  const diff = m1 - m2
  const df = n1 + n2 - 2
  const pooled = ((n1 - 1) * v1 + (n2 - 1) * v2) / df
  const seStudent = Math.sqrt(pooled * (1 / n1 + 1 / n2))
  const tStudent = diff / seStudent
  const seWelch = Math.sqrt(v1 / n1 + v2 / n2)
  const tWelch = diff / seWelch
  const dfWelch = (v1 / n1 + v2 / n2) ** 2 / ((v1 / n1) ** 2 / (n1 - 1) + (v2 / n2) ** 2 / (n2 - 1))
  const levene = leveneTest(a, b)
  return {
    a: { n: n1, mean: m1, sd: Math.sqrt(v1) },
    b: { n: n2, mean: m2, sd: Math.sqrt(v2) },
    diff,
    student: { t: tStudent, df, pValue: tTestP(tStudent, df), ...ci(diff, seStudent, df) },
    welch: { t: tWelch, df: dfWelch, pValue: tTestP(tWelch, dfWelch), ...ci(diff, seWelch, dfWelch) },
    levene,
    // Si Levene sale significativo, las varianzas no son iguales y se reporta Welch.
    recommended: levene && levene.pValue < 0.05 ? ('welch' as const) : ('student' as const),
    cohenD: diff / Math.sqrt(pooled),
  }
}

export function pairedTTest(a: number[], b: number[]) {
  if (a.length !== b.length || a.length < 2) return null
  const d = a.map((x, i) => x - b[i])
  const n = d.length
  const md = mean(d)
  const sd = Math.sqrt(sampleVar(d))
  if (sd === 0) return null
  const se = sd / Math.sqrt(n)
  const t = md / se
  const df = n - 1
  return {
    a: { n, mean: mean(a), sd: Math.sqrt(sampleVar(a)) },
    b: { n, mean: mean(b), sd: Math.sqrt(sampleVar(b)) },
    diff: md, sdDiff: sd, t, df, pValue: tTestP(t, df), ...ci(md, se, df),
    cohenD: md / sd,
  }
}

// Cohen: 0,2 pequeño, 0,5 mediano, 0,8 grande.
export function effectSizeLevel(d: number): 'negligible' | 'small' | 'medium' | 'large' {
  const x = Math.abs(d)
  if (x >= 0.8) return 'large'
  if (x >= 0.5) return 'medium'
  if (x >= 0.2) return 'small'
  return 'negligible'
}

// ---------------------------------------------------------------------------
// U de Mann-Whitney
// ---------------------------------------------------------------------------

// Distribución exacta de U sin empates (conteo de permutaciones por recurrencia).
function exactUPValue(u: number, n1: number, n2: number): number {
  // f[i][j][k] = maneras de obtener U = k con i y j elementos; se guarda por capas.
  const max = n1 * n2
  let prev: number[][] = Array.from({ length: n2 + 1 }, () => [1])
  for (let i = 1; i <= n1; i++) {
    const cur: number[][] = [[1]]
    for (let j = 1; j <= n2; j++) {
      const size = i * j + 1
      const row = new Array<number>(size).fill(0)
      const a = cur[j - 1], b = prev[j]
      // U(i, j) = U(i, j-1) [el último es de B] o U(i-1, j) + j [el último es de A].
      for (let k = 0; k < a.length; k++) row[k] += a[k]
      for (let k = 0; k < b.length; k++) row[k + j] += b[k]
      cur.push(row)
    }
    prev = cur
  }
  const dist = prev[n2]
  const total = dist.reduce((s, x) => s + x, 0)
  let tail = 0
  for (let k = 0; k <= Math.min(u, max); k++) tail += dist[k]
  return Math.min(1, (2 * tail) / total)
}

export function mannWhitney(a: number[], b: number[]) {
  const n1 = a.length, n2 = b.length
  if (n1 < 2 || n2 < 2) return null
  const all = [...a.map(v => ({ v, g: 0 })), ...b.map(v => ({ v, g: 1 }))].sort((x, y) => x.v - y.v)
  const N = all.length
  const rank = new Array<number>(N)
  let ties = 0
  let tieGroups = 0
  for (let i = 0; i < N;) {
    let j = i
    while (j + 1 < N && all[j + 1].v === all[i].v) j++
    const r = (i + j) / 2 + 1
    for (let k = i; k <= j; k++) rank[k] = r
    const t = j - i + 1
    if (t > 1) { ties += t ** 3 - t; tieGroups++ }
    i = j + 1
  }
  const r1 = all.reduce((s, x, i) => s + (x.g === 0 ? rank[i] : 0), 0)
  const r2 = N * (N + 1) / 2 - r1
  const u1 = r1 - n1 * (n1 + 1) / 2
  const u2 = n1 * n2 - u1
  const u = Math.min(u1, u2)
  const sigma = Math.sqrt((n1 * n2 / 12) * ((N + 1) - ties / (N * (N - 1))))
  if (sigma === 0) return null
  const z = (u - n1 * n2 / 2) / sigma
  const pValue = Math.min(1, 2 * (1 - normalCdf(Math.abs(z))))
  // Exacta solo sin empates y con muestras pequeñas (como el "Sig. exacta" de SPSS).
  const exact = tieGroups === 0 && n1 * n2 <= 1600 ? exactUPValue(u, n1, n2) : null
  return {
    a: { n: n1, median: median(a), meanRank: r1 / n1, rankSum: r1 },
    b: { n: n2, median: median(b), meanRank: r2 / n2, rankSum: r2 },
    u, w: u === u1 ? r1 : r2, z, pValue, exact, ties: tieGroups > 0,
    r: Math.abs(z) / Math.sqrt(N),
  }
}

// Tamaño del efecto r = |Z| / √N (Cohen: 0,1 pequeño, 0,3 mediano, 0,5 grande).
export function rEffectLevel(r: number): 'negligible' | 'small' | 'medium' | 'large' {
  if (r >= 0.5) return 'large'
  if (r >= 0.3) return 'medium'
  if (r >= 0.1) return 'small'
  return 'negligible'
}

// ---------------------------------------------------------------------------
// Wilcoxon de rangos con signo (muestras relacionadas)
// ---------------------------------------------------------------------------

// Rangos promedio de una lista; devuelve también la suma de t³ − t de los empates.
function rankWithTies(values: number[]): { ranks: number[]; tieSum: number } {
  const order = values.map((v, i) => [v, i] as const).sort((a, b) => a[0] - b[0])
  const ranks = new Array<number>(values.length)
  let tieSum = 0
  for (let i = 0; i < order.length;) {
    let j = i
    while (j + 1 < order.length && order[j + 1][0] === order[i][0]) j++
    for (let k = i; k <= j; k++) ranks[order[k][1]] = (i + j) / 2 + 1
    const t = j - i + 1
    if (t > 1) tieSum += t ** 3 - t
    i = j + 1
  }
  return { ranks, tieSum }
}

// p exacta bilateral de T (suma de rangos positivos) sin empates: cuenta los 2^n signos posibles.
function exactSignedRankP(t: number, n: number): number {
  const max = (n * (n + 1)) / 2
  const dist = new Array<number>(max + 1).fill(0)
  dist[0] = 1
  for (let k = 1; k <= n; k++) for (let s = max; s >= k; s--) dist[s] += dist[s - k]
  const total = 2 ** n
  let tail = 0
  for (let s = 0; s <= Math.min(t, max); s++) tail += dist[s]
  return Math.min(1, (2 * tail) / total)
}

// Diferencias = segunda medición − primera (como SPSS). Las diferencias en cero se descartan.
export function wilcoxonSignedRank(before: number[], after: number[]) {
  if (before.length !== after.length || before.length < 2) return null
  // Redondeo para que 4,9 − 4,1 y 3,9 − 3,1 cuenten como la misma diferencia.
  const diffs = after.map((x, i) => Math.round((x - before[i]) * 1e9) / 1e9)
  const nonZero = diffs.filter(d => d !== 0)
  const n = nonZero.length
  if (n < 1) return null
  const { ranks, tieSum } = rankWithTies(nonZero.map(Math.abs))
  let wPlus = 0, wMinus = 0, pos = 0, neg = 0
  nonZero.forEach((d, i) => {
    if (d > 0) { wPlus += ranks[i]; pos++ } else { wMinus += ranks[i]; neg++ }
  })
  const t = Math.min(wPlus, wMinus)
  const sigma = Math.sqrt((n * (n + 1) * (2 * n + 1)) / 24 - tieSum / 48)
  const z = sigma > 0 ? (t - (n * (n + 1)) / 4) / sigma : 0
  const pValue = sigma > 0 ? Math.min(1, 2 * (1 - normalCdf(Math.abs(z)))) : 1
  const zeros = diffs.length - n
  const exact = tieSum === 0 && zeros === 0 && n <= 50 ? exactSignedRankP(t, n) : null
  return {
    n: before.length, used: n, zeros, positive: pos, negative: neg,
    meanRankPositive: pos ? wPlus / pos : 0, meanRankNegative: neg ? wMinus / neg : 0,
    wPlus, wMinus, t, z, pValue, exact, ties: tieSum > 0,
    medianBefore: median(before), medianAfter: median(after),
    r: Math.abs(z) / Math.sqrt(n),
  }
}

// ---------------------------------------------------------------------------
// Kruskal-Wallis (3 o más grupos independientes)
// ---------------------------------------------------------------------------

export function kruskalWallis(groups: number[][]) {
  const valid = groups.filter(g => g.length > 0)
  if (valid.length < 2 || valid.some(g => g.length < 2)) return null
  const all = valid.flat()
  const N = all.length
  const { ranks, tieSum } = rankWithTies(all)
  let offset = 0
  const stats = valid.map(g => {
    const r = ranks.slice(offset, offset + g.length)
    offset += g.length
    const sum = r.reduce((a, b) => a + b, 0)
    return { n: g.length, median: median(g), rankSum: sum, meanRank: sum / g.length }
  })
  const correction = 1 - tieSum / (N ** 3 - N)
  if (correction <= 0) return null
  const h = ((12 / (N * (N + 1))) * stats.reduce((s, g) => s + g.rankSum ** 2 / g.n, 0) - 3 * (N + 1)) / correction
  const df = valid.length - 1
  // Comparaciones por pares de Dunn con corrección de Bonferroni (como SPSS).
  const pairs: { a: number; b: number; z: number; p: number; adjusted: number }[] = []
  const m = (valid.length * (valid.length - 1)) / 2
  const base = (N * (N + 1)) / 12 - tieSum / (12 * (N - 1))
  for (let i = 0; i < stats.length; i++) {
    for (let j = i + 1; j < stats.length; j++) {
      const se = Math.sqrt(base * (1 / stats[i].n + 1 / stats[j].n))
      const z = se > 0 ? (stats[i].meanRank - stats[j].meanRank) / se : 0
      const p = Math.min(1, 2 * (1 - normalCdf(Math.abs(z))))
      pairs.push({ a: i, b: j, z, p, adjusted: Math.min(1, p * m) })
    }
  }
  return { h, df, pValue: chiSquareP(h, df), n: N, groups: stats, epsilon2: h / (N - 1), pairs, ties: tieSum > 0 }
}

// ---------------------------------------------------------------------------
// ANOVA de un factor
// ---------------------------------------------------------------------------

export function oneWayAnova(groups: number[][]) {
  const valid = groups.filter(g => g.length > 0)
  const k = valid.length
  if (k < 2 || valid.some(g => g.length < 2)) return null
  const N = valid.reduce((s, g) => s + g.length, 0)
  const grand = mean(valid.flat())
  const stats = valid.map(g => ({ n: g.length, mean: mean(g), sd: Math.sqrt(sampleVar(g)), variance: sampleVar(g) }))
  const ssBetween = stats.reduce((s, g) => s + g.n * (g.mean - grand) ** 2, 0)
  const ssWithin = valid.reduce((s, g, i) => s + g.reduce((a, x) => a + (x - stats[i].mean) ** 2, 0), 0)
  if (ssWithin === 0) return null
  const dfB = k - 1, dfW = N - k
  const msB = ssBetween / dfB, msW = ssWithin / dfW
  const f = msB / msW

  // Levene centrado en la media para k grupos.
  const z = valid.map((g, i) => g.map(x => Math.abs(x - stats[i].mean)))
  const zMeans = z.map(mean)
  const zGrand = mean(z.flat())
  const levB = z.reduce((s, g, i) => s + g.length * (zMeans[i] - zGrand) ** 2, 0)
  const levW = z.reduce((s, g, i) => s + g.reduce((a, x) => a + (x - zMeans[i]) ** 2, 0), 0)
  const leveneF = levW > 0 ? (levB / dfB) / (levW / dfW) : null
  const levene = leveneF != null ? { f: leveneF, df1: dfB, df2: dfW, pValue: fTestP(leveneF, dfB, dfW) } : null

  // ANOVA de Welch, para cuando las varianzas no son iguales.
  let welch: { f: number; df1: number; df2: number; pValue: number } | null = null
  if (stats.every(g => g.variance > 0)) {
    const w = stats.map(g => g.n / g.variance)
    const W = w.reduce((a, b) => a + b, 0)
    const mw = stats.reduce((s, g, i) => s + w[i] * g.mean, 0) / W
    const a = stats.reduce((s, g, i) => s + w[i] * (g.mean - mw) ** 2, 0) / (k - 1)
    const lambda = stats.reduce((s, g, i) => s + (1 - w[i] / W) ** 2 / (g.n - 1), 0)
    const b = 1 + ((2 * (k - 2)) / (k * k - 1)) * lambda
    const wf = a / b
    const df2 = (k * k - 1) / (3 * lambda)
    welch = { f: wf, df1: k - 1, df2, pValue: fTestP(wf, k - 1, df2) }
  }

  // Comparaciones por pares con t y la varianza combinada, corregidas por Bonferroni.
  const m = (k * (k - 1)) / 2
  const pairs: { a: number; b: number; diff: number; p: number; adjusted: number }[] = []
  for (let i = 0; i < k; i++) {
    for (let j = i + 1; j < k; j++) {
      const diff = stats[i].mean - stats[j].mean
      const se = Math.sqrt(msW * (1 / stats[i].n + 1 / stats[j].n))
      const p = se > 0 ? tTestP(diff / se, dfW) : 1
      pairs.push({ a: i, b: j, diff, p, adjusted: Math.min(1, p * m) })
    }
  }

  return {
    groups: stats, ssBetween, ssWithin, dfB, dfW, msB, msW, f, pValue: fTestP(f, dfB, dfW),
    eta2: ssBetween / (ssBetween + ssWithin), levene, welch,
    recommended: levene && levene.pValue < 0.05 && welch ? ('welch' as const) : ('classic' as const),
    pairs, n: N,
  }
}

// Eta cuadrado (Cohen): 0,01 pequeño, 0,06 mediano, 0,14 grande. Sirve también para épsilon².
export function eta2Level(e: number): 'negligible' | 'small' | 'medium' | 'large' {
  if (e >= 0.14) return 'large'
  if (e >= 0.06) return 'medium'
  if (e >= 0.01) return 'small'
  return 'negligible'
}

// ---------------------------------------------------------------------------
// W de Kendall (concordancia entre jueces)
// ---------------------------------------------------------------------------

// Filas = ítems (lo que se ordena), columnas = jueces. Cada juez convierte sus
// notas en rangos (con empates promediados), como SPSS.
export function kendallW(rows: number[][]) {
  const n = rows.length
  const k = rows[0]?.length ?? 0
  if (n < 3 || k < 2) return null
  const rankSums = new Array<number>(n).fill(0)
  let tieTotal = 0
  for (let j = 0; j < k; j++) {
    const { ranks, tieSum } = rankWithTies(rows.map(r => r[j]))
    ranks.forEach((r, i) => { rankSums[i] += r })
    tieTotal += tieSum
  }
  const avg = rankSums.reduce((a, b) => a + b, 0) / n
  const s = rankSums.reduce((acc, r) => acc + (r - avg) ** 2, 0)
  const denom = k ** 2 * (n ** 3 - n) - k * tieTotal
  if (denom <= 0) return null
  const w = (12 * s) / denom
  const chi2 = k * (n - 1) * w
  const df = n - 1
  return {
    w, chi2, df, pValue: chiSquareP(chi2, df), items: n, judges: k, s,
    meanRanks: rankSums.map(r => r / k), ties: tieTotal > 0,
  }
}

// Escala de Landis y Koch, la que suelen citar las tesis para W.
export function kendallLevel(w: number): 'poor' | 'slight' | 'fair' | 'moderate' | 'substantial' | 'almostPerfect' {
  if (w <= 0) return 'poor'
  if (w <= 0.2) return 'slight'
  if (w <= 0.4) return 'fair'
  if (w <= 0.6) return 'moderate'
  if (w <= 0.8) return 'substantial'
  return 'almostPerfect'
}

// ---------------------------------------------------------------------------
// Prueba binomial para juicio de expertos
// ---------------------------------------------------------------------------

function binomCoef(n: number, k: number): number {
  let c = 1
  for (let i = 1; i <= k; i++) c = (c * (n - k + i)) / i
  return c
}

// P(X >= x) con X ~ Binomial(n, 0,5).
export function binomialUpperP(x: number, n: number): number {
  let p = 0
  for (let i = x; i <= n; i++) p += binomCoef(n, i)
  return Math.min(1, p / 2 ** n)
}

// p bilateral exacta con proporción de prueba 0,5 (la que da SPSS).
export function binomialTwoSidedP(x: number, n: number): number {
  const tail = Math.min(x, n - x)
  let p = 0
  for (let i = 0; i <= tail; i++) p += binomCoef(n, i)
  return Math.min(1, (2 * p) / 2 ** n)
}

// Filas = ítems, columnas = jueces; 1 = de acuerdo (sí), 0 = en desacuerdo (no).
export function binomialJudges(rows: number[][]) {
  if (rows.length < 1 || (rows[0]?.length ?? 0) < 2) return { error: 'size' as const }
  if (rows.some(r => r.some(v => v !== 0 && v !== 1))) return { error: 'binary' as const }
  const judges = rows[0].length
  const items = rows.map(r => {
    const agree = r.reduce((a, b) => a + b, 0)
    return { agree, pct: agree / judges, p: binomialUpperP(agree, judges), pTwo: binomialTwoSidedP(agree, judges) }
  })
  const totalAgree = items.reduce((a, i) => a + i.agree, 0)
  const total = rows.length * judges
  return {
    judges,
    items,
    meanP: items.reduce((a, i) => a + i.p, 0) / items.length,
    totalAgree, total, totalPct: totalAgree / total,
    totalP: binomialUpperP(totalAgree, total),
  }
}
