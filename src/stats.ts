// Cálculos de las herramientas gratis (/herramientas, /tools). Todo corre en el
// navegador y no se guarda nada. Las fórmulas son las estándar de libro:
// muestra y margen de error para una proporción (p = 0,5 si no se sabe), NPS,
// CSAT, prueba z de dos proporciones y alfa de Cronbach con varianza muestral.

export const Z_BY_CONFIDENCE: Record<string, number> = { '90': 1.645, '95': 1.96, '99': 2.576 }

const zFor = (confidence: string) => Z_BY_CONFIDENCE[confidence] ?? 1.96

// Respuestas necesarias para un margen de error dado. Con población conocida se
// aplica la corrección por población finita.
export function sampleSize(population: number | null, confidence: string, marginPct: number, p = 0.5): number | null {
  if (!(marginPct > 0) || !(p > 0 && p < 1)) return null
  if (population != null && !(population >= 1)) return null
  const z = zFor(confidence)
  const e = marginPct / 100
  const n0 = (z * z * p * (1 - p)) / (e * e)
  if (population == null) return Math.ceil(n0)
  return Math.min(Math.ceil(population), Math.ceil(n0 / (1 + (n0 - 1) / population)))
}

// Invitaciones a enviar para llegar a `needed` respuestas con la tasa esperada.
export function invitationsNeeded(needed: number, responseRatePct: number): number | null {
  if (!(responseRatePct > 0 && responseRatePct <= 100)) return null
  return Math.ceil(needed / (responseRatePct / 100))
}

// Margen de error (en puntos porcentuales) de una muestra ya recogida.
export function marginOfError(n: number, population: number | null, confidence: string, p = 0.5): number | null {
  if (!(n >= 1) || !(p > 0 && p < 1)) return null
  if (population != null && (!(population >= 1) || n > population)) return null
  const z = zFor(confidence)
  let e = z * Math.sqrt((p * (1 - p)) / n)
  if (population != null && population > 1) e *= Math.sqrt((population - n) / (population - 1))
  return e * 100
}

export function npsScore(promoters: number, passives: number, detractors: number) {
  const total = promoters + passives + detractors
  if (!(total > 0) || [promoters, passives, detractors].some(v => v < 0)) return null
  return {
    total,
    score: Math.round(((promoters - detractors) / total) * 1000) / 10,
    promotersPct: (promoters / total) * 100,
    passivesPct: (passives / total) * 100,
    detractorsPct: (detractors / total) * 100,
  }
}

// counts[0] = cuántos marcaron 1 ... counts[4] = cuántos marcaron 5.
export function csatScore(counts: number[]) {
  const total = counts.reduce((a, b) => a + b, 0)
  if (!(total > 0) || counts.some(v => v < 0)) return null
  const satisfied = (counts[3] ?? 0) + (counts[4] ?? 0)
  const average = counts.reduce((sum, c, i) => sum + c * (i + 1), 0) / total
  return { total, csatPct: (satisfied / total) * 100, average }
}

// CES (Customer Effort Score). counts[0] = cuántos marcaron 1 ... counts[k-1] =
// cuántos marcaron el valor más alto. Escala de 1 a 7 (la más usada) o de 1 a 5,
// donde el valor alto significa "me fue fácil". En 1-7 cuentan como fácil 5, 6 y
// 7 y como difícil 1, 2 y 3; en 1-5, 4 y 5 contra 1 y 2.
export function cesScore(counts: number[]) {
  const total = counts.reduce((a, b) => a + b, 0)
  const k = counts.length
  if (!(total > 0) || counts.some(v => v < 0) || (k !== 5 && k !== 7)) return null
  const top = k === 7 ? 3 : 2
  const easy = counts.slice(k - top).reduce((a, b) => a + b, 0)
  const hard = counts.slice(0, top).reduce((a, b) => a + b, 0)
  const average = counts.reduce((sum, c, i) => sum + c * (i + 1), 0) / total
  return { total, average, easyPct: (easy / total) * 100, hardPct: (hard / total) * 100 }
}

export function responseRates(sent: number, started: number | null, completed: number) {
  if (!(sent > 0) || completed < 0 || completed > sent) return null
  if (started != null && (started < completed || started > sent)) return null
  return {
    responseRatePct: (completed / sent) * 100,
    completionRatePct: started != null && started > 0 ? (completed / started) * 100 : null,
    startRatePct: started != null ? (started / sent) * 100 : null,
  }
}

// Función de distribución normal estándar (aproximación de Abramowitz y Stegun 7.1.26).
export function normalCdf(x: number): number {
  const t = 1 / (1 + 0.3275911 * Math.abs(x) / Math.SQRT2)
  const y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-(x * x) / 2)
  return x >= 0 ? (1 + y) / 2 : (1 - y) / 2
}

// ¿La diferencia entre dos porcentajes es real o puede ser azar? Prueba z de dos
// proporciones (bilateral, con proporción combinada).
export function twoProportionTest(x1: number, n1: number, x2: number, n2: number, confidence: string) {
  if (!(n1 >= 1 && n2 >= 1) || x1 < 0 || x2 < 0 || x1 > n1 || x2 > n2) return null
  const p1 = x1 / n1
  const p2 = x2 / n2
  const pooled = (x1 + x2) / (n1 + n2)
  const se = Math.sqrt(pooled * (1 - pooled) * (1 / n1 + 1 / n2))
  if (se === 0) return { p1: p1 * 100, p2: p2 * 100, diff: (p1 - p2) * 100, z: 0, pValue: 1, significant: false }
  const z = (p1 - p2) / se
  const pValue = Math.min(1, 2 * (1 - normalCdf(Math.abs(z))))
  const alpha = 1 - Number(confidence) / 100
  return { p1: p1 * 100, p2: p2 * 100, diff: (p1 - p2) * 100, z, pValue, significant: pValue < alpha }
}

// Tabla pegada desde Excel o Google Sheets: una fila por persona, una columna por
// ítem. Acepta tabulaciones, punto y coma, comas o espacios, y salta una fila de
// encabezados si la hay.
export function parseScoreTable(text: string): { rows: number[][]; skippedHeader: boolean } | { error: 'empty' | 'columns' | 'notNumber' } {
  const lines = text.split(/\r?\n/).map(l => l.trim()).filter(Boolean)
  if (lines.length === 0) return { error: 'empty' }
  const sample = lines[0]
  const sep = sample.includes('\t') ? /\t/ : sample.includes(';') ? /;/ : sample.includes(',') && !/\d,\d+\s*$/.test(sample.split(',').pop() ?? '') ? /,/ : /\s+/
  const decimalComma = sep.source !== ','
  const parseCell = (c: string) => Number((decimalComma ? c.replace(',', '.') : c).trim())
  let skippedHeader = false
  const rows: number[][] = []
  for (const [i, line] of lines.entries()) {
    const cells = line.split(sep).map(c => c.trim()).filter(c => c !== '')
    const values = cells.map(parseCell)
    if (values.some(v => Number.isNaN(v))) {
      if (i === 0) { skippedHeader = true; continue }
      return { error: 'notNumber' }
    }
    rows.push(values)
  }
  if (rows.length === 0) return { error: 'empty' }
  const width = rows[0].length
  if (rows.some(r => r.length !== width)) return { error: 'columns' }
  return { rows, skippedHeader }
}

const sampleVariance = (values: number[]) => {
  const mean = values.reduce((a, b) => a + b, 0) / values.length
  return values.reduce((sum, v) => sum + (v - mean) ** 2, 0) / (values.length - 1)
}

function alphaOf(rows: number[][], items: number[]): number | null {
  const k = items.length
  if (k < 2) return null
  const itemVar = items.reduce((sum, j) => sum + sampleVariance(rows.map(r => r[j])), 0)
  const totalVar = sampleVariance(rows.map(r => items.reduce((s, j) => s + r[j], 0)))
  if (totalVar === 0) return null
  return (k / (k - 1)) * (1 - itemVar / totalVar)
}

export function cronbachAlpha(rows: number[][]) {
  if (rows.length < 2 || (rows[0]?.length ?? 0) < 2) return null
  const items = rows[0].map((_, j) => j)
  const alpha = alphaOf(rows, items)
  if (alpha == null) return null
  const ifDeleted = rows[0].length > 2 ? items.map(j => alphaOf(rows, items.filter(i => i !== j))) : []
  return { alpha, items: items.length, respondents: rows.length, ifDeleted }
}

// Escala habitual (George y Mallery) para leer el alfa.
export function alphaLevel(alpha: number): 'excellent' | 'good' | 'acceptable' | 'questionable' | 'poor' | 'unacceptable' {
  if (alpha >= 0.9) return 'excellent'
  if (alpha >= 0.8) return 'good'
  if (alpha >= 0.7) return 'acceptable'
  if (alpha >= 0.6) return 'questionable'
  if (alpha >= 0.5) return 'poor'
  return 'unacceptable'
}

// ---------------------------------------------------------------------------
// Herramientas para tesis: V de Aiken, baremos, chi-cuadrado y correlación.
// ---------------------------------------------------------------------------

// Logaritmo de la función gamma (aproximación de Lanczos).
function logGamma(x: number): number {
  const c = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5]
  let y = x
  const tmp = x + 5.5 - (x + 0.5) * Math.log(x + 5.5)
  let ser = 1.000000000190015
  for (const ci of c) ser += ci / ++y
  return -tmp + Math.log((2.5066282746310005 * ser) / x)
}

// Gamma incompleta regularizada superior Q(a, x) (Numerical Recipes).
function gammaQ(a: number, x: number): number {
  if (x <= 0) return 1
  const gln = logGamma(a)
  if (x < a + 1) {
    let sum = 1 / a
    let del = sum
    let ap = a
    for (let n = 0; n < 500; n++) {
      ap += 1
      del *= x / ap
      sum += del
      if (Math.abs(del) < Math.abs(sum) * 1e-14) break
    }
    return 1 - sum * Math.exp(-x + a * Math.log(x) - gln)
  }
  let b = x + 1 - a
  let c = 1 / 1e-300
  let d = 1 / b
  let h = d
  for (let i = 1; i < 500; i++) {
    const an = -i * (i - a)
    b += 2
    d = an * d + b
    if (Math.abs(d) < 1e-300) d = 1e-300
    c = b + an / c
    if (Math.abs(c) < 1e-300) c = 1e-300
    d = 1 / d
    const del = d * c
    h *= del
    if (Math.abs(del - 1) < 1e-14) break
  }
  return Math.exp(-x + a * Math.log(x) - gln) * h
}

// Beta incompleta regularizada I_x(a, b) (fracción continua de Numerical Recipes).
function betaI(a: number, b: number, x: number): number {
  if (x <= 0) return 0
  if (x >= 1) return 1
  const front = Math.exp(logGamma(a + b) - logGamma(a) - logGamma(b) + a * Math.log(x) + b * Math.log(1 - x))
  const cf = (aa: number, bb: number, xx: number) => {
    let c = 1
    let d = 1 - ((aa + bb) * xx) / (aa + 1)
    if (Math.abs(d) < 1e-300) d = 1e-300
    d = 1 / d
    let h = d
    for (let m = 1; m < 500; m++) {
      const m2 = 2 * m
      let num = (m * (bb - m) * xx) / ((aa + m2 - 1) * (aa + m2))
      d = 1 + num * d; if (Math.abs(d) < 1e-300) d = 1e-300
      c = 1 + num / c; if (Math.abs(c) < 1e-300) c = 1e-300
      d = 1 / d; h *= d * c
      num = (-(aa + m) * (aa + bb + m) * xx) / ((aa + m2) * (aa + m2 + 1))
      d = 1 + num * d; if (Math.abs(d) < 1e-300) d = 1e-300
      c = 1 + num / c; if (Math.abs(c) < 1e-300) c = 1e-300
      d = 1 / d
      const del = d * c
      h *= del
      if (Math.abs(del - 1) < 1e-14) break
    }
    return h
  }
  return x < (a + 1) / (a + b + 2) ? (front * cf(a, b, x)) / a : 1 - (front * cf(b, a, 1 - x)) / b
}

// p bilateral de una t de Student con df grados de libertad.
export function tTestP(t: number, df: number): number {
  return Math.min(1, betaI(df / 2, 0.5, df / (df + t * t)))
}

// p de una F de Snedecor con (df1, df2) grados de libertad (cola superior).
export function fTestP(f: number, df1: number, df2: number): number {
  if (!(f > 0)) return 1
  return Math.min(1, Math.max(0, betaI(df2 / 2, df1 / 2, df2 / (df2 + df1 * f))))
}

// p de una chi-cuadrado con df grados de libertad (cola superior).
export function chiSquareP(chi2: number, df: number): number {
  return Math.min(1, Math.max(0, gammaQ(df / 2, chi2 / 2)))
}

// Tabla pegada con etiquetas opcionales: salta una fila de encabezados y una
// primera columna de nombres si no son números.
export function parseLabeledTable(text: string):
  | { rows: number[][]; rowLabels: string[] | null; colLabels: string[] | null }
  | { error: 'empty' | 'columns' | 'notNumber' } {
  const lines = text.split(/\r?\n/).map(l => l.trim()).filter(Boolean)
  if (lines.length === 0) return { error: 'empty' }
  const first = lines[0]
  const sep = first.includes('\t') ? /\t/ : first.includes(';') ? /;/ : first.includes(',') && !/^\s*\d+,\d+\s*$/.test(first) ? /,/ : /\s+/
  const decimalComma = sep.source !== ','
  const toNum = (c: string) => (c === '' ? NaN : Number((decimalComma ? c.replace(',', '.') : c).trim()))
  const cells = lines.map(l => l.split(sep).map(c => c.trim()))
  const dataStart = cells[0].some(c => c !== '' && Number.isNaN(toNum(c))) && cells.length > 1 ? 1 : 0
  const body = cells.slice(dataStart)
  // Primera columna de nombres, o de números de ítem si el encabezado lo dice.
  const idHeader = dataStart === 1 && /^(ítem|item|items|ítems|n|nº|n°|no\.?|#|pregunta|reactivo|id)$/i.test(cells[0][0] ?? '')
  const hasRowLabels = idHeader || body.every(r => Number.isNaN(toNum(r[0] ?? '')))
  const rows: number[][] = []
  for (const r of body) {
    const values = (hasRowLabels ? r.slice(1) : r).filter(c => c !== '').map(toNum)
    if (values.length === 0 || values.some(v => Number.isNaN(v))) return { error: 'notNumber' }
    rows.push(values)
  }
  if (rows.length === 0) return { error: 'empty' }
  if (rows.some(r => r.length !== rows[0].length)) return { error: 'columns' }
  let colLabels: string[] | null = null
  if (dataStart === 1) {
    const head = cells[0].filter(c => c !== '')
    colLabels = head.length === rows[0].length + 1 && hasRowLabels ? head.slice(1) : head.slice(-rows[0].length)
  }
  return { rows, rowLabels: hasRowLabels ? body.map(r => r[0]) : null, colLabels }
}

// V de Aiken por ítem. ratings: una fila por ítem, una columna por juez, con
// valores entre `low` y `high`. Intervalo de confianza por el método score de
// Penfield y Giacobbi (2004).
export function aikenV(ratings: number[], low: number, high: number, confidence = '95') {
  const n = ratings.length
  const k = high - low
  if (n < 1 || !(k > 0) || ratings.some(r => r < low || r > high)) return null
  const s = ratings.reduce((sum, r) => sum + (r - low), 0)
  const v = s / (n * k)
  const z = zFor(confidence)
  const nk = n * k
  const root = z * Math.sqrt(4 * nk * v * (1 - v) + z * z)
  const den = 2 * (nk + z * z)
  return { v, lower: (2 * nk * v + z * z - root) / den, upper: (2 * nk * v + z * z + root) / den, judges: n }
}

// Baremos: rangos de igual amplitud para el puntaje total y el promedio.
export function likertLevels(items: number, low: number, high: number, levels: number) {
  if (!(items >= 1) || !Number.isInteger(items) || !(high > low) || !Number.isInteger(low) || !Number.isInteger(high) || ![2, 3, 4, 5].includes(levels)) return null
  const minT = items * low
  const maxT = items * high
  if (maxT - minT + 1 < levels) return null
  const width = (maxT - minT) / levels
  const avgWidth = (high - low) / levels
  return Array.from({ length: levels }, (_, i) => {
    const from = i === 0 ? minT : Math.floor(minT + i * width) + 1
    const to = i === levels - 1 ? maxT : Math.floor(minT + (i + 1) * width)
    return { from, to, avgFrom: low + i * avgWidth, avgTo: low + (i + 1) * avgWidth }
  })
}

// Prueba chi-cuadrado de independencia sobre una tabla de frecuencias.
export function chiSquareTest(table: number[][], confidence = '95') {
  const r = table.length
  const c = table[0]?.length ?? 0
  if (r < 2 || c < 2 || table.some(row => row.length !== c || row.some(v => v < 0 || !Number.isFinite(v)))) return null
  const rowSums = table.map(row => row.reduce((a, b) => a + b, 0))
  const colSums = table[0].map((_, j) => table.reduce((s, row) => s + row[j], 0))
  const total = rowSums.reduce((a, b) => a + b, 0)
  if (!(total > 0) || rowSums.some(v => v === 0) || colSums.some(v => v === 0)) return null
  let chi2 = 0
  let lowExpected = 0
  table.forEach((row, i) => row.forEach((obs, j) => {
    const exp = (rowSums[i] * colSums[j]) / total
    if (exp < 5) lowExpected++
    chi2 += (obs - exp) ** 2 / exp
  }))
  const df = (r - 1) * (c - 1)
  const pValue = chiSquareP(chi2, df)
  const cramerV = Math.sqrt(chi2 / (total * (Math.min(r, c) - 1)))
  return {
    chi2, df, pValue, cramerV, total,
    significant: pValue < 1 - Number(confidence) / 100,
    lowExpectedPct: (lowExpected / (r * c)) * 100,
  }
}

function pearsonR(x: number[], y: number[]): number | null {
  const n = x.length
  const mx = x.reduce((a, b) => a + b, 0) / n
  const my = y.reduce((a, b) => a + b, 0) / n
  let sxy = 0, sxx = 0, syy = 0
  for (let i = 0; i < n; i++) {
    sxy += (x[i] - mx) * (y[i] - my)
    sxx += (x[i] - mx) ** 2
    syy += (y[i] - my) ** 2
  }
  if (sxx === 0 || syy === 0) return null
  return sxy / Math.sqrt(sxx * syy)
}

// Rangos promedio (los empates reciben el promedio de sus posiciones).
function ranks(values: number[]): number[] {
  const order = values.map((v, i) => [v, i] as const).sort((a, b) => a[0] - b[0])
  const out = new Array<number>(values.length)
  for (let i = 0; i < order.length;) {
    let j = i
    while (j + 1 < order.length && order[j + 1][0] === order[i][0]) j++
    const rank = (i + j) / 2 + 1
    for (let k = i; k <= j; k++) out[order[k][1]] = rank
    i = j + 1
  }
  return out
}

// Correlación de Pearson o Spearman con su p bilateral (t con n − 2 gl).
export function correlation(x: number[], y: number[], method: 'pearson' | 'spearman') {
  const n = x.length
  if (n < 3 || y.length !== n) return null
  const r = method === 'pearson' ? pearsonR(x, y) : pearsonR(ranks(x), ranks(y))
  if (r == null) return null
  const df = n - 2
  const pValue = Math.abs(r) >= 1 ? 0 : tTestP(r * Math.sqrt(df / (1 - r * r)), df)
  return { r, pValue, n }
}

// Escala habitual en tesis para la fuerza de una correlación.
export function correlationStrength(r: number): 'veryWeak' | 'weak' | 'moderate' | 'strong' | 'veryStrong' {
  const a = Math.abs(r)
  if (a >= 0.8) return 'veryStrong'
  if (a >= 0.6) return 'strong'
  if (a >= 0.4) return 'moderate'
  if (a >= 0.2) return 'weak'
  return 'veryWeak'
}

// Fuerza de la asociación según la V de Cramér (Cohen, para tablas pequeñas).
export function cramerStrength(v: number): 'negligible' | 'weak' | 'moderate' | 'strong' {
  if (v >= 0.5) return 'strong'
  if (v >= 0.3) return 'moderate'
  if (v >= 0.1) return 'weak'
  return 'negligible'
}

// ---------------------------------------------------------------------------
// Herramientas para empresas: costo de rotación, churn/CLV, Van Westendorp y
// matriz importancia-satisfacción.
// ---------------------------------------------------------------------------

// Costo de cada salida = reclutamiento + capacitación + productividad perdida
// mientras la persona nueva aprende (sueldo × meses × 50 %). Estimación
// conservadora: no incluye el tiempo con el puesto vacío ni clientes perdidos.
export function turnoverCost(input: { salary: number; recruiting: number; training: number; rampMonths: number; exits: number; headcount: number | null }) {
  const { salary, recruiting, training, rampMonths, exits, headcount } = input
  if ([salary, recruiting, training, rampMonths, exits].some(v => !(v >= 0)) || !(salary > 0)) return null
  if (headcount != null && !(headcount > 0)) return null
  const ramp = salary * rampMonths * 0.5
  const perExit = recruiting + training + ramp
  return {
    perExit,
    parts: { recruiting, training, ramp },
    months: perExit / salary,
    annual: perExit * exits,
    turnoverRatePct: headcount != null ? (exits / headcount) * 100 : null,
  }
}

// Churn mensual y valor de vida del cliente con churn constante.
export function churnClv(start: number, lost: number, revenue: number, marginPct: number) {
  if (!(start > 0) || !(lost >= 0) || lost > start || !(revenue >= 0) || !(marginPct > 0 && marginPct <= 100)) return null
  const churn = lost / start
  const lifetime = churn > 0 ? 1 / churn : null
  const clvFor = (c: number) => (c > 0 ? (revenue * (marginPct / 100)) / c : null)
  const better = churn > 0.01 ? clvFor(churn - 0.01) : null
  return {
    churnPct: churn * 100,
    retentionPct: (1 - churn) * 100,
    lifetimeMonths: lifetime,
    clv: clvFor(churn),
    clvIfBetter: better,
    annualLost: (1 - (1 - churn) ** 12) * start,
  }
}

// Van Westendorp. rows: [demasiado barato, barato, caro, demasiado caro].
export function vanWestendorp(rows: number[][]) {
  const valid = rows.filter(r => r.length === 4 && r.every(v => v >= 0) && r[0] <= r[1] && r[1] <= r[2] && r[2] <= r[3])
  const n = valid.length
  if (n < 5) return { error: 'empty' as const, dropped: rows.length - n }
  const prices = [...new Set(valid.flat())].sort((a, b) => a - b)
  const share = (pred: (r: number[]) => boolean) => valid.filter(pred).length / n
  const curves = prices.map(p => ({
    price: p,
    tooCheap: share(r => r[0] >= p),
    cheap: share(r => r[1] >= p),
    expensive: share(r => r[2] <= p),
    tooExpensive: share(r => r[3] <= p),
  }))
  type Key = 'tooCheap' | 'cheap' | 'expensive' | 'tooExpensive'
  const value = (c: (typeof curves)[number], k: Key | 'notCheap' | 'notExpensive') =>
    k === 'notCheap' ? 1 - c.cheap : k === 'notExpensive' ? 1 - c.expensive : c[k]
  // Primer cruce entre dos curvas (con interpolación lineal entre precios).
  const cross = (a: Key | 'notCheap' | 'notExpensive', b: Key | 'notCheap' | 'notExpensive'): number | null => {
    for (let i = 0; i < curves.length; i++) {
      const d = value(curves[i], a) - value(curves[i], b)
      if (d === 0) return curves[i].price
      if (i > 0) {
        const prev = value(curves[i - 1], a) - value(curves[i - 1], b)
        if (prev !== 0 && Math.sign(prev) !== Math.sign(d)) {
          const t = prev / (prev - d)
          return curves[i - 1].price + t * (curves[i].price - curves[i - 1].price)
        }
      }
    }
    return null
  }
  return {
    n,
    dropped: rows.length - n,
    curves,
    pmc: cross('tooCheap', 'notCheap'),
    pme: cross('tooExpensive', 'notExpensive'),
    opp: cross('tooCheap', 'tooExpensive'),
    ipp: cross('cheap', 'expensive'),
  }
}

// Matriz importancia-satisfacción: cuadrantes cortados en los promedios.
export function importanceMatrix(items: { name: string; importance: number; satisfaction: number }[]) {
  if (items.length < 3) return null
  const meanI = items.reduce((s, i) => s + i.importance, 0) / items.length
  const meanS = items.reduce((s, i) => s + i.satisfaction, 0) / items.length
  const quadrant = (i: { importance: number; satisfaction: number }): 'fix' | 'keep' | 'low' | 'over' =>
    i.importance >= meanI ? (i.satisfaction < meanS ? 'fix' : 'keep') : (i.satisfaction < meanS ? 'low' : 'over')
  const rows = items.map(i => ({ ...i, gap: i.importance - i.satisfaction, quadrant: quadrant(i) }))
  return { meanImportance: meanI, meanSatisfaction: meanS, rows, byGap: [...rows].sort((a, b) => b.gap - a.gap) }
}

// ¿Cambió el NPS entre dos grupos (meses, sucursales)? Prueba z con la varianza
// del NPS: 100² × (pP + pD − (pP − pD)²) / n. También estima cuántas respuestas
// por grupo harían falta para detectar una diferencia de ese tamaño (95 %, 80 %).
export function npsCompare(a: [number, number, number], b: [number, number, number]) {
  const one = ([p, pa, d]: [number, number, number]) => {
    const n = p + pa + d
    if ([p, pa, d].some(v => !(v >= 0)) || n < 1) return null
    const pp = p / n, pd = d / n
    const v1 = 10000 * (pp + pd - (pp - pd) ** 2)
    return { n, nps: 100 * (pp - pd), v1, variance: v1 / n }
  }
  const A = one(a), B = one(b)
  if (!A || !B) return null
  const diff = B.nps - A.nps
  const se = Math.sqrt(A.variance + B.variance)
  const pValue = se === 0 ? (diff === 0 ? 1 : 0) : Math.min(1, 2 * (1 - normalCdf(Math.abs(diff / se))))
  const needed = diff !== 0 ? Math.ceil(((1.96 + 0.8416) ** 2 * (A.v1 + B.v1)) / diff ** 2) : null
  return { a: A, b: B, diff, low: diff - 1.96 * se, high: diff + 1.96 * se, pValue, significant: pValue < 0.05, needed }
}

// Lista de números pegada (una por línea o separados por comas, espacios o
// punto y coma). Acepta coma decimal cuando no se usa como separador.
export function parseNumberList(text: string): { values: number[]; skipped: number } {
  const t = text.trim()
  if (!t) return { values: [], skipped: 0 }
  // "4, 5, 3" o "4,5,3" en una línea: la coma separa. Con saltos de línea,
  // tabulaciones, punto y coma o espacios: la coma es decimal ("4,5").
  const commaSeparates = /,\s/.test(t) || !/[\s;]/.test(t)
  const parts = t.split(commaSeparates ? /[\s;,]+/ : /[\s;]+/)
  let skipped = 0
  const values: number[] = []
  for (const raw of parts) {
    if (!raw) continue
    const v = Number(commaSeparates ? raw : raw.replace(',', '.'))
    if (Number.isFinite(v)) values.push(v)
    else skipped++
  }
  return { values, skipped }
}

const mean = (x: number[]) => x.reduce((s, v) => s + v, 0) / x.length
const sampleSd = (x: number[]) => {
  const m = mean(x)
  return Math.sqrt(x.reduce((s, v) => s + (v - m) ** 2, 0) / (x.length - 1))
}

// t crítica bilateral (p = 1 − confianza) por bisección sobre tTestP.
export function tCritical(df: number, alpha = 0.05) {
  let lo = 0, hi = 1000
  for (let i = 0; i < 100; i++) {
    const mid = (lo + hi) / 2
    if (tTestP(mid, df) > alpha) lo = mid
    else hi = mid
  }
  return (lo + hi) / 2
}

// Dos promedios: prueba t de Welch (no supone varianzas iguales).
export function meanCompare(a: number[], b: number[]) {
  if (a.length < 2 || b.length < 2) return null
  const mA = mean(a), mB = mean(b), sA = sampleSd(a), sB = sampleSd(b)
  const vA = sA ** 2 / a.length, vB = sB ** 2 / b.length
  const diff = mB - mA
  const se = Math.sqrt(vA + vB)
  if (se === 0) return { a: { n: a.length, mean: mA, sd: sA }, b: { n: b.length, mean: mB, sd: sB }, diff, low: diff, high: diff, pValue: diff === 0 ? 1 : 0, significant: diff !== 0, df: a.length + b.length - 2 }
  const df = (vA + vB) ** 2 / (vA ** 2 / (a.length - 1) + vB ** 2 / (b.length - 1))
  const pValue = tTestP(diff / se, df)
  const tc = tCritical(df)
  return { a: { n: a.length, mean: mA, sd: sA }, b: { n: b.length, mean: mB, sd: sB }, diff, low: diff - tc * se, high: diff + tc * se, pValue, significant: pValue < 0.05, df }
}

// Ranking ponderado (promedio bayesiano): ajustada = (n × prom + m × C) / (n + m),
// con C el promedio general ponderado por respuestas y m por defecto la mediana
// de respuestas. "move" = puestos que sube (+) o baja (−) frente al ranking crudo.
export function weightedRanking(items: { name: string; average: number; n: number }[], weight: number | null) {
  const valid = items.filter(i => Number.isFinite(i.average) && i.n > 0)
  if (valid.length < 2) return null
  const total = valid.reduce((s, i) => s + i.n, 0)
  const overall = valid.reduce((s, i) => s + i.average * i.n, 0) / total
  const ns = valid.map(i => i.n).sort((x, y) => x - y)
  const mid = Math.floor(ns.length / 2)
  const median = ns.length % 2 ? ns[mid] : (ns[mid - 1] + ns[mid]) / 2
  const m = weight != null && weight >= 0 ? weight : median
  const indexed = valid.map((i, idx) => ({ ...i, idx }))
  const rawRank = new Map([...indexed].sort((x, y) => y.average - x.average).map((r, k) => [r.idx, k + 1]))
  const rows = indexed
    .map(i => ({ ...i, adjusted: (i.n * i.average + m * overall) / (i.n + m), few: i.n < m / 2 }))
    .sort((x, y) => y.adjusted - x.adjusted)
    .map((r, k) => ({ name: r.name, average: r.average, n: r.n, adjusted: r.adjusted, few: r.few, rank: k + 1, move: rawRank.get(r.idx)! - (k + 1) }))
  return { overall, weight: m, median, rows }
}
