import { SHAPIRO_MAX_N, lilliefors as ksThesis, normalQuantile, shapiroWilk as swThesis } from './psychometrics.js'

// KR-20 / KR-21 con análisis de ítems, y pruebas de normalidad (Shapiro-Wilk y
// Kolmogorov-Smirnov con corrección de Lilliefors). Todo corre en el navegador.

const mean = (v: number[]) => v.reduce((a, b) => a + b, 0) / v.length
// Varianza poblacional (÷ n): es la que usa la fórmula clásica del KR-20 con p·q.
const popVariance = (v: number[]) => {
  const m = mean(v)
  return v.reduce((s, x) => s + (x - m) ** 2, 0) / v.length
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
  return sxx === 0 || syy === 0 ? null : sxy / Math.sqrt(sxx * syy)
}

function kr20Of(rows: number[][], items: number[]): number | null {
  const k = items.length
  if (k < 2) return null
  const totals = rows.map(r => items.reduce((s, j) => s + r[j], 0))
  const variance = popVariance(totals)
  if (variance === 0) return null
  const sumPq = items.reduce((s, j) => {
    const p = mean(rows.map(r => r[j]))
    return s + p * (1 - p)
  }, 0)
  return (k / (k - 1)) * (1 - sumPq / variance)
}

// Filas = personas, columnas = ítems, cada celda 0 (incorrecto / no) o 1 (correcto / sí).
export function kr20(rows: number[][]) {
  if (rows.length < 2 || (rows[0]?.length ?? 0) < 2) return { error: 'size' as const }
  if (rows.some(r => r.some(v => v !== 0 && v !== 1))) return { error: 'binary' as const }
  const k = rows[0].length
  const items = rows[0].map((_, j) => j)
  const value = kr20Of(rows, items)
  if (value == null) return { error: 'variance' as const }
  const totals = rows.map(r => r.reduce((a, b) => a + b, 0))
  const m = mean(totals)
  const variance = popVariance(totals)
  const kr21 = (k / (k - 1)) * (1 - (m * (k - m)) / (k * variance))
  const itemStats = items.map(j => {
    const scores = rows.map(r => r[j])
    const rest = rows.map((r, i) => totals[i] - r[j])
    return {
      difficulty: mean(scores),
      // Correlación ítem-resto (biserial puntual corregida): cuánto separa el
      // ítem a quienes sacan más puntaje en el resto de la prueba.
      discrimination: pearson(scores, rest),
      ifDeleted: k > 2 ? kr20Of(rows, items.filter(i => i !== j)) : null,
    }
  })
  return { kr20: value, kr21, items: k, respondents: rows.length, mean: m, itemStats }
}

export function difficultyLevel(p: number): 'veryEasy' | 'easy' | 'medium' | 'hard' | 'veryHard' {
  if (p > 0.9) return 'veryEasy'
  if (p > 0.7) return 'easy'
  if (p >= 0.3) return 'medium'
  if (p >= 0.1) return 'hard'
  return 'veryHard'
}

// Ebel: 0,40 o más muy bueno, 0,30 a 0,39 bueno, 0,20 a 0,29 revisar, menos de 0,20 descartar o rehacer.
export function discriminationLevel(d: number): 'veryGood' | 'good' | 'review' | 'poor' {
  if (d >= 0.4) return 'veryGood'
  if (d >= 0.3) return 'good'
  if (d >= 0.2) return 'review'
  return 'poor'
}

// Shapiro-Wilk y Kolmogorov-Smirnov (Lilliefors): se reusan los del Modo Tesis para
// que la herramienta gratis y la app den exactamente el mismo resultado.
export { normalQuantile }

export function shapiroWilk(values: number[]): { w: number; p: number } | null {
  const r = swThesis(values)
  return r ? { w: r.w, p: r.pValue } : null
}

export function lilliefors(values: number[]): { d: number; p: number } | null {
  const r = ksThesis(values)
  return r ? { d: r.d, p: r.pValue } : null
}

// Resumen para decidir: con 50 datos o menos se reporta Shapiro-Wilk; con más,
// Kolmogorov-Smirnov (la regla habitual en tesis). Se muestran las dos.
export function normalityTest(values: number[], alpha = 0.05) {
  if (values.length < 3) return null
  const sw = shapiroWilk(values)
  const ks = lilliefors(values)
  if (!sw && !ks) return null
  const recommended: 'sw' | 'ks' = values.length <= SHAPIRO_MAX_N || !ks ? 'sw' : 'ks'
  const p = recommended === 'sw' ? sw?.p : ks?.p
  const avg = mean(values)
  const sd = Math.sqrt(values.reduce((s, v) => s + (v - avg) ** 2, 0) / (values.length - 1))
  return { n: values.length, mean: avg, sd, sw, ks, recommended, normal: p != null ? p > alpha : null }
}
