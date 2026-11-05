import { describe, expect, it } from 'vitest'
import { contentValidity, cvcLevel, fleissKappa, independentTTest, interRaterAgreement, kappaLevel, kruskalWallis, mannWhitney, mcdonaldOmega, oneWayAnova, pairedTTest, parseCategoryTable, wilcoxonSignedRank, kendallW, binomialJudges, binomialUpperP } from '../src/thesis.js'

// Valores de referencia calculados con SciPy 1.18 (t, Levene, Mann-Whitney) y con
// un ajuste minres independiente (omega).
const OMEGA_DATA = [[3, 3, 2, 3, 4, 2], [4, 3, 3, 5, 4, 2], [3, 3, 3, 3, 3, 3], [3, 2, 3, 2, 3, 1], [2, 2, 3, 3, 1, 2], [3, 1, 2, 2, 4, 3], [3, 3, 3, 4, 3, 3], [4, 4, 4, 3, 4, 4], [3, 3, 3, 3, 2, 4], [2, 3, 1, 3, 1, 1], [3, 3, 3, 5, 2, 3], [3, 4, 3, 3, 4, 4], [2, 3, 3, 2, 3, 2], [3, 2, 2, 2, 2, 1], [2, 3, 1, 4, 1, 4], [3, 4, 4, 2, 5, 5], [2, 2, 2, 1, 3, 2], [3, 2, 2, 2, 4, 3], [2, 2, 1, 1, 1, 2], [2, 2, 1, 1, 4, 2], [1, 2, 3, 1, 2, 1], [2, 3, 3, 3, 2, 3], [2, 2, 1, 1, 4, 2], [3, 3, 3, 3, 4, 3], [3, 3, 4, 4, 3, 3], [2, 4, 4, 3, 3, 4], [1, 2, 1, 2, 1, 2], [3, 3, 4, 4, 2, 1], [3, 2, 3, 4, 1, 1], [3, 3, 3, 3, 2, 2], [2, 1, 1, 2, 2, 2], [2, 2, 2, 2, 3, 2], [2, 3, 4, 1, 3, 2], [2, 1, 2, 3, 2, 2], [4, 5, 4, 2, 3, 2], [1, 2, 4, 2, 1, 2], [4, 3, 3, 3, 3, 4], [4, 4, 3, 4, 3, 5], [2, 2, 3, 2, 4, 4], [3, 4, 3, 1, 3, 3]]

describe('omega de McDonald', () => {
  it('coincide con un ajuste minres de un factor', () => {
    const r = mcdonaldOmega(OMEGA_DATA)
    if ('error' in r) throw new Error(r.error)
    expect(r.omega).toBeCloseTo(0.76486, 4)
    expect(r.items.map(i => i.loading)).toEqual([0.744221, 0.782579, 0.585676, 0.402504, 0.430485, 0.575772].map(x => expect.closeTo(x, 3)))
    expect(r.items.every(i => i.ifDeleted != null)).toBe(true)
  })
  it('pide al menos 3 ítems y varianza', () => {
    expect(mcdonaldOmega([[1, 2], [2, 3], [3, 4]])).toEqual({ error: 'size' })
    expect(mcdonaldOmega([[1, 2, 3], [1, 3, 4], [1, 4, 2]])).toEqual({ error: 'variance' })
  })
})

describe('CVC de Hernández-Nieto', () => {
  it('calcula CVCi, Pei y el total', () => {
    const r = contentValidity([[5, 5, 4], [4, 4, 4], [3, 5, 4]], 5)
    if ('error' in r) throw new Error(r.error)
    expect(r.pe).toBeCloseTo(1 / 27, 10)
    expect(r.items[0].cvc).toBeCloseTo(14 / 15, 10)
    expect(r.items[0].corrected).toBeCloseTo(14 / 15 - 1 / 27, 10)
    expect(r.total).toBeCloseTo((14 / 15 + 12 / 15 + 12 / 15) / 3 - 1 / 27, 10)
    expect(cvcLevel(0.95)).toBe('excellent')
    expect(cvcLevel(0.75)).toBe('acceptable')
  })
  it('rechaza notas fuera de la escala', () => {
    expect(contentValidity([[6, 5]], 5)).toEqual({ error: 'range' })
  })
})

describe('kappa', () => {
  it('Cohen: ejemplo clásico de 50 casos da 0,40', () => {
    const rows = [
      ...Array(20).fill(['Sí', 'Sí']), ...Array(5).fill(['Sí', 'No']), ...Array(10).fill(['No', 'Sí']), ...Array(15).fill(['No', 'No']),
    ]
    const r = interRaterAgreement(rows)
    if ('error' in r || r.method !== 'cohen') throw new Error('cohen')
    expect(r.kappa).toBeCloseTo(0.4, 10)
    expect(r.po).toBeCloseTo(0.7, 10)
    expect(r.pValue).toBeLessThan(0.01)
    expect(kappaLevel(r.kappa)).toBe('fair')
  })
  it('Fleiss: ejemplo de Fleiss (1971) da 0,210', () => {
    const counts = [[0, 0, 0, 0, 14], [0, 2, 6, 4, 2], [0, 0, 3, 5, 6], [0, 3, 9, 2, 0], [2, 2, 8, 1, 1], [7, 7, 0, 0, 0], [3, 2, 6, 3, 0], [2, 5, 3, 2, 2], [6, 5, 2, 1, 0], [0, 2, 2, 3, 7]]
    const rows = counts.map(c => c.flatMap((n, j) => Array(n).fill(String(j + 1))))
    const r = fleissKappa(rows)!
    expect(r.kappa).toBeCloseTo(0.20993, 4)
    expect(r.observed).toBeCloseTo(0.37802, 4)
    expect(r.expected).toBeCloseTo(0.21276, 4)
  })
  it('lee categorías con texto y salta el encabezado', () => {
    const p = parseCategoryTable('Juez 1\tJuez 2\nAlto\tAlto\nBajo\tAlto\nBajo\tBajo')
    expect(p).toEqual({ rows: [['Alto', 'Alto'], ['Bajo', 'Alto'], ['Bajo', 'Bajo']], skippedHeader: true })
  })
})

describe('t de Student', () => {
  const a = [12, 15, 11, 18, 14, 16, 13, 17, 15, 14, 19, 12]
  const b = [10, 13, 9, 12, 11, 14, 10, 12, 13, 9, 11]
  it('muestras independientes, Welch y Levene como SciPy', () => {
    const r = independentTTest(a, b)!
    expect(r.student.t).toBeCloseTo(3.786085, 5)
    expect(r.student.pValue).toBeCloseTo(0.00108225, 6)
    expect(r.welch.t).toBeCloseTo(3.851778, 5)
    expect(r.welch.df).toBeCloseTo(19.342329, 4)
    expect(r.welch.pValue).toBeCloseTo(0.00104609, 6)
    expect(r.levene!.f).toBeCloseTo(1.630673, 5)
    expect(r.levene!.pValue).toBeCloseTo(0.21554966, 5)
    expect(r.recommended).toBe('student')
  })
  it('muestras relacionadas', () => {
    const r = pairedTTest([20, 22, 19, 24, 21, 23, 18, 25], [18, 21, 19, 20, 19, 22, 17, 21])!
    expect(r.t).toBeCloseTo(3.638034, 5)
    expect(r.pValue).toBeCloseTo(0.00830874, 6)
  })
})

describe('U de Mann-Whitney', () => {
  it('con empates: aproximación normal como SciPy', () => {
    const r = mannWhitney([12, 15, 11, 18, 14, 16, 13, 17, 15, 14, 19, 12], [10, 13, 9, 12, 11, 14, 10, 12, 13, 9, 11])!
    expect(r.u).toBe(17)
    expect(r.pValue).toBeCloseTo(0.00240902, 6)
    expect(r.exact).toBeNull()
  })
  it('sin empates: p exacta', () => {
    const r = mannWhitney([3.1, 4.5, 2.2, 5.8, 4.1, 3.9, 6.2], [1.5, 2.7, 3.3, 0.9, 2.1, 1.8])!
    expect(r.u).toBe(3)
    expect(r.exact).toBeCloseTo(0.00815851, 7)
    expect(r.pValue).toBeCloseTo(0.01012799, 6)
  })
})

describe('Wilcoxon, Kruskal-Wallis y ANOVA (referencias de SciPy y pingouin)', () => {
  it('Wilcoxon con empates: asintótica sin corrección', () => {
    const r = wilcoxonSignedRank([18, 21, 19, 20, 19, 22, 17, 21, 16, 20, 23, 18], [20, 22, 19, 24, 21, 23, 18, 25, 19, 21, 22, 22])!
    expect(r.t).toBe(3)
    expect(r.zeros).toBe(1)
    expect(r.pValue).toBeCloseTo(0.00691532, 6)
    expect(r.exact).toBeNull()
  })
  it('Wilcoxon sin empates: p exacta', () => {
    const r = wilcoxonSignedRank([3.1, 4.5, 2.2, 5.8, 4.1, 3.9, 6.2, 5.0], [3.9, 5.6, 2.0, 7.4, 5.0, 5.6, 5.9, 6.8])!
    expect(r.t).toBe(3)
    expect(r.exact).toBeCloseTo(0.0390625, 8)
    expect(r.pValue).toBeCloseTo(0.0356919, 6)
  })
  const g = [[12, 15, 11, 18, 14, 16, 13], [10, 13, 9, 12, 11, 14, 10, 12], [15, 17, 16, 19, 14, 18, 20]]
  it('Kruskal-Wallis con corrección por empates', () => {
    const r = kruskalWallis(g)!
    expect(r.h).toBeCloseTo(12.667615, 5)
    expect(r.pValue).toBeCloseTo(0.00177526, 6)
    expect(r.pairs).toHaveLength(3)
  })
  it('ANOVA clásico, Levene y Welch', () => {
    const r = oneWayAnova(g)!
    expect(r.f).toBeCloseTo(13.57529, 4)
    expect(r.pValue).toBeCloseTo(0.000218, 6)
    expect(r.eta2).toBeCloseTo(0.588304, 5)
    expect(r.levene!.f).toBeCloseTo(0.4228657, 5)
    expect(r.levene!.pValue).toBeCloseTo(0.661183, 5)
    expect(r.welch!.f).toBeCloseTo(14.958199, 4)
    expect(r.welch!.df2).toBeCloseTo(11.874602, 4)
    expect(r.welch!.pValue).toBeCloseTo(0.000569, 5)
    expect(r.recommended).toBe('classic')
  })
})

describe('kendallW', () => {
  it('matches scipy with ties (6 ítems, 3 jueces)', () => {
    const r = kendallW([[4, 5, 4], [3, 3, 2], [5, 4, 5], [2, 2, 3], [4, 4, 4], [1, 2, 1]])!
    expect(r.w).toBeCloseTo(0.900990, 5)
    expect(r.chi2).toBeCloseTo(13.514851, 5)
    expect(r.df).toBe(5)
    expect(r.pValue).toBeCloseTo(0.019003, 5)
    expect(r.ties).toBe(true)
  })
  it('gives W = 1 when all judges rank the same', () => {
    expect(kendallW([[1, 1], [2, 2], [3, 3], [4, 4]])!.w).toBeCloseTo(1, 10)
  })
  it('needs 3 items and 2 judges', () => {
    expect(kendallW([[1, 2], [2, 1]])).toBeNull()
  })
})

describe('binomialJudges', () => {
  it('gives the exact one-tailed p per item', () => {
    const r = binomialJudges([[1, 1, 1, 1, 1], [1, 1, 1, 1, 0], [1, 1, 1, 1, 1], [1, 1, 1, 1, 1], [1, 1, 1, 1, 1], [1, 1, 1, 1, 1]])
    if ('error' in r) throw new Error()
    expect(r.items[0].p).toBeCloseTo(0.03125, 10)
    expect(r.items[1].p).toBeCloseTo(0.1875, 10)
    expect(r.items[1].pTwo).toBeCloseTo(0.375, 10)
    expect(r.totalAgree).toBe(29)
    expect(r.totalPct).toBeCloseTo(29 / 30, 10)
  })
  it('matches scipy for the global test', () => {
    expect(binomialUpperP(28, 30)).toBeCloseTo(4.33996e-7, 10)
  })
  it('rejects values other than 0 and 1', () => {
    expect(binomialJudges([[1, 2], [0, 1]])).toEqual({ error: 'binary' })
  })
})
