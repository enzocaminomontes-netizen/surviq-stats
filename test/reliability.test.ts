import { describe, expect, it } from 'vitest'
import { cronbachAlpha } from '../src/stats.js'
import { kr20, lilliefors, normalityTest, normalQuantile, shapiroWilk } from '../src/reliability.js'

const DATA = [
  [1, 1, 1, 0, 1], [1, 0, 1, 0, 0], [1, 1, 1, 1, 1], [0, 0, 1, 0, 0],
  [1, 1, 0, 1, 1], [1, 0, 1, 0, 1], [0, 1, 0, 0, 0], [1, 1, 1, 1, 0],
]

describe('kr20', () => {
  it('equals Cronbach alpha on 0/1 data', () => {
    const r = kr20(DATA)
    if ('error' in r) throw new Error(r.error)
    expect(r.kr20).toBeCloseTo(cronbachAlpha(DATA)!.alpha, 10)
    expect(r.itemStats[0].difficulty).toBeCloseTo(0.75)
    expect(r.kr21).toBeLessThanOrEqual(r.kr20 + 1e-9)
  })
  it('rejects non-binary data and zero variance', () => {
    expect(kr20([[1, 2], [0, 1]])).toEqual({ error: 'binary' })
    expect(kr20([[1, 1], [1, 1]])).toEqual({ error: 'variance' })
  })
})

describe('normality', () => {
  it('normalQuantile inverts the normal', () => {
    expect(normalQuantile(0.975)).toBeCloseTo(1.959964, 5)
    expect(normalQuantile(0.5)).toBe(0)
  })
  it('matches SciPy shapiro', () => {
    // scipy.stats.shapiro([2.1, 3.4, 1.9, 5.6, 4.4, 3.0, 2.8, 3.9, 4.1, 3.3])
    const r = shapiroWilk([2.1, 3.4, 1.9, 5.6, 4.4, 3.0, 2.8, 3.9, 4.1, 3.3])!
    expect(r.w).toBeCloseTo(0.97139, 4)
    expect(r.p).toBeCloseTo(0.90343, 3)
    const skewed = shapiroWilk([1, 1, 1, 2, 2, 3, 4, 8, 15, 40, 1, 2, 1, 1, 3])!
    expect(skewed.p).toBeLessThan(0.001)
  })
  it('lilliefors flags skewed data', () => {
    expect(lilliefors([1, 1, 1, 2, 2, 3, 4, 8, 15, 40, 1, 2, 1, 1, 3])!.p).toBeLessThan(0.01)
  })
  it('recommends SW up to 50 and KS above', () => {
    const small = Array.from({ length: 30 }, (_, i) => normalQuantile((i + 0.5) / 30))
    const big = Array.from({ length: 80 }, (_, i) => normalQuantile((i + 0.5) / 80))
    expect(normalityTest(small)!.recommended).toBe('sw')
    expect(normalityTest(small)!.normal).toBe(true)
    expect(normalityTest(big)!.recommended).toBe('ks')
  })
})
