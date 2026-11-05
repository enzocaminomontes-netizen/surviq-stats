import { describe, expect, it } from 'vitest'
import { churnClv, importanceMatrix, turnoverCost, vanWestendorp, aikenV, alphaLevel, cesScore, chiSquareP, chiSquareTest, correlation, correlationStrength, likertLevels, parseLabeledTable, tTestP, cronbachAlpha, csatScore, invitationsNeeded, marginOfError, normalCdf, npsScore, parseScoreTable, responseRates, sampleSize, twoProportionTest } from '../src/stats.js'

describe('herramientas gratis', () => {
  it('tamaño de muestra', () => {
    expect(sampleSize(null, '95', 5)).toBe(385)
    expect(sampleSize(1000, '95', 5)).toBe(278)
    expect(sampleSize(50, '99', 1)).toBe(50)
    expect(sampleSize(null, '95', 0)).toBeNull()
    expect(invitationsNeeded(385, 20)).toBe(1925)
  })
  it('margen de error', () => {
    expect(marginOfError(385, null, '95')!).toBeCloseTo(4.99, 2)
    expect(marginOfError(278, 1000, '95')!).toBeCloseTo(4.99, 1)
    expect(marginOfError(10, 5, '95')).toBeNull()
  })
  it('NPS y CSAT', () => {
    expect(npsScore(50, 30, 20)!.score).toBe(30)
    expect(npsScore(0, 0, 0)).toBeNull()
    const c = csatScore([1, 1, 2, 3, 3])!
    expect(c.csatPct).toBe(60)
    expect(c.average).toBeCloseTo(3.6)
  })
  it('tasa de respuesta', () => {
    const r = responseRates(200, 80, 50)!
    expect(r.responseRatePct).toBe(25)
    expect(r.completionRatePct).toBe(62.5)
    expect(responseRates(10, 3, 5)).toBeNull()
  })
  it('dos proporciones', () => {
    expect(normalCdf(1.96)).toBeCloseTo(0.975, 3)
    const r = twoProportionTest(60, 100, 45, 100, '95')!
    expect(r.z).toBeCloseTo(2.12, 2)
    expect(r.pValue).toBeCloseTo(0.034, 2)
    expect(r.significant).toBe(true)
    expect(twoProportionTest(52, 100, 50, 100, '95')!.significant).toBe(false)
  })
  it('alfa de Cronbach', () => {
    const parsed = parseScoreTable('P1\tP2\tP3\n4\t5\t4\n3\t3\t2\n5\t5\t5\n2\t3\t2\n4\t4\t5')
    expect('rows' in parsed && parsed.skippedHeader).toBe(true)
    if (!('rows' in parsed)) throw new Error()
    const a = cronbachAlpha(parsed.rows)!
    expect(a.alpha).toBeCloseTo(0.9344, 3)
    expect(a.ifDeleted).toHaveLength(3)
    expect(alphaLevel(a.alpha)).toBe('excellent')
    expect(parseScoreTable('1;2,5\n3;4')).toEqual({ rows: [[1, 2.5], [3, 4]], skippedHeader: false })
    expect(parseScoreTable('1,2,3\n4,5,6')).toEqual({ rows: [[1, 2, 3], [4, 5, 6]], skippedHeader: false })
    expect(parseScoreTable('1 2\n3')).toEqual({ error: 'columns' })
  })
  it('CES', () => {
    const r = cesScore([1, 1, 2, 2, 4, 6, 4])!
    expect(r.total).toBe(20)
    expect(r.easyPct).toBe(70)
    expect(r.hardPct).toBe(20)
    expect(r.average).toBeCloseTo(5.05, 5)
    expect(cesScore([0, 0, 1, 2, 2])!.easyPct).toBe(80)
    expect(cesScore([0, 0, 0])).toBeNull()
    expect(cesScore([0, 0, 0, 0, 0])).toBeNull()
  })
  it('V de Aiken', () => {
    const r = aikenV([4, 4, 3, 4, 4], 1, 4)!
    expect(r.v).toBeCloseTo(0.9333, 4)
    expect(r.lower).toBeCloseTo(0.701829, 4)
    expect(r.upper).toBeCloseTo(0.988133, 4)
    expect(aikenV([5, 4], 1, 4)).toBeNull()
  })
  it('baremos', () => {
    expect(likertLevels(10, 1, 5, 3)!.map(l => [l.from, l.to])).toEqual([[10, 23], [24, 36], [37, 50]])
    const avg = likertLevels(10, 1, 5, 3)!
    expect(avg[0].avgTo).toBeCloseTo(2.3333, 4)
    expect(likertLevels(1, 1, 2, 5)).toBeNull()
  })
  it('chi-cuadrado', () => {
    const r = chiSquareTest([[20, 15, 5], [10, 25, 25]])!
    expect(r.chi2).toBeCloseTo(15.7986, 3)
    expect(r.df).toBe(2)
    expect(r.pValue).toBeCloseTo(0.000371, 5)
    expect(chiSquareP(7.81, 3)).toBeCloseTo(0.050106, 5)
    expect(chiSquareTest([[1, 2]])).toBeNull()
  })
  it('correlación', () => {
    const x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    const y = [2, 1, 4, 3, 7, 8, 6, 9, 10, 8]
    const p = correlation(x, y, 'pearson')!
    expect(p.r).toBeCloseTo(0.893994, 5)
    expect(p.pValue).toBeCloseTo(0.000485, 5)
    const s = correlation(x, y, 'spearman')!
    expect(s.r).toBeCloseTo(0.881463, 5)
    expect(s.pValue).toBeCloseTo(0.000747, 5)
    expect(tTestP(2.5, 10)).toBeCloseTo(0.031447, 5)
    expect(correlationStrength(-0.45)).toBe('moderate')
  })
  it('tabla con etiquetas', () => {
    const t = parseLabeledTable('\tSí\tNo\nHombres\t20\t15\nMujeres\t10\t25')
    expect(t).toEqual({ rows: [[20, 15], [10, 25]], rowLabels: ['Hombres', 'Mujeres'], colLabels: ['Sí', 'No'] })
    const u = parseLabeledTable('Ítem\tJ1\tJ2\n1\t4\t3')
    expect(u).toEqual({ rows: [[4, 3]], rowLabels: ['1'], colLabels: ['J1', 'J2'] })
  })
  it('costo de rotación', () => {
    const r = turnoverCost({ salary: 2000, recruiting: 800, training: 500, rampMonths: 3, exits: 10, headcount: 50 })!
    expect(r.perExit).toBe(4300)
    expect(r.annual).toBe(43000)
    expect(r.turnoverRatePct).toBe(20)
    expect(r.months).toBeCloseTo(2.15, 5)
    expect(turnoverCost({ salary: 0, recruiting: 0, training: 0, rampMonths: 0, exits: 0, headcount: null })).toBeNull()
  })
  it('churn y CLV', () => {
    const r = churnClv(200, 10, 50, 60)!
    expect(r.churnPct).toBe(5)
    expect(r.lifetimeMonths).toBe(20)
    expect(r.clv).toBeCloseTo(600, 6)
    expect(r.clvIfBetter).toBeCloseTo(750, 6)
    expect(churnClv(10, 11, 5, 50)).toBeNull()
  })
  it('Van Westendorp', () => {
    const rows = [[5, 10, 20, 30], [6, 12, 22, 35], [4, 8, 18, 25], [5, 10, 25, 40], [7, 14, 24, 32], [3, 6, 15, 20], [20, 10, 5, 1]]
    const r = vanWestendorp(rows)
    if ('error' in r) throw new Error('no')
    expect(r.n).toBe(6)
    expect(r.dropped).toBe(1)
    expect(r.opp).not.toBeNull()
    expect(r.pmc! <= r.pme!).toBe(true)
    expect(r.opp! > r.pmc! && r.opp! < r.pme!).toBe(true)
    expect('error' in vanWestendorp([[1, 2, 3, 4]])).toBe(true)
  })
  it('matriz importancia-satisfacción', () => {
    const r = importanceMatrix([
      { name: 'A', importance: 4.8, satisfaction: 3.1 },
      { name: 'B', importance: 4.5, satisfaction: 4.6 },
      { name: 'C', importance: 3.0, satisfaction: 2.9 },
      { name: 'D', importance: 2.8, satisfaction: 4.7 },
    ])!
    expect(r.rows.map(x => x.quadrant)).toEqual(['fix', 'keep', 'low', 'over'])
    expect(r.byGap[0].name).toBe('A')
  })
})

import { npsCompare, meanCompare, parseNumberList, weightedRanking } from '../src/stats.js'

describe('comparison tools', () => {
  it('npsCompare', () => {
    const r = npsCompare([40, 30, 30], [50, 30, 20])!
    expect(r.a.nps).toBeCloseTo(10)
    expect(r.b.nps).toBeCloseTo(30)
    expect(r.diff).toBeCloseTo(20)
    // var A = 10000*(0.7-0.01)/100 = 69 ; var B = 10000*(0.7-0.09)/100 = 61 ; se = sqrt(130)
    expect(r.high - r.diff).toBeCloseTo(1.96 * Math.sqrt(130), 3)
    expect(r.pValue).toBeCloseTo(0.0794, 3)
    expect(r.significant).toBe(false)
    expect(npsCompare([0, 0, 0], [1, 0, 0])).toBeNull()
  })
  it('meanCompare matches scipy welch', () => {
    const r = meanCompare([4, 5, 3, 4, 5, 4, 3, 5], [3, 2, 4, 3, 3, 2, 4, 3, 2])!
    expect(r.diff).toBeCloseTo(2.888889 - 4.125, 5)
    expect(r.pValue).toBeCloseTo(0.006988, 4)
    expect(r.df).toBeCloseTo(14.4737, 3)
  })
  it('parseNumberList', () => {
    expect(parseNumberList('4\n5\n3,5').values).toEqual([4, 5, 3.5])
    expect(parseNumberList('4,5,3').values).toEqual([4, 5, 3])
    expect(parseNumberList('4, 5, x').values).toEqual([4, 5])
    expect(parseNumberList('4, 5, x').skipped).toBe(1)
  })
  it('weightedRanking', () => {
    const r = weightedRanking([{ name: 'A', average: 5, n: 3 }, { name: 'B', average: 4.7, n: 200 }, { name: 'C', average: 4.2, n: 80 }], null)!
    expect(r.weight).toBe(80)
    expect(r.rows[0].name).toBe('B')
    expect(r.rows.find(x => x.name === 'A')!.move).toBe(-1)
    expect(r.rows.find(x => x.name === 'A')!.few).toBe(true)
  })
})
