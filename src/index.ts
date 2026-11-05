// surviq-stats: survey and thesis statistics, the same formulas behind the free
// calculators at https://surviqapp.com/en/tools. No dependencies.
export * from './stats.js'
export * from './thesis.js'
export { kr20, difficultyLevel, discriminationLevel, normalityTest } from './reliability.js'
export {
  mean, sd, median, normalQuantile, reliability, omegaTotal, shapiroWilk, lilliefors, normality, SHAPIRO_MAX_N,
  levelCuts, displayCuts, levelIndex,
} from './psychometrics.js'
export type { Reliability, ItemReliability, Normality, LevelMethod } from './psychometrics.js'
