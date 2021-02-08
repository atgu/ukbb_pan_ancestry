
export enum PopulationCode {
  AFR = "AFR",
  AMR = "AMR",
  CSA = "CSA",
  EAS = "EAS",
  EUR = "EUR",
  MID = "MID",
}
export const populationColorMapping = new Map<PopulationCode, string>([
  [PopulationCode.AFR, "#941494"],
  [PopulationCode.AMR, "#ed1e24"],
  [PopulationCode.CSA, "#FF9912"],
  [PopulationCode.EAS, "#108c44"],
  [PopulationCode.EUR, "#6aa5cd"],
  [PopulationCode.MID, "#EEA9B8"],
])

// These are the only populations that actually occur in the data set:
export const commonPopulations = [...populationColorMapping.keys()]

export type PerPopulationMetrics = Map<PopulationCode, number>
