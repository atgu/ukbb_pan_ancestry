import min from "lodash/min"
import max from "lodash/max"
import { PerPopulationExtremums } from "./PhenotypeFilters";
import { PerPopulationMetricsVisibility, PerPopulationRangeFilter } from "./phenotypesReducer"
import { commonPopulations, PopulationCode } from "./populations";
import { Datum } from "./types"
import identity from "lodash/identity"


export const maxSaigeHeritabilityValue = 1
const clampedLambdaGcMin = 0.5
const clampedLambdaGcMax = 1.5

// Used to collect all values for each population of a paticular metric for aggregation:
type PerPopulationAccumulators = {
  [K in PopulationCode]: number[]
}
const getInitialAccumulators = (): PerPopulationAccumulators => ({
  [PopulationCode.AFR]: [],
  [PopulationCode.AMR]: [],
  [PopulationCode.CSA]: [],
  [PopulationCode.EAS]: [],
  [PopulationCode.EUR]: [],
  [PopulationCode.MID]: [],
})
const getBarChartExtremums = (inputs: {
  data: Datum[]
  metric_prefix: string
  perPopulationMetricsVisibilities: PerPopulationMetricsVisibility,
}) => {
  const {data, metric_prefix, perPopulationMetricsVisibilities } = inputs;
  const nonEurPopulation = commonPopulations.filter(pop => pop !== PopulationCode.EUR)
  const allNonEurValues = []
  const allEurValues = []
  const perPopulationValues = getInitialAccumulators()
  for (const row of data) {
    for (const population of nonEurPopulation) {
      const nonEurValue = row[`${metric_prefix}_${population}`]
      if (typeof nonEurValue === "number") {
        if (perPopulationMetricsVisibilities[population] === true) {
          allNonEurValues.push(nonEurValue)
        }
        perPopulationValues[population].push(nonEurValue)
      }
    }
    const eurValue = row[`${metric_prefix}_${PopulationCode.EUR}`]
    if (typeof eurValue === "number") {
      if (perPopulationMetricsVisibilities[PopulationCode.EUR] === true) {
        allEurValues.push(eurValue)
      }
      perPopulationValues[PopulationCode.EUR].push(eurValue)
    }
  }

  const numVisiblePopulations = Object.values(perPopulationMetricsVisibilities).filter(identity).length
  const isOnlyEurVisible = (numVisiblePopulations === 1) && perPopulationMetricsVisibilities[PopulationCode.EUR] === true
  let threshold: number
  if (isOnlyEurVisible === true) {
    threshold = max(allEurValues)!
  } else {
    threshold = max(allNonEurValues)!
  }
  const perPopulationExtremums = Object.fromEntries(
    Object.entries(perPopulationValues).map(([population, values]) => ([
      population , {min: min(values)!, max: max(values)!}
    ]))
  ) as PerPopulationExtremums
  return {
    threshold, perPopulationExtremums,
  }
}

interface Inputs {
  data: Datum[]
  perPopulationMetricsVisibilities: PerPopulationMetricsVisibility,
  nCasesFilters: PerPopulationRangeFilter,
  nControlsFilters: PerPopulationRangeFilter,
  saigeHeritabilityFilters: PerPopulationRangeFilter,
  lambdaGcFilters: PerPopulationRangeFilter
}
export const determineExtremums =  (inputs: Inputs) => {
  const {
    data,
    perPopulationMetricsVisibilities,
    nCasesFilters,
    nControlsFilters,
    saigeHeritabilityFilters,
    lambdaGcFilters,
  } = inputs;
  const {
    threshold: globalNCasesUpperThreshold,
    perPopulationExtremums: perPopulationNCasesExtremums,
  } = getBarChartExtremums({
    data,
    metric_prefix: "n_cases",
    perPopulationMetricsVisibilities,
  })
  const {
    threshold: globalNControlsUpperThreshold,
    perPopulationExtremums: perPopulationNControlsExtremums,
  } = getBarChartExtremums({
    data,
    metric_prefix: "n_controls",
    perPopulationMetricsVisibilities,
    perPopulationRangeFilters: nControlsFilters,
  })
  const {
    perPopulationExtremums: perPopulationSaigeHeritabilityExtremums,
  } = getBarChartExtremums({
    data,
    metric_prefix: "saige_heritability",
    perPopulationMetricsVisibilities,
    perPopulationRangeFilters: saigeHeritabilityFilters,
  })


  const allLambdaGcValues = []
  const perPopulationLambdaGcValues = getInitialAccumulators()
  for (const row of data) {
    for (const population of commonPopulations) {
      const lambdaGcValue = row[`lambda_gc_${population}`]
      if (typeof lambdaGcValue === "number") {
        perPopulationLambdaGcValues[population].push(lambdaGcValue)
      }
      if (perPopulationMetricsVisibilities[population] === true) {
        if (typeof lambdaGcValue === "number") {
          allLambdaGcValues.push(lambdaGcValue)
        }
      }
    }
  }

  const getMetricValueForPopulationExtractor =
    (metricPrefix: string, population: PopulationCode) => (row: Datum) => row[`${metricPrefix}${population}`]

  let filteredByPopulationRangeFilters: Datum[] = data
  for (const [population, filterValue] of Object.entries(nCasesFilters)) {
    if (filterValue !== undefined) {
      const extractor = getMetricValueForPopulationExtractor("n_cases_", population as PopulationCode)
      filteredByPopulationRangeFilters = filteredByPopulationRangeFilters.filter(elem => {
        const extractedValue = extractor(elem)
        return (
          typeof extractedValue === "number" &&
          extractedValue >= filterValue.min &&
          extractedValue <= filterValue.max
        )
      })
    }
  }
  for (const [population, filterValue] of Object.entries(nControlsFilters)) {
    if (filterValue !== undefined) {
      const extractor = getMetricValueForPopulationExtractor("n_controls_", population as PopulationCode)
      filteredByPopulationRangeFilters = filteredByPopulationRangeFilters.filter(elem => {
        const extractedValue = extractor(elem)
        return (
          typeof extractedValue === "number" &&
          extractedValue >= filterValue.min &&
          extractedValue <= filterValue.max
        )
      })
    }
  }
  for (const [population, filterValue] of Object.entries(saigeHeritabilityFilters)) {
    if (filterValue !== undefined) {
      const extractor = getMetricValueForPopulationExtractor("saige_heritability_", population as PopulationCode)
      filteredByPopulationRangeFilters = filteredByPopulationRangeFilters.filter(elem => {
        const extractedValue = extractor(elem)
        return (
          typeof extractedValue === "number" &&
          extractedValue >= filterValue.min &&
          extractedValue <= filterValue.max
        )
      })
    }
  }
  for (const [population, filterValue] of Object.entries(lambdaGcFilters)) {
    if (filterValue !== undefined) {
      const extractor = getMetricValueForPopulationExtractor("lambda_gc_", population as PopulationCode)
      filteredByPopulationRangeFilters = filteredByPopulationRangeFilters.filter(elem => {
        const extractedValue = extractor(elem)
        return (
          typeof extractedValue === "number" &&
          extractedValue >= filterValue.min &&
          extractedValue <= filterValue.max
        )
      })
    }
  }
  const perPopulationLambdaGcExtremums = Object.fromEntries(
    Object.entries(perPopulationLambdaGcValues).map(([population, values]) => {
      const rawMin = min(values)!
      const rawMax = max(values)!
      const displayedMin = (rawMin < clampedLambdaGcMin) ? clampedLambdaGcMin : rawMin
      const displayedMax = (rawMax > clampedLambdaGcMax) ? clampedLambdaGcMax : rawMax
      return [population, {min: displayedMin, max: displayedMax}]
    })
  ) as PerPopulationExtremums
  const rawGlobalLambdaGcMin = min(allLambdaGcValues)!
  const rawGlobalLambdaGcMax = max(allLambdaGcValues)!
  return {
    globalNCasesUpperThreshold,
    perPopulationNCasesExtremums,
    globalNControlsUpperThreshold,
    perPopulationNControlsExtremums,
    perPopulationSaigeHeritabilityExtremums,

    globalLambdaGc: {
      min: (rawGlobalLambdaGcMin < clampedLambdaGcMin) ? clampedLambdaGcMin : rawGlobalLambdaGcMin,
      max: (rawGlobalLambdaGcMax > clampedLambdaGcMax) ? clampedLambdaGcMax : rawGlobalLambdaGcMax,
    },
    perPopulationLambdaGcExtremums,
    filteredByPopulationRangeFilters

  }

}
