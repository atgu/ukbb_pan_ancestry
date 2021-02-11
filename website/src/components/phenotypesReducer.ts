import { PopulationCode } from "./populations"

// Each key determines whether a particular population should be shown or not across all metrics:
export type PerPopulationMetricsVisibility = {
  [K in PopulationCode]: boolean
}

export enum ColumnGroupName {
  Analysis = "analysis",
  Downloads = "downloads",
  Description = "description",
  NCases = "nCases",
  NControls = "nControls",
  SaigeHeritability = "saigeHeritability",
  LambdaGc = "lambdaGc",
  Md5 = "md5",
}
// Each key determines whether a group of related columns should be shown:
export type ColumnGroupVisibility = {
  [K in ColumnGroupName]: boolean
}

// Each key determines the value of the range filter for a particular population:
type PerPopulationRangeFilter = {
  [K in PopulationCode]: undefined | [number, number]
}

interface State {
  columnGroupVisibilities: ColumnGroupVisibility
  perPopulationMetricsVisibilities: PerPopulationMetricsVisibility
  // nCasesFilters: PerPopulationRangeFilter
  // nControlsFilters: PerPopulationRangeFilter
  // saigeHeritabilityFilters: PerPopulationRangeFilter
  // lambdaGcFilters: PerPopulationRangeFilter
}

export const initialState: State = {
  columnGroupVisibilities: {
    [ColumnGroupName.Analysis]: true,
    [ColumnGroupName.Downloads]: false,
    [ColumnGroupName.Description]: true,
    [ColumnGroupName.NCases]: false,
    [ColumnGroupName.NControls]: false,
    [ColumnGroupName.SaigeHeritability]: false,
    [ColumnGroupName.LambdaGc]: false,
    [ColumnGroupName.Md5]: false,
  },
  perPopulationMetricsVisibilities: {
    [PopulationCode.AFR]: true,
    [PopulationCode.AMR]: true,
    [PopulationCode.CSA]: true,
    [PopulationCode.EAS]: true,
    [PopulationCode.EUR]: true,
    [PopulationCode.MID]: true,
  }
}

export enum ActionType {
  SET_COLUMN_GROUP_VISIBILITY,
  SET_POPULATION_METRICS_VISIBILITY,
}

type Action = {
  type: ActionType.SET_COLUMN_GROUP_VISIBILITY,
  payload: {columnGroup: keyof ColumnGroupVisibility, isVisible: boolean}
} | {
  type: ActionType.SET_POPULATION_METRICS_VISIBILITY,
  payload: {population: PopulationCode, isVisible: boolean}
}

export const reducer = (prevState: State, action: Action): State => {
  let nextState: State
  if (action.type === ActionType.SET_COLUMN_GROUP_VISIBILITY) {
    const {columnGroup, isVisible} = action.payload;
    const newValue = (prevState.columnGroupVisibilities[columnGroup] !== isVisible) ? {
      ...prevState.columnGroupVisibilities,
      [columnGroup]: isVisible,
    } : prevState.columnGroupVisibilities
    nextState = {
      ...prevState,
      columnGroupVisibilities: newValue
    }
  } else if (action.type === ActionType.SET_POPULATION_METRICS_VISIBILITY) {
    const {population, isVisible} = action.payload;
    const newValue = (prevState.perPopulationMetricsVisibilities[population] !== isVisible) ? {
      ...prevState.perPopulationMetricsVisibilities,
      [population]: isVisible,
    } : prevState.perPopulationMetricsVisibilities
    nextState = {
      ...prevState,
      perPopulationMetricsVisibilities: newValue
    }
  } else {
    nextState = prevState
  }
  return nextState
}
