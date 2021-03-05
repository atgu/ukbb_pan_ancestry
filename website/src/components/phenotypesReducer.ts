import { PopulationCode } from "./populations"

// Each key determines whether a particular population should be shown or not across all metrics:
export type PerPopulationMetricsVisibility = {
  [K in PopulationCode]: boolean
}

// Each key determines the value of the range filter for a particular population:
export type PerPopulationRangeFilter = {
  [K in PopulationCode]: RangeFilterValue
}

export type RangeFilterValue = undefined | {min: number, max: number}

// All metrics that allow range filtering:
export enum RangeFilterMetric {
  NCases = "n_cases",
  NControls = "n_controls",
  SaigeHeritability = "saige_heritability",
  LambdaGc = "lambda_gc",
}

export enum ColumnGroupName {
  Downloads = "downloads",
  Populations = "populations",
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

const withAllFiltersDisabled: PerPopulationRangeFilter = {
  [PopulationCode.AFR]: undefined,
  [PopulationCode.AMR]: undefined,
  [PopulationCode.CSA]: undefined,
  [PopulationCode.EAS]: undefined,
  [PopulationCode.EUR]: undefined,
  [PopulationCode.MID]: undefined,
}

interface State {
  columnGroupVisibilities: ColumnGroupVisibility
  perPopulationMetricsVisibilities: PerPopulationMetricsVisibility
  traitTypeFilterValue: undefined | string
  sexFilterValue: undefined | string
  descriptionFilterValue: undefined | string
  // Range filter settings:
  [RangeFilterMetric.NCases]: PerPopulationRangeFilter
  [RangeFilterMetric.NControls]: PerPopulationRangeFilter
  [RangeFilterMetric.SaigeHeritability]: PerPopulationRangeFilter
  [RangeFilterMetric.LambdaGc]: PerPopulationRangeFilter,
}

export const initialState: State = {
  traitTypeFilterValue: undefined,
  sexFilterValue: undefined,
  descriptionFilterValue: undefined,
  columnGroupVisibilities: {
    [ColumnGroupName.Downloads]: true,
    [ColumnGroupName.Populations]: false,
    [ColumnGroupName.NCases]: true,
    [ColumnGroupName.NControls]: true,
    [ColumnGroupName.SaigeHeritability]: true,
    [ColumnGroupName.LambdaGc]: true,
    [ColumnGroupName.Md5]: false,
  },
  perPopulationMetricsVisibilities: {
    [PopulationCode.AFR]: true,
    [PopulationCode.AMR]: true,
    [PopulationCode.CSA]: true,
    [PopulationCode.EAS]: true,
    [PopulationCode.EUR]: true,
    [PopulationCode.MID]: true,
  },
  [RangeFilterMetric.NCases]: {
    ...withAllFiltersDisabled
  },
  [RangeFilterMetric.NControls]: {
    ...withAllFiltersDisabled
  },
  [RangeFilterMetric.SaigeHeritability]: {
    ...withAllFiltersDisabled
  },
  [RangeFilterMetric.LambdaGc]: {
    ...withAllFiltersDisabled
  },
}

export enum ActionType {
  SET_COLUMN_GROUP_VISIBILITY,
  SET_POPULATION_METRICS_VISIBILITY,
  UPDATE_FILTER_ONE_POPULATION,
  DISABLE_FILTER_ONE_POPULATION,
  UPDATE_SEX_FILTER,
  UPDATE_TRAIT_TYPE_FILTER,
  UPDATE_DESCRIPTION_FILTER,
}

type Action = {
  type: ActionType.SET_COLUMN_GROUP_VISIBILITY,
  payload: {columnGroup: keyof ColumnGroupVisibility, isVisible: boolean}
} | {
  type: ActionType.SET_POPULATION_METRICS_VISIBILITY,
  payload: {population: PopulationCode, isVisible: boolean}
} | {
  type: ActionType.UPDATE_FILTER_ONE_POPULATION,
  payload: {population: PopulationCode, metric: RangeFilterMetric, min: number, max: number}
} | {
  type: ActionType.DISABLE_FILTER_ONE_POPULATION
  payload: {population: PopulationCode, metric: RangeFilterMetric}
} | {
  type: ActionType.UPDATE_SEX_FILTER,
  payload: {value: State["sexFilterValue"]}
} | {
  type: ActionType.UPDATE_TRAIT_TYPE_FILTER,
  payload: {value: State["traitTypeFilterValue"]}
} | {
  type: ActionType.UPDATE_DESCRIPTION_FILTER,
  payload: {value: State["descriptionFilterValue"]}
}

export const reducer = (prevState: State, action: Action): State => {
  let nextState: State
  if (action.type === ActionType.SET_COLUMN_GROUP_VISIBILITY) {
    const {columnGroup, isVisible} = action.payload;
    const hasColumnGroupVisibilityChanged = prevState.columnGroupVisibilities[columnGroup] !== isVisible
    const newValueForColumnGroup = (hasColumnGroupVisibilityChanged === true) ? {
      ...prevState.columnGroupVisibilities,
      [columnGroup]: isVisible,
    } : prevState.columnGroupVisibilities

    const tempNextState: State = {
      ...prevState,
      columnGroupVisibilities: newValueForColumnGroup
    }
    // If a column group has been turned off, make sure to disable all range filters that affect that column group:
    if (hasColumnGroupVisibilityChanged === true && columnGroup === ColumnGroupName.NCases && isVisible === false) {
      nextState = {
        ...tempNextState,
        [RangeFilterMetric.NCases]: {...withAllFiltersDisabled }
      }
    } else if (hasColumnGroupVisibilityChanged === true && columnGroup === ColumnGroupName.NControls && isVisible === false) {
      nextState = {
        ...tempNextState,
        [RangeFilterMetric.NControls]: {...withAllFiltersDisabled }
      }
    } else if (hasColumnGroupVisibilityChanged === true && columnGroup === ColumnGroupName.SaigeHeritability && isVisible === false) {
      nextState = {
        ...tempNextState,
        [RangeFilterMetric.SaigeHeritability]: {...withAllFiltersDisabled }
      }
    } else if (hasColumnGroupVisibilityChanged === true && columnGroup === ColumnGroupName.LambdaGc && isVisible === false) {
      nextState = {
        ...tempNextState,
        [RangeFilterMetric.LambdaGc]: {...withAllFiltersDisabled }
      }
    } else {
      nextState = tempNextState
    }
  } else if (action.type === ActionType.SET_POPULATION_METRICS_VISIBILITY) {
    const {population, isVisible} = action.payload;
    const newValueForPopulation = (prevState.perPopulationMetricsVisibilities[population] !== isVisible) ? {
      ...prevState.perPopulationMetricsVisibilities,
      [population]: isVisible,
    } : prevState.perPopulationMetricsVisibilities
    nextState = {
      ...prevState,
      perPopulationMetricsVisibilities: newValueForPopulation
    }
  } else if (action.type === ActionType.UPDATE_FILTER_ONE_POPULATION) {
    const {population, metric, min, max} = action.payload;
    const currentFilterValue = prevState[metric][population]
    const willFilterValueChange = (currentFilterValue === undefined) ||
      (currentFilterValue.min !== min) || (currentFilterValue.max !== currentFilterValue.max);
    const newValueForMetric = willFilterValueChange ? {
      ...prevState[metric],
      [population]: {min, max}
    } : prevState[metric]
    nextState = {
      ...prevState,
      [metric]: newValueForMetric
    }
  } else if (action.type === ActionType.DISABLE_FILTER_ONE_POPULATION) {
    const {population, metric} = action.payload;
    const newValueForMetric = (prevState[metric][population] !== undefined) ? {
      ...prevState[metric],
      [population]: undefined
    } : prevState[metric]
    nextState = {
      ...prevState,
      [metric]: newValueForMetric
    }
  } else if (action.type === ActionType.UPDATE_SEX_FILTER) {
    nextState = {
      ...prevState,
      sexFilterValue: action.payload.value,
    }
  } else if (action.type === ActionType.UPDATE_TRAIT_TYPE_FILTER) {
    nextState = {
      ...prevState,
      traitTypeFilterValue: action.payload.value,
    }
  } else if (action.type === ActionType.UPDATE_DESCRIPTION_FILTER) {
    nextState = {
      ...prevState,
      descriptionFilterValue: action.payload.value,
    }
  } else {
    nextState = prevState
  }
  return nextState
}
