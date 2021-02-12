import { Box, FormControl, FormControlLabel, Switch, TextField } from "@material-ui/core"
import React, { ChangeEvent, useCallback, useMemo, useState } from "react"
import {  useAsyncDebounce } from "react-table";
import { RangeFilterMetric, RangeFilterValue } from "./phenotypesReducer";
import { PopulationCode } from "./populations";

interface Props {
  label: string
  population: PopulationCode
  globalMinValue: number
  globalMaxValue: number
  metric: RangeFilterMetric
  filterValue: RangeFilterValue
  disableFilter: (args: {metric: RangeFilterMetric, population: PopulationCode}) => void
  updateFilter: (args: {metric: RangeFilterMetric, population: PopulationCode, min: number, max: number}) => void
}

export const NumberRangeColumnFilter = (props: Props) => {
  const {
    label, globalMinValue, globalMaxValue, metric,
    filterValue, disableFilter, updateFilter, population
  } = props;
  const [minVal, setMinVal] = useState(() => {
    return (filterValue === undefined) ? globalMinValue : filterValue.min
  })
  const [maxVal, setMaxVal] = useState(() => {
    return (filterValue === undefined) ? globalMaxValue : filterValue.max
  })
  const debouncedSetFilterValue = useAsyncDebounce((minAndMax: {min: number, max: number}) => {
    updateFilter({population, metric, min: minAndMax.min, max: minAndMax.max})
  }, 250)
  const handleMinChange = useCallback( (event: ChangeEvent<HTMLInputElement>) => {
    const value = Number(event.target.value)
    setMinVal(value)
    if (filterValue !== undefined) {
      // We can safely assume this because if `filterValue` is undefined, it
      // means the input fields will not be shown at all.
      debouncedSetFilterValue({min: value, max: filterValue.max})
    }
  }, [filterValue])
  const handleMaxChange = useCallback( (event: ChangeEvent<HTMLInputElement>) => {
    const value = Number(event.target.value)
    setMaxVal(value)
    if (filterValue !== undefined) {
      debouncedSetFilterValue({min: filterValue.min, max: value})
    }
  }, [filterValue])
  let content: React.ReactNode
  if (filterValue === undefined) {
    content = null
  } else {
    content = (
      <Box display="flex" flexDirection="row">
        <TextField
          value={minVal}
          type="number"
          defaultValue={minVal}
          InputLabelProps={{
            shrink: true,
          }}
          onChange={handleMinChange}
        />
        <div style={{padding: "10px"}}>to</div>
        <TextField
          type="number"
          value={maxVal}
          InputLabelProps={{
            shrink: true,
          }}
          onChange={handleMaxChange}
          style={{marginRight: "10px"}}
        />
      </Box>
    )
  }
  const handleSwitchChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    const isChecked = event.target.checked
    if (isChecked) {
      updateFilter({population, metric, min: globalMinValue, max: globalMaxValue})
    } else {
      disableFilter({population, metric})
    }
  }
  return (
    <Box display="flex" m={1}>
      <FormControl component="fieldset">
        <FormControlLabel
          control={
            <Switch
              checked={filterValue !== undefined}
              onChange={handleSwitchChange}
            />
          }
          label={label}
        />
        {content}
      </FormControl>
    </Box>
  )
}
