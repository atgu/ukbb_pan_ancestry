import { Box, Button, FormControl, FormControlLabel, FormGroup, FormLabel, makeStyles, Slider, Switch, TextField } from "@material-ui/core"
import React, { ChangeEvent, useCallback, useEffect, useMemo, useState } from "react"
import min from "lodash/min"
import max from "lodash/max"
import {  useAsyncDebounce } from "react-table";

const useStyles = makeStyles(() => ({

}))


export const NumberRangeColumnFilter = ({column}) => {
  const {filterValue, setFilter, preFilteredRows, id, } = column;
  const extremums = useMemo(() => {
    const allValues = preFilteredRows.map(row => row.values[id])
    return [ min(allValues), max(allValues) ] as const
  }, [preFilteredRows, id])
  const [minVal, setMinVal] = useState(() => {
    if (filterValue === undefined) {
      return extremums[0]
    } else {
      return filterValue[0]
    }
  })
  const [maxVal, setMaxVal] = useState(() => {
    if (filterValue === undefined) {
      return extremums[1]
    } else {
      return filterValue[1]
    }
  })
  const debouncedSetFilterValue = useAsyncDebounce((filterValue) => {
    setFilter(filterValue)
  }, 250)
  const handleMinChange = useCallback( (event: ChangeEvent<HTMLInputElement>) => {
    const value = Number(event.target.value)
    const newFilter = [value, filterValue[1]]
    setMinVal(value)
    debouncedSetFilterValue(newFilter)
  }, [setFilter, filterValue])
  const handleMaxChange = useCallback( (event: ChangeEvent<HTMLInputElement>) => {
    const value = Number(event.target.value)
    const newFilter = [filterValue[0], value]
    setMaxVal(value)
    debouncedSetFilterValue(newFilter)
  }, [setFilter, filterValue])
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
      setFilter(extremums)
    } else {
      setFilter(undefined)
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
          label={column.Header}
        />
        {content}
      </FormControl>
    </Box>
  )
}
