import { Button, Slider, TextField } from "@material-ui/core"
import React, { ChangeEvent, useCallback, useEffect, useMemo, useState } from "react"
import min from "lodash/min"
import max from "lodash/max"
import {  useAsyncDebounce } from "react-table";


export const NumberRangeColumnFilter = ({column: {filterValue, setFilter, preFilteredRows, id}}) => {
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
  if (filterValue === undefined) {
    return (
      <Button size="small" variant="outlined" onClick={() => setFilter(extremums)}>Enable Filter</Button>
    )
  } else {
    return (
      <div style={{display: "flex", alignItems: "center"}}>
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
        <Button size="small" variant="outlined" onClick={() => setFilter(undefined)}>Disable</Button>
      </div>
    )
  }
}
