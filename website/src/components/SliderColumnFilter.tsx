import { Slider } from "@material-ui/core"
import React, { useCallback, useEffect, useMemo } from "react"
import min from "lodash/min"
import max from "lodash/max"


export const SliderColumnFilter = ({column: {filterValue, setFilter, preFilteredRows, id}}) => {
  const extremums = useMemo(() => {
    const allValues = preFilteredRows.map(row => row.values[id])
    return [ min(allValues), max(allValues) ] as const
  }, [preFilteredRows, id])
  const handleChange = useCallback((_: unknown, newValue) => setFilter(newValue), [setFilter])
  useEffect(() => {
    if (filterValue === undefined) {
      setFilter(extremums)
    }
  })
  if (filterValue === undefined) {
    return null
  } else {
    const [minVal, maxVal] = extremums;
    return (
      <Slider
        type="range"
        min={minVal}
        max={maxVal}
        value={filterValue}
        onChange={handleChange}
      />
    )
  }
}
