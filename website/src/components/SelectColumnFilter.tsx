import { MenuItem, Select } from "@material-ui/core";
import React, { useEffect, useMemo } from "react"
import {FilterProps} from "react-table";
import {Datum} from "./types";

export const SelectColumnFilter = (props: FilterProps<Datum>) => {
  const {column} = props;
  const {filterValue, setFilter, preFilteredRows, id} = column;
  const options = useMemo(() => {
    const allOptions: string = preFilteredRows.filter(row => !!row.values[id]).map(row => row.values[id])
    const distinctOptions = new Set(allOptions)
    return [...distinctOptions.values()].sort()
  }, [id, preFilteredRows])
  const optionElems = options.map((option, index) => (
    <MenuItem key={index} value={option}>{option}</MenuItem>
  ))
  const effectiveValue = (filterValue === undefined) ? "" : filterValue
  return (
    <Select
      value={effectiveValue}
      displayEmpty={true}
      onChange={e => setFilter(e.target.value)}
    >
      <MenuItem value="">All</MenuItem>
      {optionElems}
    </Select>
  )
}

