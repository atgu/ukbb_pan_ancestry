import React, { useMemo } from "react"
import {FilterProps} from "react-table";
import {Datum} from "./types";

export const SelectColumnFilter = (props: FilterProps<Datum>) => {
  const {column, state} = props;
  const {filterValue, setFilter, preFilteredRows, id} = column;
  const options = useMemo(() => {
    const allOptions: string = preFilteredRows.map(row => row.values[id])
    const distinctOptions = new Set(allOptions)
    return [...distinctOptions.values()].sort()
  }, [id, preFilteredRows])
  const optionElems = options.map((option, index) => (
    <option key={index} value={option}>{option}</option>
  ))
  return (
    <select
      value={filterValue}
      onChange={e => setFilter(e.target.value || undefined)}
    >
      <option value="">All</option>
      {optionElems}
    </select>
  )
}

