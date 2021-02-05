import { TextField } from "@material-ui/core"
import React from "react"


export const TextColumnFilter = ({column}) => {
  const {filterValue, preFilteredRows, setFilter} = column;
  const count = preFilteredRows.length
  return (
    <TextField
      type="search"
      label={column.Header}
      fullWidth={true}
      value={filterValue || ""}
      onChange={e => setFilter(e.target.value)}
      placeholder={`Search ${count} records`}
    />
  )
}
