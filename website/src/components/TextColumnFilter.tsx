import { TextField } from "@material-ui/core"
import React from "react"


export const TextColumnFilter = ({column: {filterValue, preFilteredRows, setFilter}}) => {
  const count = preFilteredRows.length
  return (
    <TextField
      type="search"
      size="small"
      margin="dense"
      fullWidth={true}
      value={filterValue || ""}
      onChange={e => setFilter(e.target.value)}
      placeholder={`Search ${count} records`}
    />
  )
}
