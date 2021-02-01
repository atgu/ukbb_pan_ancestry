import React from "react"


export const DefaultColumnFilter = ({column: {filterValue, preFilteredRows, setFilter}}) => {
  const count = preFilteredRows.length
  return (
    <input
      value={filterValue || ""}
      onChange={e => setFilter(e.target.value)}
      placeholder={`Search ${count} records`}
    />
  )
}
