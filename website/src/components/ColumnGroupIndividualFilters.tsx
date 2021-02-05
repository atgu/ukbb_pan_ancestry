import React from "react"
import { ColumnInstance } from "react-table"
import { Datum } from "./types"

interface Props {
  columnGroup: ColumnInstance<Datum>
  // Do not render filters for these column ids:
  excludedColumnIds: string[]
}
export const ColumnGroupIndividualFilters = ({columnGroup, excludedColumnIds}: Props) => {
  const filters = columnGroup.columns.filter(col => !excludedColumnIds.includes(col.id)).map((col, index) => col.canFilter ? (
    <div key={index}>
      {col.render("Filter")}
    </div>
  ) : null)
  return (
    <>
      {filters}
    </>
  )
}
