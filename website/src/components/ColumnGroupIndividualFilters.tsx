import React from "react"
import { ColumnInstance } from "react-table"
import { Datum } from "./types"

interface Props {
  columnGroup: ColumnInstance<Datum>
}
export const ColumnGroupIndividualFilters = ({columnGroup}: Props) => {
  const filters = columnGroup.columns.map((col, index) => col.canFilter ? (
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
