import React from "react"
import { ColumnInstance } from "react-table"
import { Datum } from "./types"

interface Props {
  columns: ColumnInstance<Datum>[]
}
export const ColumnGroupIndividualFilters = ({columns}: Props): React.ReactNode => {
  const filters = columns.map((col, index) => col.canFilter ? (
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
