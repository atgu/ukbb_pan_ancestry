export const rangeFilterFunction = (rows, columnIds, filterValue: undefined | [number, number]) => {
  if (filterValue === undefined) {
    return rows
  } else {
    const id = columnIds[0]
    const [minVal, maxVal] = filterValue
    const filteredRows = rows.filter(row => {
      const rowValue = row.values[id]
      return (typeof rowValue === "number") && rowValue >= minVal && rowValue <= maxVal
    })
    return filteredRows
  }
}

