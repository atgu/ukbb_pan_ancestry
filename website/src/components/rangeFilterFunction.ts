export const rangeFilterFunction = (rows, columnIds, filterValue: undefined | [number, number]) => {
  let result
  if (filterValue === undefined) {
    result = rows
  } else {
    const id = columnIds[0]
    const [minVal, maxVal] = filterValue
    const filteredRows = rows.filter(row => {
      const rowValue = row.values[id]
      return (typeof rowValue === "number") && rowValue >= minVal && rowValue <= maxVal
    })
    result = filteredRows
  }
  return result
}

