export const customIncludesFilterFn = (rows, id, filterValue) => {
  if (!!filterValue) {
    return rows.filter(row =>  {
      const rowValue = row.values[id]
      return !!rowValue ? rowValue.includes(filterValue) : true
    })
  } else {
    return rows
  }
}
