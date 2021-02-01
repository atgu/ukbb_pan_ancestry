import {matchSorter} from "match-sorter"

export const fuzzyTextFilterFunction = (rows, id, filterValue) => {
  if (!!filterValue) {
    return matchSorter(rows, filterValue, {
      keys: [row => row.values[id]]
    })
  } else {
    return rows
  }
}

