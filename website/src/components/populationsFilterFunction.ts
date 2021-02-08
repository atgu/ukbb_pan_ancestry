import identity from "lodash/identity"

export enum Operator {
  AND = "AND",
  OR = "OR"
}

export interface FilterValue {
  operator: Operator
  populations: string[]
}

export const populationsFilterFunction = (rows, columnIds, filterValue: FilterValue | undefined) => {
  if (filterValue === undefined) {
    return rows;
  } else {
    if (filterValue.populations.length === 0) {
      return rows
    } else {
      const id = columnIds[0]
      const predicates = filterValue.populations.map(population => {
        return (populationsInRow: string[]) => (populationsInRow.includes(population))
      })
      return rows.filter(row => {
        const populationsInRow = row.values[id].split(",")
        const predicateEvaluations = predicates.map(predicate => predicate(populationsInRow))
        return (filterValue.operator === Operator.AND) ? predicateEvaluations.every(identity) : predicateEvaluations.some(identity)
      })
    }
  }
}
