import React from 'react'

export interface Option {
  value: string
  label: string
}
interface Props {
  options: Option[]
  setFilterValue: (value: string | undefined) => void
  filterValue: string | undefined
  label: string
  className?: string
}
export const DropdownFilter = (props: Props) => {
  const { options, setFilterValue, filterValue, label } = props

  const optionElems = options.map((option, index) => (
    <option key={index} value={option.value}>
      {option.label}
    </option>
  ))
  const effectiveValue = filterValue === undefined ? '' : filterValue
  return (
    <>
      <label>{label}</label>

      <select
        style={{ height: 30 }}
        name={label}
        id={label}
        value={effectiveValue}
        onChange={(e) => {
          const newValue = e.target.value === '' ? undefined : e.target.value
          setFilterValue(newValue)
        }}
      >
        {optionElems}
      </select>
    </>
  )
}
