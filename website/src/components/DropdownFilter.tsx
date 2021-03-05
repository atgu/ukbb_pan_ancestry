import { MenuItem, TextField } from "@material-ui/core";
import React from "react"

export interface Option {
  value: string, label: string
}
interface Props {
  options: Option[]
  setFilterValue: (value: string | undefined) => void
  filterValue: string | undefined
  label: string
  className?: string
}
export const DropdownFilter = (props: Props) => {
  const {options, setFilterValue, filterValue, label, className} = props;

  const optionElems = options.map((option, index) => (
    <MenuItem key={index} value={option.value}>{option.label}</MenuItem>
  ))
  const effectiveValue = (filterValue === undefined) ? "" : filterValue
  return (
      <TextField
        select={true}
        fullWidth={true}
        label={label}
        value={effectiveValue}
        SelectProps={{
          displayEmpty: true,
          renderValue: value => {return (value === "") ? value : options.find(option => option.value === value)!.label},
        }}
        onChange={e => {
          const newValue = (e.target.value === "") ? undefined : e.target.value
          setFilterValue(newValue)
        }}
        className={className}
      >
        <MenuItem value="">All</MenuItem>
        {optionElems}
      </TextField>
  )
}

