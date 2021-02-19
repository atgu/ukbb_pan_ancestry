import { makeStyles, TextField } from "@material-ui/core"
import React, { Ref, useState } from "react"
import {  useAsyncDebounce } from "react-table";

const useStyles = makeStyles(() => ({
  labelUnfocused: {
    whiteSpace: "nowrap",
    overflow: "hidden",
    textOverflow: "ellipsis",
    maxWidth: "100%"
  },
  labelFocused: {
    whiteSpace: "unset",
    overflow: "unset",
    textOverflow: "unset",
    maxWidth: "unset",
  }
}))

interface Props {
  filterValue: string | undefined
  setFilterValue: (value: string | undefined) => void
  label: string
}

export const BaseTextFilter = React.forwardRef((props: Props, ref: Ref<HTMLInputElement>) => {
  const {filterValue, setFilterValue, label} = props;
  const classes = useStyles()
  const [value, setValue] = useState(filterValue)
  const onChange = useAsyncDebounce(val => {
    const newValue = (val === "") ? undefined : val
    setFilterValue(newValue)
  }, 250)

  return (
    <TextField
      type="search"
      fullWidth={true}
      value={value}
      onChange={e => {
        setValue(e.target.value)
        onChange(e.target.value)
      }}
      inputRef={ref}
      InputLabelProps={{
        classes: {
          root: classes.labelUnfocused,
          focused: classes.labelFocused,
        }
      }}
      label={label}
    />

  )
})
